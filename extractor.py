"""
fasta_region_extractor.py — Extract sub-sequences from a multi-FASTA file
using a coordinate manifest (CSV/TSV).  See README.md for full documentation.
"""

from __future__ import annotations

import csv
import logging
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CoordEntry:
    """
    Immutable container for a single sequence's cut specification.

    Attributes
    ----------
    seq_id : str
        Sequence identifier as read from the coordinate file (version-stripped
        during lookup, not stored here).
    start : int or None
        1-based inclusive start position, or None when not required by the mode.
    end : int or None
        1-based inclusive end position, or None when not required by the mode.
    mode : str
        Slicing mode: 'range', 'after', or 'before'.
    """

    seq_id: str
    start: Optional[int]
    end: Optional[int]
    mode: str


# ---------------------------------------------------------------------------
# Coordinate file parser
# ---------------------------------------------------------------------------

# Bytes read by csv.Sniffer to detect the delimiter.  1 MiB is sufficient for
# even wide files and avoids the misidentification risk of the default 1 KiB.
_SNIFFER_READ_SIZE = 1_048_576

# Fields that indicate a header row rather than a data row.
_KNOWN_HEADER_TOKENS = {"id", "seq_id", "sequence_id", "name"}


def _strip_version(seq_id: str) -> str:
    """
    Remove a trailing version suffix (e.g. '.1', '.2') from a sequence ID.

    This allows coordinate file entries recorded without a version suffix to
    match FASTA headers that carry one, and vice-versa.

    Parameters
    ----------
    seq_id : str
        Raw sequence identifier.

    Returns
    -------
    str
        Identifier with the last dot-delimited numeric token removed, or the
        original string if no such suffix exists.

    Examples
    --------
    >>> _strip_version("XP_015630497.1")
    'XP_015630497'
    >>> _strip_version("AT1G01010")
    'AT1G01010'
    """
    parts = seq_id.rsplit(".", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return seq_id


def parse_coords_file(file_path: str) -> Dict[str, CoordEntry]:
    """
    Parse a CSV or TSV coordinate manifest into a mapping keyed by sequence ID.

    The delimiter is detected automatically via :class:`csv.Sniffer`.  If
    detection fails the parser falls back to comma.

    Duplicate sequence IDs in the file raise a :class:`ValueError` rather than
    silently overwriting earlier entries, because two conflicting cut
    instructions for the same ID almost always indicate a mistake in the
    manifest.

    Parameters
    ----------
    file_path : str
        Path to the coordinate file.

    Returns
    -------
    dict
        Mapping of ``seq_id`` → :class:`CoordEntry`.  Both the original ID and
        its version-stripped form are stored as keys so that lookup succeeds
        regardless of which form appears in the FASTA.

    Raises
    ------
    FileNotFoundError
        If *file_path* does not exist.
    ValueError
        If a duplicate sequence ID is encountered.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Coordinate file not found: {path}")

    coords_map: Dict[str, CoordEntry] = {}
    skipped = 0

    with open(path, "r", newline="") as fh:
        sample = fh.read(_SNIFFER_READ_SIZE)
        fh.seek(0)

        try:
            dialect = csv.Sniffer().sniff(sample)
        except csv.Error:
            logger.warning(
                "Could not detect delimiter in coordinate file; defaulting to comma."
            )
            dialect = csv.excel  # comma-delimited fallback

        reader = csv.reader(fh, dialect)

        for line_num, row in enumerate(reader, start=1):
            # Skip blank lines and comment lines.
            if not row or row[0].strip().startswith("#"):
                continue

            # Detect and skip a header row.  A header row is identified by its
            # first field being a known header token.  Only line 1 is checked
            # because a header token appearing mid-file is more likely a data
            # error than an intentional header.
            if line_num == 1 and row[0].strip().lower() in _KNOWN_HEADER_TOKENS:
                logger.info("Header row detected in coordinate file — skipping.")
                continue

            try:
                seq_id = row[0].strip()
                s_val = row[1].strip()
                e_val = row[2].strip()
                mode = row[3].strip().lower()

                start: Optional[int] = int(s_val) if s_val and s_val != "0" else None
                end: Optional[int] = int(e_val) if e_val and e_val != "0" else None

            except (ValueError, IndexError):
                logger.warning(
                    f"Line {line_num}: malformed row {row!r} — skipping."
                )
                skipped += 1
                continue

            entry = CoordEntry(seq_id=seq_id, start=start, end=end, mode=mode)

            # Register under both the original ID and the version-stripped form
            # so that lookup works regardless of whether the FASTA uses versions.
            keys_to_register = {seq_id, _strip_version(seq_id)}
            for key in keys_to_register:
                if key in coords_map:
                    raise ValueError(
                        f"Duplicate sequence ID '{key}' detected in coordinate "
                        f"file (line {line_num}).  Each ID must appear only once."
                    )
                coords_map[key] = entry

    logger.info(
        f"Loaded coordinates for {len(coords_map)} sequence ID(s) "
        f"({skipped} malformed line(s) skipped)."
    )
    return coords_map


# ---------------------------------------------------------------------------
# Sequence slicing
# ---------------------------------------------------------------------------

class SequenceCutter:
    """
    Applies coordinate-based slicing to Biopython :class:`~Bio.SeqRecord.SeqRecord`
    objects.

    All public methods are static; this class acts as a namespace rather than
    a stateful object.
    """

    # Supported slicing modes.
    VALID_MODES = {"range", "after", "before"}

    @staticmethod
    def slice_record(
        record: SeqRecord,
        entry: CoordEntry,
    ) -> Optional[SeqRecord]:
        """
        Extract a sub-sequence from *record* according to *entry*.

        Coordinate validation is performed before slicing.  Errors produce a
        logged warning and return ``None`` rather than raising, so that a single
        bad entry does not abort processing of the entire file.

        Parameters
        ----------
        record : SeqRecord
            Source sequence record.
        entry : CoordEntry
            Cut specification for this sequence.

        Returns
        -------
        SeqRecord or None
            A new record whose ID is ``<original_id>:<start>-<end>`` and whose
            description is empty (preventing Biopython from repeating the
            original description in the FASTA header).  Returns ``None`` when
            the entry is invalid or the resulting slice is empty.

        Notes
        -----
        Slicing uses 1-based inclusive coordinates converted to Python's
        0-based half-open notation:

        - ``range`` : ``seq[start-1 : end]``
        - ``after``  : ``seq[start-1 :]``
        - ``before`` : ``seq[: end]``   — ``[:end]`` is equivalent to positions
          1 through *end* inclusive, matching biological convention without an
          adjustment because the upper bound is already exclusive in Python.
        """
        seq_len = len(record.seq)
        mode = entry.mode

        # --- Mode validation ------------------------------------------------
        if mode not in SequenceCutter.VALID_MODES:
            logger.warning(
                f"{record.id}: unknown mode '{mode}'. "
                f"Valid modes are: {', '.join(sorted(SequenceCutter.VALID_MODES))}. "
                f"Skipping."
            )
            return None

        # --- Coordinate requirements per mode --------------------------------
        if mode == "range" and (entry.start is None or entry.end is None):
            logger.warning(
                f"{record.id}: mode 'range' requires both start and end. Skipping."
            )
            return None

        if mode == "after" and entry.start is None:
            logger.warning(
                f"{record.id}: mode 'after' requires a start coordinate. Skipping."
            )
            return None

        if mode == "before" and entry.end is None:
            logger.warning(
                f"{record.id}: mode 'before' requires an end coordinate. Skipping."
            )
            return None

        # --- range-mode specific validation ----------------------------------
        if mode == "range":
            assert entry.start is not None and entry.end is not None  # type narrowing
            if entry.start > entry.end:
                logger.warning(
                    f"{record.id}: start ({entry.start}) is greater than "
                    f"end ({entry.end}) in 'range' mode — coordinates may be "
                    f"inverted. Skipping."
                )
                return None

        # --- Bounds-against-sequence-length check ----------------------------
        start = entry.start
        end = entry.end

        if start is not None and start > seq_len:
            logger.warning(
                f"{record.id}: start coordinate ({start}) exceeds sequence "
                f"length ({seq_len}). Skipping."
            )
            return None

        if end is not None and end > seq_len:
            logger.warning(
                f"{record.id}: end coordinate ({end}) exceeds sequence "
                f"length ({seq_len}); will be clamped to {seq_len}."
            )
            end = seq_len

        # --- Slicing ---------------------------------------------------------
        if mode == "range":
            assert start is not None and end is not None
            new_seq = record.seq[start - 1 : end]
            coord_label = f"{start}-{end}"

        elif mode == "after":
            assert start is not None
            new_seq = record.seq[start - 1 :]
            coord_label = f"{start}-{seq_len}"

        else:  # mode == "before"
            assert end is not None
            # [:end] in Python returns positions 1 through end (1-based inclusive),
            # which is correct without a -1 adjustment because Python's upper
            # bound is already exclusive.
            new_seq = record.seq[:end]
            coord_label = f"1-{end}"

        # --- Empty-slice guard -----------------------------------------------
        if len(new_seq) == 0:
            logger.warning(
                f"{record.id}: resulting slice is empty (coords: {coord_label}, "
                f"seq_len: {seq_len}). Skipping."
            )
            return None

        # --- Build output record ---------------------------------------------
        # description="" prevents Biopython from appending the original
        # description text after the new ID in the FASTA header line.
        return SeqRecord(
            new_seq,
            id=f"{record.id}:{coord_label}",
            description="",
        )


# ---------------------------------------------------------------------------
# FASTA processing
# ---------------------------------------------------------------------------

def process_fasta(
    input_path: Path,
    coords_map: Dict[str, CoordEntry],
) -> Tuple[List[SeqRecord], int, int, int]:
    """
    Iterate over a FASTA file and apply coordinate-based slicing to matched
    sequences.

    Sequences whose IDs are absent from *coords_map* are silently skipped
    (they are not part of the requested extraction).  Duplicate IDs within the
    FASTA are detected and trigger a warning for every occurrence beyond the
    first; all copies are still processed.

    Parameters
    ----------
    input_path : Path
        Path to the input FASTA file.
    coords_map : dict
        Mapping returned by :func:`parse_coords_file`.

    Returns
    -------
    tuple
        ``(records, n_matched, n_written, n_failed)`` where:

        - ``records``   – list of output :class:`~Bio.SeqRecord.SeqRecord` objects
        - ``n_matched`` – number of FASTA sequences matched to a coord entry
        - ``n_written`` – number of sequences successfully sliced and retained
        - ``n_failed``  – number of matched sequences that failed slicing
    """
    records: List[SeqRecord] = []
    seen_fasta_ids: Dict[str, int] = {}  # id → first-occurrence line number
    n_matched = 0
    n_failed = 0

    for record in SeqIO.parse(input_path, "fasta"):
        clean_id = record.id.strip()
        stripped_id = _strip_version(clean_id)

        # Duplicate-ID detection within the FASTA.
        for candidate in (clean_id, stripped_id):
            if candidate in seen_fasta_ids:
                logger.warning(
                    f"Duplicate FASTA ID '{clean_id}' encountered. "
                    f"All copies will be processed, but downstream tools may "
                    f"behave unexpectedly."
                )
                break
        seen_fasta_ids[clean_id] = 1

        # Look up by original ID first, then version-stripped.
        entry = coords_map.get(clean_id) or coords_map.get(stripped_id)
        if entry is None:
            continue

        n_matched += 1
        result = SequenceCutter.slice_record(record, entry)

        if result is not None:
            records.append(result)
        else:
            n_failed += 1

    n_written = len(records)
    return records, n_matched, n_written, n_failed


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """
    Construct and return the command-line argument parser.

    Returns
    -------
    argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog="fasta_region_extractor",
        description=(
            "Extract specific regions from sequences in a multi-FASTA file "
            "using a CSV/TSV coordinate manifest."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Coordinate file format (no header row):\n"
            "  ID, Start, End, Mode\n\n"
            "Modes:\n"
            "  range  – extract positions start..end (both required)\n"
            "  after  – extract from start to end of sequence\n"
            "  before – extract from position 1 to end\n\n"
            "Use 0 or leave empty for unused coordinates.\n"
            "Lines starting with '#' are treated as comments.\n\n"
            "Example:\n"
            "  XP_015630497.1,21,0,after\n"
            "  AT1G01010,1,300,range\n"
        ),
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        metavar="FASTA",
        help="Input multi-FASTA file.",
    )
    parser.add_argument(
        "-c", "--coords",
        required=True,
        metavar="CSV/TSV",
        help="Coordinate manifest (CSV or TSV; format: ID, Start, End, Mode).",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        metavar="FASTA",
        help="Output FASTA file for extracted sub-sequences.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """
    Entry point for the command-line interface.

    Parameters
    ----------
    argv : list of str, optional
        Argument list.  Defaults to ``sys.argv[1:]`` when ``None``.

    Returns
    -------
    int
        Exit code: 0 on success, 1 on error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    input_path = Path(args.input)
    output_path = Path(args.output)

    # --- Validate input paths ------------------------------------------------
    if not input_path.exists():
        logger.error(f"Input FASTA not found: {input_path}")
        return 1

    # --- Load coordinates ----------------------------------------------------
    try:
        coords_map = parse_coords_file(args.coords)
    except FileNotFoundError as exc:
        logger.error(str(exc))
        return 1
    except ValueError as exc:
        logger.error(str(exc))
        return 1

    if not coords_map:
        logger.error("No valid coordinate entries loaded. Aborting.")
        return 1

    # --- Process FASTA -------------------------------------------------------
    logger.info(f"Processing: {input_path}")
    records, n_matched, n_written, n_failed = process_fasta(input_path, coords_map)

    # --- Report unmatched coord entries --------------------------------------
    # Any ID in the coord file that was never seen in the FASTA is reported so
    # the user can detect typos or mismatched file pairs.
    matched_ids = {r.id.split(":")[0] for r in records}
    coord_ids = {entry.seq_id for entry in coords_map.values()}
    unmatched = coord_ids - {_strip_version(i) for i in matched_ids} - matched_ids
    if unmatched:
        logger.warning(
            f"{len(unmatched)} coordinate ID(s) had no matching FASTA sequence: "
            + ", ".join(sorted(unmatched))
        )

    # --- Write output --------------------------------------------------------
    if records:
        SeqIO.write(records, output_path, "fasta")
        logger.info(
            f"Done.  Matched: {n_matched} | Written: {n_written} | "
            f"Failed: {n_failed} | Output: {output_path}"
        )
        return 0
    else:
        logger.warning(
            "No sequences written.  Check that sequence IDs in the FASTA match "
            "those in the coordinate file."
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
