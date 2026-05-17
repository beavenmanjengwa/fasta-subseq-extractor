# fasta_region_extractor

Extract specific sub-sequences from a multi-FASTA file using a coordinate
manifest supplied as CSV or TSV.

Typical use cases include signal-peptide removal, domain extraction, and
N- or C-terminal trimming prior to structural or functional analysis.

---

## Requirements

- Python >= 3.8
- Biopython >= 1.79

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Usage

```bash
python fasta_region_extractor.py \
    -i proteins.fasta \
    -c coords.csv \
    -o trimmed.fasta
```

| Flag | Description |
|---|---|
| `-i` / `--input` | Input multi-FASTA file |
| `-c` / `--coords` | Coordinate manifest (CSV or TSV) |
| `-o` / `--output` | Output FASTA file |

---

## Coordinate file format

The coordinate file has **no header row**.  Each data row specifies one
sequence:

```
ID, Start, End, Mode
```

| Column | Description |
|---|---|
| `ID` | Sequence identifier matching the FASTA header (first whitespace-delimited token). Version suffixes (`.1`, `.2`) are matched automatically. |
| `Start` | 1-based inclusive start position. Use `0` or leave empty when the mode does not require it. |
| `End` | 1-based inclusive end position. Use `0` or leave empty when the mode does not require it. |
| `Mode` | `range`, `after`, or `before` (see below). |

Lines beginning with `#` are treated as comments and ignored.  A detected
header row (first field matching `id`, `seq_id`, `sequence_id`, or `name`) is
skipped automatically with a log notice.

Both CSV (comma-separated) and TSV (tab-separated) are accepted; the delimiter
is detected automatically.

### Modes

| Mode | Extracts | Requires |
|---|---|---|
| `range` | Positions `start` through `end` (inclusive) | Both `start` and `end` |
| `after` | From `start` to the end of the sequence | `start` only |
| `before` | From position 1 to `end` (inclusive) | `end` only |

All coordinates are **1-based and inclusive**, matching standard biological
convention.

### Example coordinate file

```
# Signal peptide removal (keep from position 21 onward)
XP_015630497.1,21,0,after

# Domain extraction
AT1G01010,50,300,range

# N-terminal fragment
Osmotin_WsP2,0,180,before
```

---

## Output

Each output FASTA record is assigned the ID `<original_id>:<start>-<end>`.
The original description is dropped to keep headers clean.

```
>XP_015630497.1:21-512
MKTIIALSYIFCLVFA...
```

---

## Validation and error handling

The following conditions are detected and reported:

| Condition | Behaviour |
|---|---|
| Inverted range (`start > end`) | Warning; sequence skipped |
| `start` or `end` beyond sequence length | Warning; `end` clamped or sequence skipped |
| Empty resulting slice | Warning; sequence skipped |
| Duplicate ID in coordinate file | Error; run aborted |
| Duplicate ID in FASTA | Warning; all copies processed |
| Coordinate ID with no matching FASTA sequence | Warning; listed by name |
| Malformed coordinate file row | Warning; row skipped |

---

## Slicing semantics

Internally, 1-based inclusive coordinates are converted to Python's 0-based
half-open slice notation:

| Mode | Python slice |
|---|---|
| `range` | `seq[start-1 : end]` |
| `after` | `seq[start-1 :]` |
| `before` | `seq[: end]` |

`before` mode does not apply a `-1` adjustment to `end` because Python's upper
bound is already exclusive, so `seq[:50]` returns positions 1–50 (1-based
inclusive) without further correction.

---

## Author

Beaven Manjengwa
