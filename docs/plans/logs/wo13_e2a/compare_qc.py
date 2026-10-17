"""Two ctx qc output folders, table by table: equal bytes, else equal as sets of rows.

Prints the D-061 counts of each: read_length flags, runs flagged, flags by metric.
    python compare_qc.py <folder a> <folder b>
"""

import sys
from pathlib import Path

import polars as pl


def tables(folder: Path) -> dict[str, Path]:
    """Each CSV by its description: the name after the dataset slug."""
    return {p.name.split("_tier1_", 1)[1]: p for p in sorted(folder.glob("*_tier1_*.csv"))}


def counts(folder: Path) -> str:
    """The D-061 numbers of one run."""
    found = tables(folder)
    flags = pl.read_csv(found["flags_v1.csv"])
    status = pl.read_csv(found["status_v1.csv"])
    reads = (pl.read_csv(found["read_flags_v1.csv"]) if "read_flags_v1.csv" in found
             else flags.head(0))  # fmt: skip
    length = flags.filter(pl.col("metric") == "read_length").height + reads.filter(
        pl.col("metric") == "read_length").height  # fmt: skip
    flagged = status.filter(pl.col("qc_status") != "pass").height
    return (f"read_length flags {length}; runs flagged {flagged} of {status.height}; "
            f"sample flags {flags.height}")  # fmt: skip


a, b = (Path(p) for p in sys.argv[1:3])
left, right = tables(a), tables(b)
print(f"a {a}\n  {counts(a)}\nb {b}\n  {counts(b)}")
for name in sorted(left.keys() | right.keys()):
    if name not in left or name not in right:
        print(f"{name}: only in {'a' if name in left else 'b'}")
        continue
    if left[name].read_bytes() == right[name].read_bytes():
        print(f"{name}: bytes equal")
        continue
    x, y = pl.read_csv(left[name]), pl.read_csv(right[name])
    same = x.columns == y.columns and x.sort(x.columns).equals(y.sort(y.columns))
    print(f"{name}: bytes differ; rows equal as sets {same}")
