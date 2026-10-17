"""Package wo13-de: every parquet table of two M5 output folders, compared value by value.

A table in both folders must be equal to the bit (nulls and NaNs included); one in only one
folder is listed. Run: bash env.sh python compare_tables.py <05_de before> <05_de after>
"""

import sys
from pathlib import Path

import polars as pl

before, after = Path(sys.argv[1]), Path(sys.argv[2])
names = {p.name: p for p in before.glob("*.parquet")}
later = {p.name: p for p in after.glob("*.parquet")}
equal, differ = 0, []
for name in sorted(names.keys() & later.keys()):
    a, b = pl.read_parquet(names[name]), pl.read_parquet(later[name])
    if a.equals(b, null_equal=True):
        equal += 1
    else:
        differ.append(name)
print(f"tables in both {len(names.keys() & later.keys())}; equal {equal}; differ {len(differ)}")
for name in differ:
    print("  differs:", name)
for name in sorted(names.keys() - later.keys()):
    print("  only before:", name)
for name in sorted(later.keys() - names.keys()):
    print("  only after:", name)
# Every figure's numbers are its CSV; they are compared byte for byte.
csvs = {p.name: p for p in before.glob("*.csv")}
csvs_after = {p.name: p for p in after.glob("*.csv")}
both = sorted(csvs.keys() & csvs_after.keys())
same = [n for n in both if csvs[n].read_bytes() == csvs_after[n].read_bytes()]
print(f"csv in both {len(both)}; byte-equal {len(same)}")
for name in sorted(set(both) - set(same)):
    print("  csv differs:", name)
