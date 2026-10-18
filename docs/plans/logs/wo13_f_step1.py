"""wo13-f: Step-1's M5 rerun in a scratch root against ~/ctx_tier1's 2026-10-16 run (read only).

Compares the metrics value by value and every CSV of the same name byte for byte, and lists
the figures that now carry a spec.
"""

import sys
from pathlib import Path

import polars as pl

new_dir, old_dir = Path(sys.argv[1]), Path.home() / "ctx_tier1" / "outputs" / "05_de"
day = "2026-10-16"
new = {p.name[10:]: p for p in new_dir.glob("*.csv")}
old = {p.name[10:]: p for p in old_dir.glob(f"{day}_*.csv")}
a = pl.read_csv(next(p for k, p in new.items() if k.endswith("_metrics_v1.csv")))
b = pl.read_csv(next(p for k, p in old.items() if k.endswith("_metrics_v1.csv")))
joined = a.join(b, on="metric", how="full", suffix="_old")
equal = joined.filter(
    (pl.col("value") == pl.col("value_old"))
    | (pl.col("value").is_nan() & pl.col("value_old").is_nan())
)
print(f"metrics: {a.height} new, {b.height} recorded, {equal.height} equal")
same = [k for k in sorted(set(new) & set(old)) if new[k].read_bytes() == old[k].read_bytes()]
differ = [k for k in sorted(set(new) & set(old)) if k not in same]
print(f"CSV files: {len(same)} byte-equal, {len(differ)} differ: {differ}")
print(f"only new: {sorted(set(new) - set(old))}; only recorded: {sorted(set(old) - set(new))}")
ids = [p.name for p in new_dir.glob("*_README.md") if "## Figure" in p.read_text()]
print(f"figures with a spec: {len(ids)}: {sorted(i[11:-10] for i in ids)}")
