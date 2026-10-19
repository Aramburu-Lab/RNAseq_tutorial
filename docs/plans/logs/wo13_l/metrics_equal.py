"""Package wo13-l: every metrics table of the after roots against the before roots (main's code).

Usage: python metrics_equal.py
Prints, per table, the rows equal in value, and the metrics only one side holds.
"""

import math
from pathlib import Path

import polars as pl

ROOTS = Path(__file__).absolute().parents[4] / "outputs" / "scratch_roots"
for before in sorted(ROOTS.glob("*_before/outputs/*/*_metrics_v1.csv")):
    after = Path(str(before).replace("_before/", "_after/"))
    if not after.exists():
        print("no after table:", before.name)
        continue
    if "metric" not in pl.read_csv(before, n_rows=1).columns:
        print("not a metric table:", before.name)
        continue
    a = dict(pl.read_csv(before, infer_schema_length=0).select("metric", "value").rows())
    b = dict(pl.read_csv(after, infer_schema_length=0).select("metric", "value").rows())

    def same(x: str | None, y: str | None) -> bool:
        """Equal as numbers where both are, else as text."""
        try:
            fx, fy = float(x), float(y)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            return x == y
        return fx == fy or (math.isnan(fx) and math.isnan(fy))

    equal = sum(same(a[k], b[k]) for k in a.keys() & b.keys())
    print(
        f"{before.parent.parent.parent.name[:-7]} {before.name[11:]}: {equal} of {len(a)} equal;"
        f" only after {sorted(b.keys() - a.keys())}; only before {sorted(a.keys() - b.keys())}"
    )
