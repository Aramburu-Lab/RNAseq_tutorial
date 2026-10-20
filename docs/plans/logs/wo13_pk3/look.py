"""PK3: passing ORA sets per direction in each of Step-1's 2026-10-16 ORA tables (read only)."""

import glob

import polars as pl

from ctx.de.stats import ALPHA

for path in sorted(glob.glob("/Users/rytis/ctx_tier1/outputs/06_function/2026-10-16_*ora*_v1.csv")):
    f = pl.read_csv(path, comment_prefix="#")
    if "direction" not in f.columns:
        continue
    counts = f.group_by("direction").agg(n=pl.len(), passing=(pl.col("fdr") < ALPHA).sum())
    print(path.rsplit("/", 1)[1], sorted(counts.rows()))
