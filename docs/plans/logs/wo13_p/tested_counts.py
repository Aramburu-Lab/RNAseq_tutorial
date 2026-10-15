"""Package P, GS.3: tested genes, and those under the loosest preset, in each results table.

Run: python tested_counts.py <results.parquet> [...]
"""

import sys

import polars as pl

LOOSEST = 0.1  # the loosest cutoff preset, below which main drew every gene

for name in sys.argv[1:]:
    r = pl.read_parquet(name)
    tested = r.filter(pl.col("padj").is_not_null())
    drawn_always = tested.filter(pl.col("pvalue") < LOOSEST).height
    print(
        f"{name.rsplit('/', 1)[-1]}: rows {r.height:,}, tested {tested.height:,}, "
        f"raw p < 0.1 {drawn_always:,}, other {tested.height - drawn_always:,}"
    )
