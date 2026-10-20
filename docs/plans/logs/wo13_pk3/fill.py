"""PK3: Step-1's 2026-10-16 ORA families with fewer than live.TOP passing sets (read only).

For each, the sets holding at least one call that fill it to TOP by p, with -log10 p and
-log10 FDR, the two lengths a grey bar could take.
"""

import glob
import math

import polars as pl

from ctx.de.stats import ALPHA
from ctx.function.live import TOP

for path in sorted(glob.glob("/Users/rytis/ctx_tier1/outputs/06_function/2026-10-16_*_ora_v1.csv")):
    f = pl.read_csv(path, comment_prefix="#")
    for d in ("up", "down"):
        rows = f.filter(pl.col("direction") == d).sort("p", "gene_set")
        passing = rows.filter(pl.col("fdr") < ALPHA).height
        if passing >= TOP:
            continue
        fill = rows.filter((pl.col("fdr") >= ALPHA) & (pl.col("hits") > 0)).head(TOP - passing)
        shown = [
            (r["gene_set"][9:], r["hits"], f"{-math.log10(r['p']):.2f}",
             f"{-math.log10(r['fdr']):.3f}")
            for r in fill.iter_rows(named=True)
        ]  # fmt: skip
        name = path.rsplit("tier1_", 1)[1]
        print(name, d, "passing", passing, "fill (set, hits, -log10 p, -log10 fdr)", shown)
