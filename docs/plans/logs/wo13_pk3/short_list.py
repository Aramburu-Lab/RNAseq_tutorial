"""PK3: lite.ora of 6 genes on Step-1's featured contrast (read only).

What a short list's top sets look like before choosing how to draw them: the first 6 up calls
by padj, and the first 6 tested genes by padj whatever their call, signed by raw log2 FC.
"""

import glob
import math
from pathlib import Path

import polars as pl

from ctx.function import lite, live

ROOT = "/Users/rytis/ctx_tier1"
(path,) = glob.glob(f"{ROOT}/outputs/05_de/2026-10-16_de_tier1_ndn_gcsf_results_v1.parquet")
results = pl.read_parquet(path)
sets = lite.read_gmt(Path(f"{ROOT}/resources/msigdb/h.all.v2026.1.Hs.symbols.gmt"))
tested = live.symbols(results)
ranked = results.filter(pl.col("padj").is_not_null() & pl.col("gene_name").is_not_null())
ranked = ranked.sort("padj")
for label, rows in (
    ("up calls", ranked.filter(pl.col("call") == "up").head(6)),
    ("any", ranked.head(6)),
):
    up = rows.filter(pl.col("log2FoldChange") >= 0)["gene_name"].to_list()
    down = rows.filter(pl.col("log2FoldChange") < 0)["gene_name"].to_list()
    print(label, "up", up, "down", down)
    for d, hits in (("up", up), ("down", down)):
        for r in sorted(lite.ora(hits, tested, sets), key=lambda r: (r.p, r.name))[:8]:
            line = f"{d} {r.name} hits {r.hits}/{r.size} p {r.p:.3g} fdr {r.fdr:.3g}"
            print(f"  {line} -log10 fdr {-math.log10(r.fdr):.3g}")
