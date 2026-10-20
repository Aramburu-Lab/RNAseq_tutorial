"""PK3: the link a reader follows to Step-1's M6 with 6 up calls, and Python's bars for it.

The first 6 up calls of the featured contrast by padj (short_list.py's list), as a ``#v=1`` fragment
on the scratch root's M6 page, and ``live.bars`` of the same list, to set beside the page.
"""

import glob
import sys
from pathlib import Path

import polars as pl

from ctx.function import lite, live
from ctx.function.page import set_words
from ctx.report import genes

ROOT = "/Users/rytis/ctx_tier1"
page = Path(sys.argv[1])
(path,) = glob.glob(f"{ROOT}/outputs/05_de/2026-10-16_de_tier1_ndn_gcsf_results_v1.parquet")
results = pl.read_parquet(path)
tested = results.filter(pl.col("padj").is_not_null() & pl.col("gene_name").is_not_null())
six = tested.filter(pl.col("call") == "up").sort("padj").head(6)
table = genes.table(tested["gene_id"].to_list(), tested["gene_name"].to_list(),
                    tested["log2FoldChange"].to_list())  # fmt: skip
fragment = genes.encode(genes.select("tier1:ndn_gcsf", six["gene_id"].to_list(), table))
print(page.as_uri() + fragment)
listed = Path(sys.argv[2])
listed.write_text("\n".join(six["gene_name"]) + "\n")
sets = lite.read_gmt(Path(f"{ROOT}/resources/msigdb/h.all.v2026.1.Hs.symbols.gmt"))
for r in live.bars(live.listed(listed, results, sets)).iter_rows(named=True):
    print(r["series"], set_words(r["gene_set"]), r["hits"], f"p {r['p']:.3g} fdr {r['fdr']:.3g}")
