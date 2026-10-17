"""Package Q: tracemalloc's peak inside ``ctx.de.deseq.pca`` on tcga_lung's chosen samples.

numpy registers its buffers with tracemalloc, so the traced peak and the lines that held it
at that moment are numpy's and pandas' allocations; polars' own allocator is not traced.
Usage: py.sh python docs/plans/logs/wo13_q/profile_pca.py.
"""

import tracemalloc
from pathlib import Path

import polars as pl

from ctx.de.deseq import pca
from ctx.gather import gdc

GATHER = Path.home() / "ctx_lung" / "outputs" / "01_gather"
DAY = "2026-10-16"
sheet = gdc.chosen(pl.read_parquet(GATHER / f"{DAY}_gather_tcga_lung_samplesheet_v1.parquet"))
names = sheet["sample"].to_list()
counts = pl.read_parquet(GATHER / f"{DAY}_gather_tcga_lung_counts_v1.parquet",
                         columns=["gene_id", *names])  # fmt: skip
tracemalloc.start(25)
pca(counts, names, threads=4)
current, top = tracemalloc.get_traced_memory()
print(f"traced peak {top / 1e6:.0f} MB, current {current / 1e6:.0f} MB")
snapshot = tracemalloc.take_snapshot()
for stat in snapshot.statistics("lineno")[:12]:
    print(stat)
