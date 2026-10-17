"""Package Q: where `ctx qc`'s count stage peaks on tcga_lung, by RSS sampled every 20 ms.

Each step of ``checkpoint.sample_metrics`` and ``run.run``'s reads runs in turn; the script
prints, per step, the RSS before it, the highest RSS sampled during it, and the frames'
``estimated_size``. Usage: py.sh python docs/plans/logs/wo13_q/profile_qc.py [chosen].
"""

import sys
import threading
import time
from pathlib import Path

import polars as pl
import psutil

from ctx.de.deseq import pca
from ctx.gather import gdc
from ctx.qc import checkpoint
from ctx.qc.identity import sex_check

GATHER = Path.home() / "ctx_lung" / "outputs" / "01_gather"
DAY = "2026-10-16"
PEERS = ("project",)
process = psutil.Process()
peak = [0]


def sample() -> None:
    """Record the highest RSS seen, every 20 ms, for the life of the process."""
    while True:
        peak[0] = max(peak[0], process.memory_info().rss)
        time.sleep(0.02)


def step(name: str, call):  # noqa: ANN001, ANN201 - a proof script
    """Run one step and print its RSS before and its peak during, in MB."""
    before = process.memory_info().rss
    peak[0] = before
    out = call()
    print(f"{name:28s} before {before / 1e6:8.0f} MB  peak {peak[0] / 1e6:8.0f} MB")
    return out


threading.Thread(target=sample, daemon=True).start()
path = GATHER / f"{DAY}_gather_tcga_lung_samplesheet_v1.parquet"
sheet = step("read sheet", lambda: pl.read_parquet(path))
if sys.argv[1:] == ["chosen"]:
    sheet = gdc.chosen(sheet)
names = sheet.unique("sample", keep="first", maintain_order=True)["sample"].to_list()
counts_path = GATHER / f"{DAY}_gather_tcga_lung_counts_v1.parquet"
counts = step("read counts", lambda: pl.read_parquet(counts_path, columns=["gene_id", *names]))
print(f"counts estimated_size {counts.estimated_size() / 1e6:.0f} MB, shape {counts.shape}")
genes = pl.read_parquet(GATHER / f"{DAY}_gather_tcga_lung_genes_v1.parquet")
samples = sheet.unique("sample", keep="first", maintain_order=True)
step("sex_check", lambda: sex_check(sheet, counts, genes))
found = step("pca (ctx.de.deseq)", lambda: pca(counts, names, threads=4))
step("count_metrics", lambda: checkpoint.count_metrics(counts, None))
step("replicate_spearman", lambda: checkpoint.replicate_spearman(counts, samples, PEERS))
print(f"process max RSS {process.memory_info().rss / 1e6:.0f} MB now")
