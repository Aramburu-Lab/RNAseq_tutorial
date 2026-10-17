"""Package Q: RSS per step of ``ctx.de.deseq.pca``'s body, replayed on tcga_lung's chosen samples.

The body is replayed here line for line (deseq.py at 628ad3a), only to time each step; the
module itself is another package's. RSS is sampled every 20 ms. Usage:
py.sh python docs/plans/logs/wo13_q/profile_pca_steps.py.
"""

import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import psutil
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference

from ctx.gather import gdc

GATHER = Path.home() / "ctx_lung" / "outputs" / "01_gather"
DAY = "2026-10-16"
process = psutil.Process()
peak = [0]


def sample() -> None:
    """Record the highest RSS seen, every 20 ms."""
    while True:
        peak[0] = max(peak[0], process.memory_info().rss)
        time.sleep(0.02)


def step(name: str, call):  # noqa: ANN001, ANN201 - a proof script
    """Run one step and print its RSS before and its peak during, in MB."""
    before = process.memory_info().rss
    peak[0] = before
    out = call()
    print(f"{name:24s} before {before / 1e6:8.0f} MB  peak {peak[0] / 1e6:8.0f} MB")
    return out


threading.Thread(target=sample, daemon=True).start()
sheet = gdc.chosen(pl.read_parquet(GATHER / f"{DAY}_gather_tcga_lung_samplesheet_v1.parquet"))
names = sheet["sample"].to_list()
counts = pl.read_parquet(GATHER / f"{DAY}_gather_tcga_lung_counts_v1.parquet",
                         columns=["gene_id", *names])  # fmt: skip
matrix = step("to_numpy", lambda: counts.select(names).to_numpy())
print("matrix dtype", matrix.dtype, matrix.shape)
nonzero = matrix.sum(axis=1) > 0
genes = counts["gene_id"].filter(pl.Series(nonzero)).to_list()
frame = step("pandas frame", lambda: pd.DataFrame(matrix[nonzero].T, index=names, columns=genes))
metadata = pd.DataFrame({"all": ["all"] * len(names)}, index=names)
dds = step("DeseqDataSet", lambda: DeseqDataSet(
    counts=frame, metadata=metadata, design="~1",
    inference=DefaultInference(n_cpus=4), quiet=True))  # fmt: skip
print("X dtype", dds.X.dtype, dds.X.shape)
step("vst", lambda: dds.vst(use_design=False))
for key, layer in dds.layers.items():
    print(f"layer {key}: {np.asarray(layer).dtype} {np.asarray(layer).nbytes / 1e6:.0f} MB")
