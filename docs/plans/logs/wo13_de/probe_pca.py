"""Package wo13-de: M5's blind-VST PCA step by step, RSS and traced bytes after each (S.5).

tracemalloc sees numpy's buffers (numpy reports them to it) but not polars'; RSS sees both.
``low`` passes pydeseq2's own ``low_memory=True``.

Run: bash env.sh python docs/plans/logs/wo13_de/probe_pca.py <01_gather dir> <dataset> [low]
"""

import sys
import tracemalloc
import warnings
from pathlib import Path

import psutil

from ctx import runtime

runtime.start(None, Path.cwd())
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import polars as pl  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.default_inference import DefaultInference  # noqa: E402

from ctx.dataset import load  # noqa: E402
from ctx.de.run import _subset, newest  # noqa: E402
from ctx.gather import gdc  # noqa: E402

folder, dataset, low = Path(sys.argv[1]), sys.argv[2], "low" in sys.argv[3:]
me = psutil.Process()
tracemalloc.start()


def at(label: str) -> None:
    """Print RSS, traced current and traced peak since the last call."""
    current, peak = tracemalloc.get_traced_memory()
    print(f"{label:<34}{me.memory_info().rss:>16,}{current:>16,}{peak:>16,}", flush=True)
    tracemalloc.reset_peak()


print(f"{'step':<34}{'RSS':>16}{'traced now':>16}{'traced peak':>16}")
spec = load(Path(__file__).parents[4] / "conf" / "datasets" / f"{dataset}.toml")
sheet = pl.read_parquet(newest(folder, dataset, "samplesheet"))
sheet = gdc.chosen(sheet) if "chosen" in sheet.columns else sheet
counts = pl.read_parquet(newest(folder, dataset, "counts"))
print("counts estimated_size", f"{counts.estimated_size():,}")
samples = _subset(sheet, spec.de.pca[0].filters, spec.de.pca[0].one_per)["sample"].to_list()
at(f"read, {len(samples)} samples")
matrix = counts.select(samples).to_numpy()
at(f"to_numpy {matrix.dtype} {matrix.nbytes:,}")
nonzero = matrix.sum(axis=1) > 0
genes = counts["gene_id"].filter(pl.Series(nonzero)).to_list()
frame = pd.DataFrame(matrix[nonzero].T, index=samples, columns=genes)
at(f"pandas frame {int(nonzero.sum())} genes")
metadata = pd.DataFrame({"all": ["all"] * len(samples)}, index=samples)
dds = DeseqDataSet(counts=frame, metadata=metadata, design="~1",
                   inference=DefaultInference(n_cpus=4), quiet=True, low_memory=low)  # fmt: skip
at(f"DeseqDataSet X {dds.X.dtype} {dds.X.nbytes:,}")
dds.fit_size_factors(fit_type=dds.size_factors_fit_type)
at("fit_size_factors")
dds.obsm["design_matrix_buffer"] = dds.obsm["design_matrix"].copy()
dds.obsm["design_matrix"] = pd.DataFrame(1, index=dds.obs_names, columns=["Intercept"])
dds.vst_fit_type = dds.fit_type
dds.fit_genewise_dispersions(vst=True)
at("fit_genewise_dispersions")
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="The dispersion trend curve fitting did not")
    dds._fit_parametric_dispersion_trend(vst=True)
at(f"trend {dds.vst_fit_type}")
dds.obsm["design_matrix"] = dds.obsm["design_matrix_buffer"].copy()
del dds.obsm["design_matrix_buffer"]
dds.layers["vst_counts"] = dds.vst_transform()
at("vst_transform")
values = np.asarray(dds.layers["vst_counts"])
top = np.argsort(-values.var(axis=0), kind="stable")[:500]
at("variance ranking")
print("layers", {k: f"{v.nbytes:,}" for k, v in dds.layers.items()})
print("peak RSS (getrusage)", f"{runtime.measure(0.0).peak_rss_self:,}")
