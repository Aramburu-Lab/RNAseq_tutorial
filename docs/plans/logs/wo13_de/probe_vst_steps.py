"""Package wo13-de: the traced peak of each pydeseq2 step inside M5's blind VST (S.5).

Runs ``ctx.de.deseq.pca`` on one dataset's PCA samples with tracemalloc on; each wrapped
pydeseq2 step prints the bytes traced on entry and its own peak above them.

Run: bash env.sh python docs/plans/logs/wo13_de/probe_vst_steps.py <01_gather dir> <dataset>
"""

import sys
import tracemalloc
from pathlib import Path

from ctx import runtime

runtime.start(None, Path.cwd())
import polars as pl  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.default_inference import DefaultInference  # noqa: E402

from ctx.dataset import load  # noqa: E402
from ctx.de.deseq import pca  # noqa: E402
from ctx.de.run import _subset, newest  # noqa: E402
from ctx.gather import gdc  # noqa: E402

depth = [0]


def traced(owner: type, name: str) -> None:
    """Wrap ``owner.name`` to print its entry bytes and its peak above them."""
    function = getattr(owner, name)

    def wrapper(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202 - a probe
        entry = tracemalloc.get_traced_memory()[0]
        tracemalloc.reset_peak()
        depth[0] += 1
        try:
            return function(*args, **kwargs)
        finally:
            depth[0] -= 1
            peak = tracemalloc.get_traced_memory()[1]
            print(f"{'  ' * depth[0]}{name:<28} entry {entry:>14,}  peak {peak:>14,}  "
                  f"above entry {peak - entry:>14,}", flush=True)  # fmt: skip

    setattr(owner, name, wrapper)


for step in ("fit_size_factors", "fit_genewise_dispersions", "_fit_MoM_dispersions",
             "_fit_parametric_dispersion_trend", "vst_transform"):  # fmt: skip
    traced(DeseqDataSet, step)
# fit_rough_dispersions and fit_moments_dispersions are static methods; a plain wrapper
# would bind them, so they stay inside _fit_MoM_dispersions' line.
for step in ("lin_reg_mu", "alpha_mle"):
    traced(DefaultInference, step)
folder, dataset = Path(sys.argv[1]), sys.argv[2]
spec = load(Path(__file__).parents[4] / "conf" / "datasets" / f"{dataset}.toml")
sheet = pl.read_parquet(newest(folder, dataset, "samplesheet"))
sheet = gdc.chosen(sheet) if "chosen" in sheet.columns else sheet
counts = pl.read_parquet(newest(folder, dataset, "counts"))
samples = _subset(sheet, spec.de.pca[0].filters, spec.de.pca[0].one_per)["sample"].to_list()
tracemalloc.start()
pca(counts, samples, threads=4)
print(f"pca: {len(samples)} samples; run peak {tracemalloc.get_traced_memory()[1]:,} B")
