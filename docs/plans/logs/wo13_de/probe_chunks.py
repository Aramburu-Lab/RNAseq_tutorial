"""Package wo13-de: are pydeseq2's rough and moments dispersions the same to the bit in chunks?

Each gene's estimate reads only that gene's column, so a chunked call should equal the whole
one; BLAS may still round a multi-column least squares by its width. Compares both on one
dataset's PCA samples (design ~1) and, with ``paired``, on a 9-column one-hot design.

Run: bash env.sh python docs/plans/logs/wo13_de/probe_chunks.py <01_gather dir> <dataset>
"""

import sys
from pathlib import Path

from ctx import runtime

runtime.start(None, Path.cwd())
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import polars as pl  # noqa: E402
from pydeseq2 import utils  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402

from ctx.dataset import load  # noqa: E402
from ctx.de.deseq import samples_by_genes  # noqa: E402
from ctx.de.run import _subset, newest  # noqa: E402
from ctx.gather import gdc  # noqa: E402

CHUNK = 4096
folder, dataset = Path(sys.argv[1]), sys.argv[2]
spec = load(Path(__file__).parents[4] / "conf" / "datasets" / f"{dataset}.toml")
sheet = pl.read_parquet(newest(folder, dataset, "samplesheet"))
sheet = gdc.chosen(sheet) if "chosen" in sheet.columns else sheet
counts = pl.read_parquet(newest(folder, dataset, "counts"))
samples = _subset(sheet, spec.de.pca[0].filters, spec.de.pca[0].one_per)["sample"].to_list()
nonzero = counts.select(pl.any_horizontal(pl.col(samples) > 0)).to_series().to_numpy()
metadata = pd.DataFrame({"all": ["all"] * len(samples)}, index=samples)
dds = DeseqDataSet(counts=samples_by_genes(counts, samples, nonzero), metadata=metadata,
                   design="~1", quiet=True)  # fmt: skip
dds.fit_size_factors()
normed = dds.layers["normed_counts"][:, np.arange(dds.n_vars)]
sf = dds.obs["size_factors"]
designs = {"~1": dds.obsm["design_matrix"].values}
rng = np.random.default_rng(0)
designs["one-hot, 9 columns"] = np.eye(9)[rng.integers(0, 9, len(samples))]
for name, design in designs.items():
    whole = utils.fit_rough_dispersions(normed, design)
    parts = np.concatenate([utils.fit_rough_dispersions(normed[:, i : i + CHUNK], design)
                            for i in range(0, normed.shape[1], CHUNK)])  # fmt: skip
    print(f"rough {name}: {normed.shape}, equal to the bit {np.array_equal(whole, parts)}, "
          f"differing genes {int((whole != parts).sum())}")  # fmt: skip
whole = utils.fit_moments_dispersions(normed, sf)
parts = np.concatenate([utils.fit_moments_dispersions(normed[:, i : i + CHUNK], sf)
                        for i in range(0, normed.shape[1], CHUNK)])  # fmt: skip
print(f"moments: equal to the bit {np.array_equal(whole, parts)}")
