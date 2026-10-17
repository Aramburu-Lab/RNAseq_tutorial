"""Package wo13-de: the memory layout pydeseq2's counts get, built as M5 builds them and lean.

A reduction's rounding follows the layout, so a leaner build must hand pydeseq2 the same
layout to give the same numbers. Prints each build's X flags, whether X shares our buffer,
and whether the two X are equal.

Run: bash env.sh python docs/plans/logs/wo13_de/probe_layout.py <01_gather dir> <dataset>
"""

import sys
from pathlib import Path

from ctx import runtime

runtime.start(None, Path.cwd())
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import polars as pl  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402

from ctx.de.run import newest  # noqa: E402

folder, dataset = Path(sys.argv[1]), sys.argv[2]
counts = pl.read_parquet(newest(folder, dataset, "counts"))
samples = counts.columns[1:201]
metadata = pd.DataFrame({"all": ["all"] * len(samples)}, index=samples)

matrix = counts.select(samples).to_numpy()
flags = matrix.flags
print("polars to_numpy", matrix.dtype, "C", flags.c_contiguous, "F", flags.f_contiguous)
nonzero = matrix.sum(axis=1) > 0
genes = counts["gene_id"].filter(pl.Series(nonzero)).to_list()
old = pd.DataFrame(matrix[nonzero].T, index=samples, columns=genes)
a = DeseqDataSet(counts=old, metadata=metadata, design="~1", quiet=True)
print("current X", a.X.dtype, "C", a.X.flags.c_contiguous, "F", a.X.flags.f_contiguous)

lean = counts.select(samples).filter(pl.Series(nonzero)).to_numpy().astype(np.int64).T
new = pd.DataFrame(lean, index=samples, columns=genes, copy=False)
b = DeseqDataSet(counts=new, metadata=metadata, design="~1", quiet=True)
print("lean X", b.X.dtype, "C", b.X.flags.c_contiguous, "F", b.X.flags.f_contiguous,
      "shares our buffer", np.shares_memory(b.X, lean))  # fmt: skip
print("equal", np.array_equal(a.X, b.X))
