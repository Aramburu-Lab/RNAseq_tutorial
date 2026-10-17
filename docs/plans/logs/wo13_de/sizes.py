"""Package wo13-de: the shape, dtype mix and estimated_size of each gather table M5 reads.

Run: bash env.sh python docs/plans/logs/wo13_de/sizes.py <01_gather dir> <dataset>
"""

import sys
from pathlib import Path

import polars as pl

from ctx.de.run import newest
from ctx.gather import gdc

folder, dataset = Path(sys.argv[1]), sys.argv[2]
for what in ("samplesheet", "counts", "genes"):
    frame = pl.read_parquet(newest(folder, dataset, what))
    if what == "samplesheet" and "chosen" in frame.columns:
        print("chosen", gdc.chosen(frame).height, "of", frame.height)
    dtypes = sorted({str(t) for t in frame.dtypes})
    print(what, frame.shape, dtypes, f"{frame.estimated_size():,} B")
