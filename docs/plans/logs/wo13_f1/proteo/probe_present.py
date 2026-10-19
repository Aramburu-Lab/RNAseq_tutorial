"""Read-only probe: bytes of the protein layer's presence as Boolean and as UInt8, both formats."""

import io
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

from ctx import runtime
from ctx.dataset import load
from ctx.proteo import run

root, dataset_file = Path(sys.argv[1]).expanduser(), Path(sys.argv[2])
spec = load(dataset_file)
settings = runtime.Settings("laptop", root, 4, None, time.time())
design = {
    "across_samples": run.across,
    "tumour_normal": run.tumour_normal,
    "cross_dataset": run.cross,
}
m = design[spec.proteo.design](settings, spec, run.Options(dataset_file))
present = ~np.isnan(m.protein.values)
for name, cast in (("bool", lambda c: c), ("uint8", lambda c: c.astype(np.uint8))):
    frame = pl.DataFrame(
        {
            "gene": list(m.protein.genes),
            **{s: cast(present[:, j]) for j, s in enumerate(m.protein.samples)},
        }
    )
    a, b = io.BytesIO(), io.BytesIO()
    frame.write_parquet(a)
    frame.write_csv(b)
    print(name, present.shape, "parquet", len(a.getvalue()), "csv", len(b.getvalue()))
print("missing share", 1 - present.mean())
