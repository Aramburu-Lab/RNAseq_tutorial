"""Read-only probe: shapes and bytes of the tables M7's figures would need, on a real root."""

import io
import sys
import time
from pathlib import Path

import polars as pl

from ctx import runtime
from ctx.dataset import load
from ctx.proteo import run

root, dataset_file = Path(sys.argv[1]).expanduser(), Path(sys.argv[2])
spec = load(dataset_file)
settings = runtime.Settings("laptop", root, 4, None, time.time())
options = run.Options(dataset_file)
known = run.classes(settings, spec)
if spec.proteo.table is not None:
    from ctx.proteo.study import across_table

    m = across_table(settings, spec, options)
else:
    m = {
        "across_samples": run.across,
        "cross_dataset": run.cross,
        "tumour_normal": run.tumour_normal,
    }[spec.proteo.design](settings, spec, options)


def size(frame: pl.DataFrame) -> tuple[int, int]:
    """Parquet and CSV bytes."""
    a, b = io.BytesIO(), io.BytesIO()
    frame.write_parquet(a)
    frame.write_csv(b)
    return len(a.getvalue()), len(b.getvalue())


def wide(genes, samples, values) -> pl.DataFrame:
    """Gene then one column per sample."""
    return pl.DataFrame({"gene": list(genes), **{s: values[:, j] for j, s in enumerate(samples)}})


print("design", spec.proteo.design, "genes", m.genes.height)
print(
    "protein layer",
    m.protein.values.shape,
    size(wide(m.protein.genes, m.protein.samples, m.protein.values)),
)
if "rna" in m.extra:
    n = m.extra["rna"].shape
    names = [f"s{j}" for j in range(n[1])]
    print("rna paired", n, size(wide(m.genes["gene"], names, m.extra["rna"])))
    print(
        "protein paired",
        m.extra["values"].shape,
        size(wide(m.genes["gene"], names, m.extra["values"])),
    )
if "worked" in m.extra:
    r, v, _ = m.extra["worked"]
    names = [f"s{j}" for j in range(r.shape[1])]
    print("worked rna", r.shape, size(wide(m.genes["gene"], names, r)))
if "protein_all" in m.extra:
    print("protein_all", m.extra["protein_all"][1].shape)
