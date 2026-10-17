"""Package wo13-e2a: Step-1's recount per contrast, against the n each contrast declares.

python m1_contrasts.py <gather samplesheet parquet> <dataset file>
"""

import sys
from pathlib import Path

import polars as pl

from ctx.dataset import load
from ctx.gather import composition

sheet, spec = pl.read_parquet(sys.argv[1]), load(Path(sys.argv[2]))
rows = composition.runs(sheet, composition.colours(sheet, spec))
print("colours offered:", composition.colours(sheet, spec))
every = composition.recount(rows, "batch", "runs")
by_runs = every.unique("group", maintain_order=True)
by_individuals = composition.recount(rows, "batch", "individuals").unique(
    "group", maintain_order=True
)
mono = by_runs.filter(pl.col("group") == "Monocyte hsct")["total"][0]
mono_i = by_individuals.filter(pl.col("group") == "Monocyte hsct")["total"][0]
print(f"every run: Monocyte hsct {mono} runs from {mono_i} individuals")
for contrast, kept in composition.members(sheet, spec):
    part = rows[kept]
    runs = composition.recount(part, "batch", "runs").unique("group", maintain_order=True)
    people = composition.recount(part, "batch", "individuals").unique("group", maintain_order=True)
    print(
        f"{contrast.name}: declared n {contrast.expected_n}; runs "
        f"{dict(runs.select('group', 'total').rows())}; individuals "
        f"{dict(people.select('group', 'total').rows())}"
    )
