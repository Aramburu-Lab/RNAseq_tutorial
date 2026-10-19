"""wo13-f1 M10: is analyses.cramers_v the same on the same table every call?

Usage: wo13_f1_clinical_py.sh python <this file> <patients.parquet>
Prints, per pair of the confounding map's columns, how many distinct values 50 calls return.
"""

import sys

import polars as pl

from ctx.clinical import analyses

patients = pl.read_parquet(sys.argv[1])
table = patients.with_columns(analyses.age_group(pl.col("age")).alias("age_group"),
                              pl.col("stage").cast(pl.String))  # fmt: skip
names = [c for c in ("er", "stage", "histology", "menopause", "age_group", "plate", "site")
         if c in table.columns]  # fmt: skip
for i, a in enumerate(names):
    for b in names[i + 1 :]:
        seen = {analyses.cramers_v(table, a, b)[0] for _ in range(50)}
        print(f"{a:10} {b:10} distinct values in 50 calls: {len(seen)}")
