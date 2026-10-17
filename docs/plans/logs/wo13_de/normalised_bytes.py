"""Package wo13-de: the bytes a table of every tested gene's normalised values would take (E.1b).

For each contrast of one M5 run: the tested genes (a padj) by the contrast's samples, as
log2(count / size factor + 1) from the written size-factor table, written as parquet wide
(Float64 and Float32) and long (gene, sample, value) into a scratch folder, and measured.

Run: bash env.sh python normalised_bytes.py <root> <dataset> <scratch folder>
"""

import sys
from pathlib import Path

import polars as pl

from ctx.de.normalised import WHAT, log2_normalised
from ctx.de.run import newest

root, dataset, scratch = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
scratch.mkdir(parents=True, exist_ok=True)
de, gathered = root / "outputs" / "05_de", root / "outputs" / "01_gather"
counts = pl.read_parquet(newest(gathered, dataset, "counts"))
factors = pl.read_parquet(next(de.glob(f"*_{dataset}_{WHAT}_v1.parquet")))
print(f"size-factor table {next(de.glob(f'*_{dataset}_{WHAT}_v1.parquet')).stat().st_size:,} B")
total = {"wide64": 0, "wide32": 0, "long": 0}
for (contrast,), rows in factors.group_by("contrast", maintain_order=True):
    results = pl.read_parquet(next(de.glob(f"*_{dataset}_{contrast}_results_v1.parquet")))
    tested = results.filter(pl.col("padj").is_not_null())["gene_id"].to_list()
    long = log2_normalised(counts, rows, tested)
    wide = long.pivot(on="sample", index="gene_id", values="value")
    sizes = {}
    for name, frame in (("wide64", wide), ("wide32", wide.with_columns(pl.exclude("gene_id")
                        .cast(pl.Float32))), ("long", long)):  # fmt: skip
        path = scratch / f"{dataset}_{contrast}_{name}.parquet"
        frame.write_parquet(path)
        sizes[name] = path.stat().st_size
        total[name] += sizes[name]
    print(f"{contrast}: {len(tested):,} genes x {rows.height:,} samples = {long.height:,} values; "
          + "; ".join(f"{k} {v:,} B" for k, v in sizes.items())
          + f"; estimated_size wide64 {wide.estimated_size():,} B")  # fmt: skip
print("all contrasts:", "; ".join(f"{k} {v:,} B" for k, v in total.items()))
results_bytes = sum(p.stat().st_size for p in de.glob(f"*_{dataset}_*_results_v1.parquet"))
print(f"results tables, for scale: {results_bytes:,} B")
