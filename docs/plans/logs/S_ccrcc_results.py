# ruff: noqa: T201 - a proof record, run once and kept as run
"""Package S: the ccRCC M7 verdicts and the metrics they read, from the 2026-10-15 run."""

import polars as pl

pl.Config.set_tbl_rows(80).set_fmt_str_lengths(90).set_tbl_width_chars(200)
P = "/Users/rytis/ctx_kidney/outputs/07_proteo/2026-10-15_proteo_cptac_ccrcc_"
sanity = pl.read_parquet(P + "sanity_v1.parquet")
print(sanity.select("rule", "condition", "left", "right", "verdict"))
metrics = pl.read_parquet(P + "metrics_v1.parquet")
print(metrics.filter(~pl.col("metric").str.starts_with("gene.")))
genes = pl.read_parquet(P + "genes_v1.parquet")
print(
    "genes in both layers",
    genes.height,
    "tested tumour",
    genes["tested"].sum(),
    "tested normal",
    genes["normal_tested"].sum(),
)
print(
    genes.filter(pl.col("gene") == "CA9").select(
        "gene", "n_quantified", "tested", "rho", "normal_n_quantified", "normal_tested"
    )
)
