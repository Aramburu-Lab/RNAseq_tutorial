# ruff: noqa: T201 - a proof record, run once and kept as run
"""Package S: what M2 (ctx qc) measured on the combined lung sheet, per project."""

import polars as pl

pl.Config.set_tbl_rows(40).set_tbl_width_chars(250).set_tbl_cols(20)
Q = "/Users/rytis/ctx_lung/outputs/02_qc/2026-10-15_qc_tcga_lung_"
G = "/Users/rytis/ctx_lung/outputs/01_gather/2026-10-15_gather_tcga_lung_samplesheet_v1.parquet"
metrics = pl.read_parquet(Q + "metrics_v1.parquet")
print("qc rows", metrics.height, "columns", metrics.columns)
sheet = pl.read_parquet(G)
print("sheet files", sheet.height, "chosen", sheet["chosen"].sum())
key = "sample" if "sample" in metrics.columns else metrics.columns[0]
joined = metrics.join(
    sheet.select("sample", "project", "condition", "chosen"),
    left_on=key,
    right_on="sample",
    how="left",
)
print("rows by chosen", joined.group_by("chosen").len().rows())
numeric = [c for c, t in metrics.schema.items() if t.is_numeric()]
print(
    joined.filter(pl.col("chosen") & (pl.col("condition") == "tumour"))
    .group_by("project")
    .agg(pl.len(), *(pl.col(c).median() for c in numeric))
    .sort("project")
)
print(pl.read_parquet(Q + "confounding_v1.parquet"))
