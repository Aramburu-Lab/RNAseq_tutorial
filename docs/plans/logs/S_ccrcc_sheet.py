# ruff: noqa: T201 - a proof record, run once and kept as run
"""Package S: what gather's 2026-10-15 ccRCC sheet holds, against the files' own ceilings."""

from collections import Counter

import polars as pl

G = "/Users/rytis/ctx_kidney/outputs/01_gather/2026-10-15_gather_cptac_ccrcc_"
sheet = pl.read_parquet(G + "samplesheet_v1.parquet")
aliquots = pl.read_parquet(G + "aliquots_v1.parquet")
print(
    "files",
    sheet.height,
    "cases",
    sheet["case"].n_unique(),
    "projects",
    sorted(set(sheet["project"])),
)
print("drop reasons", dict(Counter(sheet["drop_reason"])))
chosen = sheet.filter(pl.col("chosen"))
print("chosen by condition", dict(Counter(chosen["condition"])))
tumour = set(chosen.filter(pl.col("condition") == "tumour")["case"])
normal = set(chosen.filter(pl.col("condition") == "normal")["case"])
print("pairs", len(tumour & normal))
excluded = sheet.filter(pl.col("drop_reason") == "excluded")
print("excluded files", excluded.height, "cases", sorted(set(excluded["case"])))
print(
    "aliquots",
    aliquots.height,
    "cases",
    aliquots["case"].n_unique(),
    dict(Counter(aliquots["sample_type"])),
)
