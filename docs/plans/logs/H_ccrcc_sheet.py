"""Package H: what gather kept for CPTAC ccRCC, against the census ceilings."""

from collections import Counter

import polars as pl

from ctx.gather import gdc

D = "/Users/rytis/ctx_kidney/outputs/01_gather/2026-10-14_gather_cptac_ccrcc_"
sheet = pl.read_parquet(D + "samplesheet_v1.parquet")
print("sheet rows", sheet.height, "cases", sheet["case"].n_unique())
print("by condition", dict(Counter(sheet["condition"])))
print("drop reasons", dict(Counter(zip(sheet["condition"], sheet["drop_reason"], strict=True))))
kept = gdc.chosen(sheet)
print("chosen by condition", dict(Counter(kept["condition"])))
both = set(kept.filter(pl.col("condition") == "tumour")["case"]) & set(
    kept.filter(pl.col("condition") == "normal")["case"]
)
print("pairs", len(both))
for c in ("C3N-00435", "C3N-01175"):
    print(c, "in pairs", c in both)
counts = pl.read_parquet(D + "counts_v1.parquet")
print("counts", counts.height, "genes by", counts.width - 1, "files")
pooled = [
    "CPT0012900007",
    "CPT0089460006",
    "CPT0088550007",
    "CPT0066470007",
    "CPT0014370010",
    "CPT0025580007",
    "CPT0000870008",
    "CPT0014450014",
    "CPT0071150007",
]
p = sheet.filter(pl.col("sample").is_in(pooled))
print("pooled-piece aliquots", p.height, "chosen", p.filter(pl.col("drop_reason") == "").height)
