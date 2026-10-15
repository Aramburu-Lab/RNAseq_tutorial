# ruff: noqa: T201 - a proof record, run once and kept as run
"""Package S: the LUSC-against-LUAD M5 run of 2026-10-15: n, markers, PCA, library sizes."""

import os

import polars as pl

from ctx.dataset import load
from ctx.de.sanity import decide, parse

pl.Config.set_tbl_rows(40).set_fmt_str_lengths(80).set_tbl_width_chars(250).set_tbl_cols(20)
D = "/Users/rytis/ctx_lung/outputs/05_de/2026-10-15_de_tcga_lung_"
G = "/Users/rytis/ctx_lung/outputs/01_gather/2026-10-15_gather_tcga_lung_"
WT = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a8ab53e1e7f33bad4"
print(pl.read_parquet(D + "contrasts_v1.parquet"))
metrics = pl.read_parquet(D + "metrics_v1.parquet")
print(metrics.filter(pl.col("metric").str.contains("lusc_luad|permanova")))
spec = load(__import__("pathlib").Path(WT + "/conf/datasets/tcga_lung.toml"))
record = dict(metrics.select("metric", "value").rows())
for row in decide(record, parse(spec.sanity)):
    print(row["rule"], row["condition"], row["left"], row["verdict"])
sheet = pl.read_parquet(G + "samplesheet_v1.parquet")
tumours = sheet.filter(pl.col("chosen") & (pl.col("condition") == "tumour"))
print(
    tumours.group_by("project")
    .agg(
        n=pl.len(),
        library_median=pl.col("library_size").median(),
        library_q1=pl.col("library_size").quantile(0.25),
        library_q3=pl.col("library_size").quantile(0.75),
        nofeature_pct=(
            pl.col("n_nofeature") / (pl.col("n_nofeature") + pl.col("library_size")) * 100
        ).median(),
        multimapping_pct=(pl.col("n_multimapping") / pl.col("library_size") * 100).median(),
        plates=pl.col("plate").n_unique(),
        sites=pl.col("site").n_unique(),
        ffpe=(pl.col("preservation") == "FFPE").sum(),
    )
    .sort("project")
)
print(
    "plates shared by both projects",
    len(
        set(tumours.filter(pl.col("project") == "TCGA-LUAD")["plate"])
        & set(tumours.filter(pl.col("project") == "TCGA-LUSC")["plate"])
    ),
)
print("page bytes", os.path.getsize(D + "overview_v1.html"))
