"""Package M9, H.3: the lung cards' exon shares, against each tissue's GTEx median TPM.

Reads the 2026-10-14 offline build's tables in ~/ctx_lung (read only) and prints, per gene,
each exon's median share over all tissues and over the tissues where the gene's median TPM
reaches 1, HPA's detection limit and the tau floor.
"""

from pathlib import Path

import polars as pl

OUT = Path.home() / "ctx_lung/outputs/11_genecard"
for slug, symbols in (("tcga_luad", ("NKX2-1", "SFTPC", "EGFR", "MKI67")),
                      ("tcga_lusc", ("SOX2", "TP63", "KRT5", "MKI67"))):  # fmt: skip
    exons = pl.read_csv(OUT / f"2026-10-14_genecard_{slug}_exons_v1.csv")
    tissues = pl.read_csv(OUT / f"2026-10-14_genecard_{slug}_tissues_v1.csv")
    for symbol in symbols:
        x = exons.filter(pl.col("symbol") == symbol)
        tpm = tissues.filter((pl.col("symbol") == symbol) & (pl.col("sex") == "all"))
        detected = tpm.filter(pl.col("median") >= 1)["tissue"].to_list()
        print(f"{slug} {symbol}: {x['exon'].n_unique()} exons, {x['tissue'].n_unique()} tissues,"
              f" {len(detected)} at median TPM >= 1")  # fmt: skip
        every = x.group_by("exon").agg(mid_all=pl.col("share").median(),
                                       tissues_nonzero=(pl.col("share") > 0).sum())  # fmt: skip
        held = (
            x.filter(pl.col("tissue").is_in(detected))
            .group_by("exon")
            .agg(mid_detected=pl.col("share").median())
        )
        lung = x.filter(pl.col("tissue") == "Lung").select("exon", lung="share")
        table = every.join(held, on="exon", how="left").join(lung, on="exon", how="left")
        with pl.Config(tbl_rows=60):
            print(table.sort(pl.col("exon").str.slice(1).cast(pl.Int64)))
