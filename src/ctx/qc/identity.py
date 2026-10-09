"""Identity checks on counts: whether a sample's expression agrees with what is recorded of it.

The sex check moved here from ``ctx.de.run`` after the Step-1 record (CONVENTIONS, Statistics
rules): M2 judges it, M5 records it for rule R8, both through this one definition.
"""

import math

import numpy as np
import polars as pl

SEX_X = "XIST"
SEX_Y = ("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")


def sex_check(
    sheet: pl.DataFrame, counts: pl.DataFrame, genes: pl.DataFrame
) -> tuple[pl.DataFrame | None, list[tuple[str, float, str]]]:
    """Sex from XIST against five chrY genes, per sample, against the recorded sex.

    A transplanted sample is compared with its donor's sex, where one is recorded: myeloid
    cells after HSCT come from the donor.
    """
    if "sex" not in sheet.columns:
        return None, [("sex.agree", math.nan, "not computed: the sheet has no sex column")]
    ids = dict(genes.select("gene_name", "gene_id").unique("gene_name", keep="first").rows())
    if missing := [g for g in (SEX_X, *SEX_Y) if g not in ids]:
        return None, [("sex.agree", math.nan, f"not computed: no gene named {missing}")]
    samples = sheet.unique("sample", keep="first", maintain_order=True)
    names = samples["sample"].to_list()
    matrix = counts.select(names).to_numpy().astype(float)
    cpm = matrix / matrix.sum(axis=0) * 1e6
    row = {g: i for i, g in enumerate(counts["gene_id"].to_list())}
    xist = np.log2(cpm[row[ids[SEX_X]]] + 1)
    chr_y = np.log2(cpm[[row[ids[g]] for g in SEX_Y]] + 1).mean(axis=0)
    donor = samples["donor_sex"] if "donor_sex" in samples.columns else pl.Series([""] * len(names))
    expected = [
        d if d else s
        for s, d in zip(samples["sex"].fill_null(""), donor.fill_null(""), strict=True)
    ]
    frame = pl.DataFrame(
        {
            "sample": names,
            "expected": expected,
            "xist_log2cpm": xist,
            "chry_log2cpm": chr_y,
            "expression_sex": np.where(chr_y > xist, "M", "F"),
        }
    ).with_columns(
        agrees=pl.when(pl.col("expected") == "")
        .then(None)
        .otherwise(pl.col("expected") == pl.col("expression_sex"))
    )
    compared = frame.filter(pl.col("agrees").is_not_null())
    wrong = compared.filter(~pl.col("agrees"))["sample"].to_list()
    unrecorded = frame.filter(pl.col("agrees").is_null())["sample"].to_list()
    return frame, [
        (
            "sex.compared",
            float(compared.height),
            f"{len(unrecorded)} with no recorded sex: {', '.join(unrecorded)}",
        ),
        ("sex.disagree", float(len(wrong)), ", ".join(wrong) or "none"),
    ]
