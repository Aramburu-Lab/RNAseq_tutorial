"""F5: the two libraries of one tumour, rRNA-depleted against the other, by what they count.

A few TCGA tumours were sequenced twice, once after rRNA depletion. Their two STAR files are
the one place a library type is compared within a sample, so F5 reads both, though
``gdc.chosen`` takes one file per patient: it is the deliberate exception, stated on the page.
Files that the rule drops for a reason other than the tie-break (FFPE, another protocol, a
redaction) stay out. Shares are of each file's counted reads; the replication-dependent histone
genes are named by HGNC's cluster symbols, as the plan measured them.
"""

import re

import numpy as np
import polars as pl
from scipy.stats import wilcoxon

from ctx.ncrna import biotypes

RRNA_DEPLETION = "rRNA Depletion"
# HGNC symbols of the replication-dependent histone clusters: H1 and the core histones.
HISTONE = re.compile(r"^(H1-[0-6]|H2AC\d+|H2BC\d+|H3C\d+|H4C\d+)$")
CLASSES = ("histone", biotypes.SMALL, biotypes.MITO)
KEPT_REASONS = ("", "not_first")


def pairs(files: pl.DataFrame, counts: pl.DataFrame, genes: pl.DataFrame) -> pl.DataFrame:
    """Each tumour sample with an rRNA-depleted and another library: both files' shares.

    Args:
        files: Gather's whole STAR sheet, every file, with library_selection.
        counts: STAR counts, gene_id then one column per aliquot.
        genes: The annotation's record of each count row, in count order, with biotype.

    Returns:
        sample_id, library (rRNA-depleted or the other selection, as GDC names it), sample
        (the aliquot) and the percent of counted reads of each class in ``CLASSES``; empty
        where no sample has both.
    """
    usable = files.filter(pl.col("condition") == "tumour", pl.col("drop_reason").is_in(
        list(KEPT_REASONS)), pl.col("library_selection").is_not_null())  # fmt: skip
    kinds = usable.group_by("sample_id").agg(
        (pl.col("library_selection") == RRNA_DEPLETION).any().alias("rrna"),
        (pl.col("library_selection") != RRNA_DEPLETION).any().alias("other"))  # fmt: skip
    both = kinds.filter("rrna", "other")["sample_id"]
    chosen = usable.filter(pl.col("sample_id").is_in(both.implode())).sort("sample_id", "sample")
    row = {g: i for i, g in enumerate(counts["gene_id"].to_list())}
    named = genes.with_columns(histone=pl.col("gene_name").str.contains(HISTONE.pattern))
    masks = {"histone": named["histone"].to_numpy(),
             biotypes.SMALL: (named["biotype"] == biotypes.SMALL).to_numpy(),
             biotypes.MITO: (named["biotype"] == biotypes.MITO).to_numpy()}  # fmt: skip
    order = np.array([row[g] for g in genes["gene_id"]])
    matrix = counts.select(chosen["sample"].to_list()).to_numpy()[order].astype(float)
    total = matrix.sum(axis=0)
    shares = {c: 100 * matrix[m].sum(axis=0) / total for c, m in masks.items()}
    return chosen.select("sample_id", "sample", "library_selection").with_columns(
        **{c: pl.Series(shares[c]) for c in CLASSES}
    )


def compare(paired: pl.DataFrame) -> dict[str, dict[str, float]]:
    """Per class: median of each library type, rRNA-depleted higher in how many, Wilcoxon p.

    A sample with two files of one type keeps the first by aliquot ID, as the rule does.
    """
    out: dict[str, dict[str, float]] = {}
    first = paired.unique(["sample_id", "library_selection"], keep="first", maintain_order=True)
    rrna = first.filter(pl.col("library_selection") == RRNA_DEPLETION).sort("sample_id")
    other = (first.filter(pl.col("library_selection") != RRNA_DEPLETION)
             .unique("sample_id", keep="first", maintain_order=True).sort("sample_id"))  # fmt: skip
    for c in CLASSES:
        a, b = rrna[c].to_numpy(), other[c].to_numpy()
        p = float(wilcoxon(a, b).pvalue) if a.size else float("nan")
        out[c] = {"pairs": float(a.size), "higher": float((a > b).sum()),
                  "rrna": float(np.median(a)) if a.size else float("nan"),
                  "other": float(np.median(b)) if b.size else float("nan"), "p": p}  # fmt: skip
    return out


def wide(paired: pl.DataFrame) -> pl.DataFrame:
    """One row per sample: each class's share in the rRNA-depleted and in the other library."""
    first = paired.unique(["sample_id", "library_selection"], keep="first", maintain_order=True)
    kind = pl.when(pl.col("library_selection") == RRNA_DEPLETION).then(pl.lit("rrna")).otherwise(
        pl.lit("other"))  # fmt: skip
    long = first.with_columns(kind=kind).unique(["sample_id", "kind"], keep="first",
                                                maintain_order=True)  # fmt: skip
    return long.pivot(on="kind", index="sample_id", values=list(CLASSES)).sort("sample_id")


Metrics = dict[str, tuple[float, str]]


def measured(files: pl.DataFrame) -> bool:
    """Whether the sheet carries the read groups' library selection at all."""
    column = "library_selection"
    return column in files.columns and files[column].drop_nulls().len() > 0
