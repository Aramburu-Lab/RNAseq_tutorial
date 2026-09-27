"""The one up / down / not-significant rule, used by the volcano, the captions and the CSVs."""

import polars as pl

ALPHA = 0.05
LEVELS = ("up", "down", "ns")


def classify(padj: pl.Expr, lfc: pl.Expr, alpha: float = ALPHA) -> pl.Expr:
    """Call each gene up, down or ns from its adjusted p and the sign of its raw LFC.

    A missing padj (Cook's outlier or independent filtering) is ns: the gene was not tested.
    The raw LFC gives the sign because the call comes from the Wald test on the raw fit;
    shrinkage can pull an estimate to zero but never across it.

    Args:
        padj: BH-adjusted p-values.
        lfc: Raw (unshrunk) log2 fold changes.
        alpha: Significance threshold on padj, strict.

    Returns:
        An expression of ``LEVELS`` values.
    """
    called = padj.is_not_null() & padj.is_not_nan() & (padj < alpha)
    return (
        pl.when(called & (lfc > 0))
        .then(pl.lit("up"))
        .when(called & (lfc < 0))
        .then(pl.lit("down"))
        .otherwise(pl.lit("ns"))
    )
