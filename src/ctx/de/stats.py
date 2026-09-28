"""The one rule that calls a gene, and the p-value adjustments the pages offer.

numpy and scipy only at import, so the same definitions run in Pyodide for the JupyterLite
notebooks, which have no polars; ``classify`` imports polars itself, and only the CLI calls it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.stats import false_discovery_control

if TYPE_CHECKING:
    import polars as pl

ALPHA = 0.05
# A doubling or a halving, in log2 fold change.
LFC_MIN = 1.0
LEVELS = ("up", "down", "ns")
# The smallest p-value a figure draws: -log10 of 0, which pydeseq2 returns for the strongest
# effects, has no place on an axis. Anything below is drawn at 300 and the table keeps it.
P_FLOOR = 1e-300
# Results column to method, in the order the pages list them. ``padj`` is pydeseq2's own.
ADJUSTMENTS = {
    "padj": "Benjamini-Hochberg",
    "padj_by": "Benjamini-Yekutieli",
    "padj_holm": "Holm",
    "padj_bonferroni": "Bonferroni",
}


def classify(
    padj: pl.Expr, lfc: pl.Expr, alpha: float = ALPHA, lfc_min: float = LFC_MIN
) -> pl.Expr:
    """Call each gene up, down or ns from its adjusted p and its raw LFC.

    The fold change is filtered after the Wald test, not tested (DESeq2's ``lfcThreshold``):
    a call says the estimate moved at least ``lfc_min``, not that the true change provably
    exceeds it. The raw LFC is the estimate the Wald test judged, so the size filter and the
    p-value describe the same number. A missing padj (Cook's outlier or independent
    filtering) is ns: the gene was not tested.

    Args:
        padj: Adjusted p-values.
        lfc: Raw (unshrunk) log2 fold changes.
        alpha: Significance threshold on padj, strict.
        lfc_min: Smallest absolute log2 fold change called, inclusive.

    Returns:
        An expression of ``LEVELS`` values.
    """
    import polars as pl  # noqa: PLC0415 - keeps the module importable where polars is absent

    passes = padj.is_not_null() & padj.is_not_nan() & (padj < alpha)
    return (
        pl.when(passes & (lfc >= lfc_min))
        .then(pl.lit("up"))
        .when(passes & (lfc <= -lfc_min))
        .then(pl.lit("down"))
        .otherwise(pl.lit("ns"))
    )


def adjust(p: np.ndarray, method: str) -> np.ndarray:
    """Adjust raw p-values for multiple testing, over exactly the p-values given.

    BH and BY use the scipy function pydeseq2 uses for its own ``padj``, so BH here equals
    pydeseq2's column on the same genes.

    Args:
        p: Raw p-values, 1-D, each in [0, 1]; no NaN.
        method: A key of ``ADJUSTMENTS``.

    Returns:
        Adjusted p-values, in the input order.

    Raises:
        ValueError: If ``method`` is unknown or ``p`` is not 1-D in [0, 1].
    """
    if p.ndim != 1 or not np.all((p >= 0) & (p <= 1)):
        raise ValueError("p must be 1-D with every value in [0, 1] and no NaN")
    n = p.size
    if method == "padj":
        return np.asarray(false_discovery_control(p, method="bh"), dtype=float)
    if method == "padj_by":
        return np.asarray(false_discovery_control(p, method="by"), dtype=float)
    if method == "padj_bonferroni":
        return np.asarray(np.minimum(p * n, 1.0), dtype=float)
    if method == "padj_holm":
        order = np.argsort(p, kind="stable")
        # The cumulative max makes the step-down monotone: a smaller p never adjusts higher.
        stepped = np.minimum(np.maximum.accumulate(p[order] * (n - np.arange(n))), 1.0)
        adjusted = np.empty(n)
        adjusted[order] = stepped
        return adjusted
    raise ValueError(f"unknown adjustment {method!r}; one of {list(ADJUSTMENTS)}")


def p_cutoff(p: np.ndarray, adjusted: np.ndarray, alpha: float) -> float | None:
    """The largest raw p whose adjusted p is below alpha: the volcano's horizontal line.

    Exact, not an approximation, because every adjustment in ``ADJUSTMENTS`` is monotone
    non-decreasing in p: a gene passes if and only if its raw p is at most this value. It
    holds over the genes the adjustment was computed on, and nowhere else.

    Args:
        p: Raw p-values.
        adjusted: Their adjusted p-values, same order.
        alpha: Threshold on the adjusted p, strict.

    Returns:
        The cutoff on raw p, or None when no gene passes.
    """
    passing = p[adjusted < alpha]
    return float(passing.max()) if passing.size else None
