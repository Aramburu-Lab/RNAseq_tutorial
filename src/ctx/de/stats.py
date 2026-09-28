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


def p_line(basis: str, alpha: float, n_tested: int, n_passing: int) -> float | None:
    """The raw p the cutoff sits at: the volcano's horizontal line, for any basis.

    A tested gene has ``basis`` below ``alpha`` if and only if its raw p is below this value.
    The value is the bar the weakest passing gene cleared in the procedure itself (rank k of
    BH or BY, step k of Holm), or the first bar when no gene passes. So it is proportional to
    ``alpha`` while no gene crosses, and the line moves with every cutoff, where the weakest
    passing gene's own p stood still between cutoffs that call the same genes.

    Args:
        basis: ``"pvalue"`` for the raw p, or a key of ``ADJUSTMENTS``.
        alpha: The cutoff on that basis, strict.
        n_tested: Genes the adjustment counted over.
        n_passing: Of those, the genes whose ``basis`` is below ``alpha``.

    Returns:
        The line on the raw-p scale, or None when no gene was tested.

    Raises:
        ValueError: If ``basis`` is unknown.
    """
    if n_tested == 0:
        return None
    m, rank = n_tested, max(n_passing, 1)
    if basis == "pvalue":
        return alpha
    if basis == "padj":
        return rank * alpha / m
    if basis == "padj_by":
        # The constant scipy's BY multiplies by, summed as scipy sums it.
        harmonic = float(np.sum(1.0 / np.arange(1, m + 1)))
        return rank * alpha / (m * harmonic)
    if basis == "padj_holm":
        return alpha / (m - rank + 1)
    if basis == "padj_bonferroni":
        return alpha / m
    raise ValueError(f"unknown basis {basis!r}; pvalue or one of {list(ADJUSTMENTS)}")
