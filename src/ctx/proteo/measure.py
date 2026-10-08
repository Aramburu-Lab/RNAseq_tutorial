"""M7's statistics, each defined once; the order of ``per_gene`` is the plan's and is tested.

``per_gene`` runs, in this order: one row per gene (the readers, ``tables.one_per_gene``), the
quantified filter decided in exact arithmetic, a correlation on each gene's pairwise-complete
samples, then Benjamini-Hochberg over the genes with a correlation. A missing protein value is
never a zero: it leaves that sample out of that gene's test (CONVENTIONS.md, proteomics).
"""

import itertools
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
import polars as pl
from scipy import stats

from ctx.concordance import correlation, rank_correlation
from ctx.de.stats import adjust

MIN_QUANTIFIED = 0.7
METHODS = ("spearman", "pearson")
# A correlation's t statistic needs n - 2 degrees of freedom; a pair of points fixes r at 1.
MIN_PAIRS = 3
# The class every gene of no declared class falls in, and one in two declared classes.
OTHER = "other"


def enough(present: np.ndarray, total: int, minimum: float) -> np.ndarray:
    """Whether each count of quantified samples reaches the share, decided as fractions.

    ``Fraction(str(0.7))`` is 7/10 exactly, where the float 0.7 is a hair below it, so 7 of 10
    passes and 69 of 99 cannot pass by rounding (lessons.md, exact arithmetic at a band's edge).
    """
    need = Fraction(str(minimum))
    return np.array([Fraction(int(k), total) >= need for k in present], dtype=bool)


def _pair(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float]:
    """One gene's correlation and two-sided p on its complete pairs; NaN where undefined.

    The p is the t-distribution's, as ``scipy.stats.spearmanr`` and ``pearsonr`` give it.
    """
    keep = ~(np.isnan(x) | np.isnan(y))
    x, y = x[keep], y[keep]
    if x.size < MIN_PAIRS or np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan, np.nan
    r = float(rank_correlation(x, y) if method == "spearman" else correlation(x, y))
    if abs(r) == 1:
        return r, 0.0
    t = r * np.sqrt((x.size - 2) / (1 - r * r))
    return r, float(2 * stats.t.sf(abs(t), x.size - 2))


def per_gene(genes: tuple[str, ...], rna: np.ndarray, protein: np.ndarray,
             minimum: float = MIN_QUANTIFIED) -> pl.DataFrame:  # fmt: skip
    """Per-gene RNA-protein correlation across samples, both methods, BH over tested genes.

    Args:
        genes: One symbol per row, unique: the readers' one row per gene comes first.
        rna: Genes by samples, the same samples in the same order as ``protein``.
        protein: Genes by samples, NaN where missing.
        minimum: Share of samples a protein needs a value in to be tested.

    Returns:
        One row per gene: quantified count and share, whether tested, and for each method the
        estimate, its p and its BH-adjusted p (null for a gene not tested).

    Raises:
        ValueError: If the shapes differ or a gene appears twice.
    """
    if rna.shape != protein.shape or len(genes) != rna.shape[0]:
        raise ValueError(f"RNA {rna.shape}, protein {protein.shape}, {len(genes)} genes")
    if len(set(genes)) != len(genes):
        raise ValueError("a gene appears twice: collapse to one row per gene first")
    present = (~np.isnan(protein)).sum(axis=1)
    tested = enough(present, protein.shape[1], minimum)
    share = present / protein.shape[1]
    frame = pl.DataFrame({"gene": list(genes), "n_quantified": present.astype(np.int64),
                          "quantified_share": share, "tested": tested})  # fmt: skip
    for method, short in (("spearman", "rho"), ("pearson", "r")):
        found = np.full((len(genes), 2), np.nan)
        for i in np.flatnonzero(tested):
            found[i] = _pair(rna[i], protein[i], method)
        padj = np.full(len(genes), np.nan)
        defined = ~np.isnan(found[:, 1])
        if defined.any():
            padj[defined] = adjust(found[defined, 1], "padj")
        frame = frame.with_columns(
            pl.Series(short, found[:, 0]).fill_nan(None),
            pl.Series(f"{short}_p", found[:, 1]).fill_nan(None),
            pl.Series(f"{short}_padj", padj).fill_nan(None),
        )
    return frame


def membership(genes: list[str], classes: Mapping[str, set[str]]) -> tuple[list[str], int]:
    """Each gene's class, ``OTHER`` for none; a gene in two classes is left out of both.

    Returns:
        The class per gene (empty for a gene in two classes), and how many such genes.
    """
    labels, twice = [], 0
    for gene in genes:
        hits = [name for name, members in classes.items() if gene in members]
        if len(hits) > 1:
            twice += 1
        labels.append(hits[0] if len(hits) == 1 else ("" if hits else OTHER))
    return labels, twice


@dataclass(frozen=True)
class ClassTest:
    """Kruskal-Wallis across classes with Dunn's pairs, BH over the pairs.

    Attributes:
        n: Genes per class, in the order given.
        h: The H statistic, tie-corrected.
        p: Its p.
        epsilon2: H / (N - 1), the share of rank variance between classes.
        pairs: Rows of (class a, class b, z, p, BH-adjusted p).
    """

    n: dict[str, int]
    h: float
    p: float
    epsilon2: float
    pairs: tuple[tuple[str, str, float, float, float], ...]


def class_test(groups: Mapping[str, np.ndarray]) -> ClassTest:
    """Kruskal-Wallis and Dunn's test on the values of each class (Dunn 1964, tie-corrected)."""
    names = [k for k, v in groups.items() if v.size]
    values = [groups[k] for k in names]
    h, p = stats.kruskal(*values)
    pooled = np.concatenate(values)
    ranks = stats.rankdata(pooled)
    total = pooled.size
    _, ties = np.unique(pooled, return_counts=True)
    spread = total * (total + 1) / 12 - (ties**3 - ties).sum() / (12 * (total - 1))
    edges = np.cumsum([0, *(v.size for v in values)])
    mean = {k: ranks[edges[i] : edges[i + 1]].mean() for i, k in enumerate(names)}
    sizes = {k: v.size for k, v in zip(names, values, strict=True)}
    raw = []
    for a, b in itertools.combinations(names, 2):
        z = (mean[a] - mean[b]) / np.sqrt(spread * (1 / sizes[a] + 1 / sizes[b]))
        raw.append((a, b, float(z), float(2 * stats.norm.sf(abs(z)))))
    padj: list[float] = adjust(np.array([r[3] for r in raw]), "padj").tolist() if raw else []
    pairs = tuple((*r, float(q)) for r, q in zip(raw, padj, strict=True))
    return ClassTest(sizes, float(h), float(p), float(h / (total - 1)), pairs)


def shifted(higher: np.ndarray, lower: np.ndarray) -> tuple[float, float]:
    """One-sided Mann-Whitney p that ``higher`` sits above ``lower``, and the rank-biserial r.

    r = 2 U / (n1 n2) - 1 with U of ``higher``: +1 when every value of it is the larger.
    """
    result = stats.mannwhitneyu(higher, lower, alternative="greater")
    return float(result.pvalue), float(2 * result.statistic / (higher.size * lower.size) - 1)


def welch(values: np.ndarray, tested: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per row: mean tested minus mean reference over present values, and Welch's p.

    scipy gives NaN where a group has fewer than two values or both groups are constant.
    """
    a, b = values[:, tested], values[:, ~tested]
    with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lfc = np.nanmean(a, axis=1) - np.nanmean(b, axis=1)
        result = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    return np.asarray(lfc, dtype=np.float64), np.asarray(result.pvalue, dtype=np.float64)


def median(values: pl.Series) -> float:
    """The median of a column's present values, NaN when none is present."""
    present = values.drop_nulls().drop_nans() if values.dtype.is_float() else values.drop_nulls()
    return float(np.median(present.to_numpy())) if present.len() else float("nan")


def tertiles(abundance: np.ndarray) -> np.ndarray:
    """Each value's abundance tertile, 1 low to 3 high, by numpy's linear percentiles."""
    cuts = np.percentile(abundance, [100 / 3, 200 / 3])
    return np.asarray(1 + (abundance > cuts[0]).astype(int) + (abundance > cuts[1]).astype(int))
