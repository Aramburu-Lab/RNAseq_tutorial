"""Agreement between two measurements of the same genes: correlation beside Bland-Altman (D-040).

Neutral about where the two vectors come from, so M3's counters and M7's RNA against protein
share one definition. A high correlation says the two order the genes alike; Bland and
Altman's limits say how far apart they are, which is the question when one method stands in
for another (Bland and Altman 1986). Both are on the page, and no hypothesis is tested.

Every interval is nominal: the genes are a chromosome's whole annotated set, not a sample of
genes, so a page leads with the spread across samples and gives these beside it.
"""

import math
import warnings
from collections.abc import Sequence
from dataclasses import asdict, dataclass, fields

import numpy as np
import polars as pl
from numpy.typing import ArrayLike, NDArray
from scipy import stats

SEED = 42
RESAMPLES = 1000
LEVEL = 0.95
# Fisher's z has no standard error below four genes (1/sqrt(n - 3)), and three pairs have
# at most 10 distinct bootstrap resamples, so no interval is given below four.
MIN_CI = 4
# Bland and Altman's limits of agreement: bias ± 1.96 SD of the differences.
LIMIT_Z = 1.96
# The unit holding each gene's sum over every unit.
ALL = "all"

Vector = NDArray[np.float64]


@dataclass(frozen=True)
class Agreement:
    """Two measurements of the same genes compared; None where a value is undefined.

    Attributes:
        n_genes: Genes with a count in either measurement; a gene zero in both is dropped.
        n_identical: Of those, genes with the same count in both.
        spearman: Spearman's rho on the counts.
        spearman_low: Lower end of its percentile bootstrap interval over genes.
        spearman_high: Upper end of that interval.
        resamples_dropped: Bootstrap resamples with no rank correlation (all ties), left out.
        pearson: Pearson's r on log2(count + 1).
        pearson_low: Lower end of its Fisher z interval.
        pearson_high: Upper end of that interval.
        bias: Mean of d = log2(b + 1) - log2(a + 1).
        bias_low: Lower end of the bias's t interval.
        bias_high: Upper end of that interval.
        sd: Standard deviation of d, with n - 1 in the denominator.
        loa_low: bias - 1.96 sd.
        loa_high: bias + 1.96 sd.
        d_p2_5: The 2.5th percentile of d.
        d_p97_5: The 97.5th percentile of d.
        fold: 2 ** bias, a ratio of (count + 1) values (Euser et al. 2008).
    """

    n_genes: int
    n_identical: int
    spearman: float | None
    spearman_low: float | None
    spearman_high: float | None
    resamples_dropped: int
    pearson: float | None
    pearson_low: float | None
    pearson_high: float | None
    bias: float | None
    bias_low: float | None
    bias_high: float | None
    sd: float | None
    loa_low: float | None
    loa_high: float | None
    d_p2_5: float | None
    d_p97_5: float | None
    fold: float | None


SCHEMA = {f.name: pl.Int64 if f.type is int else pl.Float64 for f in fields(Agreement)}


def agreement(
    a: ArrayLike, b: ArrayLike, *, seed: int = SEED, resamples: int = RESAMPLES
) -> Agreement:
    """Compare two count vectors over the same genes, in the same order.

    Args:
        a: The first method's counts, one per gene; fractions allowed.
        b: The second method's counts for the same genes.
        seed: The bootstrap's random seed.
        resamples: Bootstrap resamples of genes for Spearman's interval.

    Returns:
        The comparison; d is b minus a.

    Raises:
        ValueError: If the vectors differ in shape or hold a negative or non-finite count.
    """
    x, y = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if x.ndim != 1 or x.shape != y.shape:
        raise ValueError(f"two vectors of one length needed, not {x.shape} and {y.shape}")
    if not (np.isfinite(x).all() and np.isfinite(y).all() and (x >= 0).all() and (y >= 0).all()):
        raise ValueError("counts must be finite and non-negative")
    # Filtered before any transform: a gene zero in one method only is a disagreement.
    kept = (x + y) > 0
    x, y = x[kept], y[kept]
    d = np.log2(y + 1) - np.log2(x + 1)
    rho = spearman(x, y, seed=seed, resamples=resamples)
    return Agreement(int(x.size), int((x == y).sum()), *rho,
                     *_pearson(np.log2(x + 1), np.log2(y + 1)), *_bland_altman(d))  # fmt: skip


def per_unit(frame: pl.DataFrame, a: str, b: str, *, unit: str, gene: str) -> pl.DataFrame:
    """One agreement per unit, and one for ``ALL``: each gene's sums over every unit.

    The pooled row compares gene sums, so its genes stay the resampled unit; stacking every
    unit's genes would count one gene once per unit and narrow the interval with each one.

    Args:
        frame: Long form, one row per gene and unit.
        a: The first method's column.
        b: The second method's column.
        unit: The unit column, usually the sample.
        gene: The gene column.

    Returns:
        The unit, then the fields of ``Agreement``, units in first-seen order and ``ALL`` last.

    Raises:
        ValueError: If a unit is named ``ALL``, or a gene appears twice in one unit.
    """
    if ALL in set(frame[unit]):
        raise ValueError(f"a {unit} is named {ALL!r}, the name of the pooled row")
    if frame.select(unit, gene).is_duplicated().any():
        raise ValueError(f"a {gene} appears twice in one {unit}")
    rows = [{unit: name, **asdict(agreement(part[a], part[b]))}
            for (name,), part in frame.group_by(unit, maintain_order=True)]  # fmt: skip
    sums = frame.group_by(gene, maintain_order=True).agg(pl.col(a).sum(), pl.col(b).sum())
    rows.append({unit: ALL, **asdict(agreement(sums[a], sums[b]))})
    return pl.DataFrame(rows, schema={unit: pl.String, **SCHEMA})


def rank_correlation(x: Vector, y: Vector, axis: int = -1) -> Vector:
    """Spearman's rho along an axis: Pearson's r on average ranks, as ``spearmanr`` defines it.

    Vectorised, so the bootstrap ranks all its resamples at once; a row of tied values has no
    rank correlation and gives NaN.
    """
    return correlation(stats.rankdata(x, axis=axis), stats.rankdata(y, axis=axis), axis)


def correlation(x: Vector, y: Vector, axis: int = -1) -> Vector:
    """Pearson's r along an axis, the one formula behind both estimates and every resample.

    Exactly 1 on identical input, where ``spearmanr`` and ``pearsonr`` of a vector with itself
    can give 0.9999999999999998; one formula also keeps an estimate inside its own interval.
    Clipped to [-1, 1], as ``pearsonr`` is, because rounding can carry an exact linear relation
    past 1, where atanh is undefined.
    """
    xc = x - x.mean(axis=axis, keepdims=True)
    yc = y - y.mean(axis=axis, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (xc * yc).sum(axis=axis) / np.sqrt((xc * xc).sum(axis=axis) * (yc * yc).sum(axis=axis))
    result: Vector = np.clip(r, -1.0, 1.0)
    return result


def _constant(values: Vector) -> bool:
    """Whether a vector has one value only, where no correlation is defined."""
    return bool(values.size == 0 or (values == values[0]).all())


def spearman(
    a: ArrayLike, b: ArrayLike, *, seed: int = SEED, resamples: int = RESAMPLES
) -> tuple[float | None, float | None, float | None, int]:
    """Spearman's rho, its percentile bootstrap interval over genes, and resamples dropped.

    Signed values are allowed, so M7 compares fold changes with it; ``agreement()`` calls it on
    counts. None where undefined: below two pairs or on a constant vector, and no interval below
    ``MIN_CI`` pairs or when every resample is tied.

    Raises:
        ValueError: If the vectors differ in shape or hold a non-finite value.
    """
    x, y = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if x.ndim != 1 or x.shape != y.shape:
        raise ValueError(f"two vectors of one length needed, not {x.shape} and {y.shape}")
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("values must be finite")
    if x.size < 2 or _constant(x) or _constant(y):  # noqa: PLR2004 - a correlation needs two
        return None, None, None, 0
    rho = float(rank_correlation(x, y))
    if x.size < MIN_CI:
        return rho, None, None, 0
    with warnings.catch_warnings():
        # scipy's quantile turns one NaN resample into a NaN interval and warns; those
        # resamples are dropped and counted below instead.
        warnings.simplefilter("ignore", stats.DegenerateDataWarning)
        result = stats.bootstrap(
            (x, y), rank_correlation, n_resamples=resamples, vectorized=True, paired=True,
            method="percentile", confidence_level=LEVEL, rng=np.random.default_rng(seed),
        )  # fmt: skip
    draws = np.asarray(result.bootstrap_distribution, dtype=np.float64)
    finite = draws[np.isfinite(draws)]
    dropped = int(draws.size - finite.size)
    if not finite.size:
        return rho, None, None, dropped
    # numpy's default linear quantile is scipy.stats.quantile's, which bootstrap itself uses.
    low, high = np.quantile(finite, [(1 - LEVEL) / 2, (1 + LEVEL) / 2])
    return rho, float(low), float(high), dropped


def _pearson(x: Vector, y: Vector) -> tuple[float | None, float | None, float | None]:
    """Pearson's r and its Fisher z interval, tanh(atanh(r) ± z / sqrt(n - 3))."""
    if x.size < 2 or _constant(x) or _constant(y):  # noqa: PLR2004 - a correlation needs two
        return None, None, None
    r = float(correlation(x, y))
    if x.size < MIN_CI:
        return r, None, None
    if abs(r) == 1:
        # atanh(±1) is infinite: the interval collapses to the point.
        return r, r, r
    half = float(stats.norm.ppf((1 + LEVEL) / 2)) / math.sqrt(x.size - 3)
    return r, math.tanh(math.atanh(r) - half), math.tanh(math.atanh(r) + half)


def _bland_altman(d: Vector) -> Sequence[float | None]:
    """Bias with its t interval, SD (ddof 1), limits, d's 2.5th and 97.5th percentiles, fold."""
    if not d.size:
        return [None] * 9
    bias = float(d.mean())
    low, high = np.percentile(d, [2.5, 97.5])
    if d.size < 2:  # noqa: PLR2004 - a spread needs two
        return [bias, None, None, None, None, None, float(low), float(high), 2**bias]
    sd = float(d.std(ddof=1))
    half = None
    if d.size >= MIN_CI:
        half = float(stats.t.ppf((1 + LEVEL) / 2, d.size - 1)) * sd / math.sqrt(d.size)
    return [
        bias, None if half is None else bias - half, None if half is None else bias + half,
        sd, bias - LIMIT_Z * sd, bias + LIMIT_Z * sd, float(low), float(high), 2**bias,
    ]  # fmt: skip
