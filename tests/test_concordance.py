"""ctx.concordance: the agreement statistics D-040 defines, each against a hand computation."""

import math
import warnings

import numpy as np
import polars as pl
import pytest
from scipy import stats

from ctx import concordance

# d = log2(b + 1) - log2(a + 1) = [0, 1, 1, 2, 6]: log2(a + 1) is 1..5, log2(b + 1) 1, 3, 4, 6, 11.
A = [1, 3, 7, 15, 31]
B = [1, 7, 15, 63, 2047]
# Ties in both, and no gene zero in both, so spearmanr sees exactly the genes agreement keeps.
TIED_A = [1, 1, 2, 2, 2, 5, 5, 5, 9, 9, 0, 3]
TIED_B = [1, 3, 2, 2, 4, 3, 8, 8, 8, 1, 6, 0]


# On these 17 kept genes spearmanr of the vector with itself is 0.9999999999999998.
SEVENTEEN = [0, 107, 108, 95, 95, 108, 86, 113, 102, 112, 102, 104, 92, 109, 94, 100, 104, 115, 0]
# Here spearmanr and pearsonr (on log2(count + 1)) of the vector with itself both miss 1.
FIVE = [61, 49, 56, 55, 40]


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("counts", [SEVENTEEN, FIVE])
def test_identical_vectors_agree_exactly(counts):
    """Identical counts give rho 1 and r 1 exactly, where scipy misses 1 by rounding.

    Mutations caught: the estimate from spearmanr (0.9999999999999998 on SEVENTEEN); Pearson's
    r from pearsonr (0.9999999999999999 on FIVE)."""
    kept = np.array([c for c in counts if c], dtype=float)
    assert stats.spearmanr(kept, kept).statistic != 1.0
    if counts is FIVE:
        assert stats.pearsonr(np.log2(kept + 1), np.log2(kept + 1)).statistic != 1.0
    result = concordance.agreement(counts, counts)
    assert result.n_genes == result.n_identical == kept.size
    assert (result.spearman, result.spearman_low, result.spearman_high) == (1.0, 1.0, 1.0)
    assert (result.pearson, result.pearson_low, result.pearson_high) == (1.0, 1.0, 1.0)
    assert (result.bias, result.sd, result.loa_low, result.loa_high, result.fold) == (
        0.0, 0.0, 0.0, 0.0, 1.0)  # fmt: skip
    assert (result.bias_low, result.bias_high) == (0.0, 0.0)


def test_rank_identical_counts_keep_rho_inside_its_interval():
    """One gene a read higher leaves every rank, so rho, its interval and every resample are 1.

    Mutation caught: the estimate from spearmanr, 0.9999999999999998 below its own [1, 1]."""
    other = list(SEVENTEEN)
    other[7] += 1
    kept = [[c for c in v if c] for v in (SEVENTEEN, other)]
    assert (stats.rankdata(kept[0]) == stats.rankdata(kept[1])).all()
    result = concordance.agreement(SEVENTEEN, other)
    assert result.n_identical == 16
    assert result.spearman_low <= result.spearman <= result.spearman_high
    assert (result.spearman, result.spearman_low, result.spearman_high) == (1.0, 1.0, 1.0)


def test_a_fourfold_count_correlates_perfectly_and_disagrees():
    """b + 1 = 4(a + 1) on every gene: log2(b + 1) is log2(a + 1) + 2, so r is 1 while the bias is
    2, the case Bland and Altman's limits exist for. Here the sum formula gives
    1.0000000000000002.

    Mutation caught: r not clipped to [-1, 1] (atanh of it is a math domain error)."""
    a = [13, 14, 32, 33, 35]
    result = concordance.agreement(a, [4 * (c + 1) - 1 for c in a])
    assert (result.spearman, result.spearman_low, result.spearman_high) == (1.0, 1.0, 1.0)
    assert (result.pearson, result.pearson_low, result.pearson_high) == (1.0, 1.0, 1.0)
    assert result.bias == pytest.approx(2.0) and result.fold == pytest.approx(4.0)
    assert result.n_identical == 0


def test_the_estimate_is_spearmanr_on_tied_data():
    """Mutation caught: Pearson on the raw values (0.281 here, against Spearman's 0.204);
    ordinal ranks in the bootstrap's statistic instead of average ones (0.308)."""
    result = concordance.agreement(TIED_A, TIED_B)
    expected = stats.spearmanr(TIED_A, TIED_B).statistic
    assert result.spearman == pytest.approx(expected, abs=1e-15)
    assert abs(stats.pearsonr(TIED_A, TIED_B).statistic - expected) > 0.05
    rows = np.array([TIED_A, TIED_B[::-1]], dtype=float)
    batch = concordance.rank_correlation(rows, np.array([TIED_B, TIED_A], dtype=float))
    assert batch[0] == pytest.approx(expected, abs=1e-12)
    assert batch[1] == pytest.approx(stats.spearmanr(TIED_B[::-1], TIED_A).statistic, abs=1e-12)


def test_pearson_interval_is_fisher_z_by_hand():
    """Mutation caught: a standard error of 1/sqrt(n) instead of 1/sqrt(n - 3)."""
    a = [3, 10, 0, 45, 7, 120, 16, 2, 60, 5, 33, 1]
    b = [5, 7, 1, 60, 3, 90, 30, 0, 41, 9, 20, 4]
    result = concordance.agreement(a, b)
    r = stats.pearsonr(np.log2(np.add(a, 1)), np.log2(np.add(b, 1))).statistic
    half = 1.959963984540054 / math.sqrt(len(a) - 3)
    assert result.pearson == pytest.approx(r, abs=1e-15)
    assert result.pearson_low == pytest.approx(math.tanh(math.atanh(r) - half), abs=1e-12)
    assert result.pearson_high == pytest.approx(math.tanh(math.atanh(r) + half), abs=1e-12)


def test_spearman_interval_covers_the_population_value():
    """1,000 bivariate normal datasets, n = 60, Pearson 0.6: the population Spearman is
    (6/pi) asin(0.3) = 0.5819. Measured coverage 0.960 at data seed 7.

    Mutation caught: a [5, 95] interval (coverage 0.920 in review 3's measurement)."""
    target = 6 / math.pi * math.asin(0.3)
    rng = np.random.default_rng(7)
    hits = 0
    for _ in range(1000):
        xy = rng.multivariate_normal([0, 0], [[1, 0.6], [0.6, 1]], size=60)
        # exp keeps the ranks and makes the values counts-like: positive, never both zero.
        result = concordance.agreement(np.exp(xy[:, 0]), np.exp(xy[:, 1]))
        hits += result.spearman_low <= target <= result.spearman_high
    assert 0.94 <= hits / 1000 <= 0.99


def test_duplicated_samples_leave_the_pooled_row_unchanged():
    """The pooled row compares gene sums, so a sample counted twice adds no genes.

    Mutation caught: resampling gene x sample points (the doubled frame has twice the points,
    and its interval narrows)."""
    rng = np.random.default_rng(3)
    genes = [f"g{i}" for i in range(40)]
    long = pl.DataFrame({
        "sample": [s for s in ("s1", "s2", "s3") for _ in genes],
        "gene_id": genes * 3,
        "a": rng.poisson(20, 120).astype(float),
        "b": rng.poisson(20, 120).astype(float),
    })  # fmt: skip
    doubled = pl.concat([long, long.with_columns(pl.col("sample") + "_copy")])
    once = concordance.per_unit(long, "a", "b", unit="sample", gene="gene_id")
    twice = concordance.per_unit(doubled, "a", "b", unit="sample", gene="gene_id")
    assert once["sample"].to_list() == ["s1", "s2", "s3", "all"]
    pooled = [t.filter(pl.col("sample") == "all").row(0, named=True) for t in (once, twice)]
    for column in ("n_genes", "spearman", "spearman_low", "spearman_high"):
        assert pooled[0][column] == pooled[1][column], column
    assert pooled[0]["n_genes"] == 40
    with pytest.raises(ValueError, match="named 'all'"):
        concordance.per_unit(long.with_columns(sample=pl.lit("all")), "a", "b", unit="sample",
                             gene="gene_id")  # fmt: skip


def test_a_gene_twice_in_one_unit_is_refused():
    """Its two rows would enter one comparison as two genes, and the pooled sum twice.

    Mutation caught: per_unit without its duplicate check."""
    frame = pl.DataFrame({"sample": ["s1", "s1", "s1", "s2"], "gene_id": ["g1", "g2", "g1", "g1"],
                          "a": [1.0, 2.0, 3.0, 4.0], "b": [1.0, 2.0, 3.0, 4.0]})  # fmt: skip
    with pytest.raises(ValueError, match="appears twice in one sample"):
        concordance.per_unit(frame, "a", "b", unit="sample", gene="gene_id")


def test_a_gene_zero_in_one_method_only_is_kept():
    """Mutation caught: a "nonzero in both" filter (n_genes 4, and the 0-against-50 gene's
    difference leaves the bias); keeping genes zero in both (n_genes 6)."""
    result = concordance.agreement([0, 10, 20, 30, 40, 0], [50, 10, 20, 30, 40, 0])
    assert result.n_genes == 5
    assert result.n_identical == 4
    assert result.bias == pytest.approx(math.log2(51) / 5)


def test_bland_altman_by_hand():
    """d = [0, 1, 1, 2, 6]: bias 2, SD sqrt(5.5) = 2.3452, limits 2 ± 1.96 SD; percentiles by
    linear interpolation 0.1 and 5.6; the bias's t interval with t(0.975, 4) = 2.7764.

    Mutations caught: ddof 0 (SD 2.0976); d as a - b (bias -2)."""
    result = concordance.agreement(A, B)
    sd = math.sqrt(5.5)
    assert result.bias == pytest.approx(2.0)
    assert result.sd == pytest.approx(sd)
    assert round(result.sd, 4) == 2.3452
    assert result.loa_low == pytest.approx(2 - 1.96 * sd)
    assert result.loa_high == pytest.approx(2 + 1.96 * sd)
    assert round(result.loa_low, 4) == -2.5966 and round(result.loa_high, 4) == 6.5966
    assert (result.d_p2_5, result.d_p97_5) == (pytest.approx(0.1), pytest.approx(5.6))
    half = 2.7764451051977934 * sd / math.sqrt(5)
    assert (result.bias_low, result.bias_high) == (pytest.approx(2 - half), pytest.approx(2 + half))
    assert result.fold == pytest.approx(4.0)


def intervals(result: concordance.Agreement) -> tuple:
    """Every interval end of a result."""
    return (result.spearman_low, result.spearman_high, result.pearson_low, result.pearson_high,
            result.bias_low, result.bias_high)  # fmt: skip


def test_no_interval_below_four_genes():
    """The size rules hold for identical counts too: no correlation from one gene or a constant
    vector, no spread from one gene, no interval from three.

    Mutations caught: MIN_CI = 3 (Fisher's z then divides by sqrt(0)); identical vectors
    answered before the size rules (rho 1 with [1, 1] from one gene, the old short-circuit)."""
    three = concordance.agreement([1, 5, 9], [2, 4, 12])
    assert three.spearman == pytest.approx(1.0) and three.pearson is not None
    assert intervals(three) == (None,) * 6
    assert three.loa_low is not None
    four = concordance.agreement([1, 5, 9, 3], [2, 4, 12, 3])
    assert None not in intervals(four)
    one = concordance.agreement([0, 7], [0, 7])
    assert (one.n_genes, one.spearman, one.pearson, one.sd, one.loa_low) == (1, None, None, None,
                                                                               None)  # fmt: skip
    assert (one.bias, one.fold) == (0.0, 1.0) and intervals(one) == (None,) * 6
    for counts in ([3, 7], [3, 7, 9]):
        same = concordance.agreement(counts, counts)
        assert (same.spearman, same.pearson, same.sd) == (1.0, 1.0, 0.0), counts
        assert intervals(same) == (None,) * 6, counts
    flat = concordance.agreement([5] * 5, [5] * 5)
    assert (flat.spearman, flat.pearson, flat.n_identical) == (None, None, 5)
    assert (flat.bias, flat.sd, flat.bias_low, flat.bias_high) == (0.0, 0.0, 0.0, 0.0)


def test_genes_zero_in_both_leave_nothing_to_compare():
    """Mutation caught: no empty guard in the Bland-Altman fields (the mean of no differences)."""
    empty = concordance.agreement([0, 0, 0], [0, 0, 0])
    assert (empty.n_genes, empty.n_identical) == (0, 0)
    assert {empty.spearman, empty.pearson, empty.bias, empty.sd, empty.fold} == {None}


def independent_interval(a: list, b: list, seed: int = 42) -> tuple[float, float, int]:
    """A percentile bootstrap written here: 1,000 index draws of shape (1000, n) from the seeded
    generator (scipy 1.18's draws), spearmanr per resample, finite draws, np.percentile."""
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    x, y = x[(x + y) > 0], y[(x + y) > 0]
    rows = np.random.default_rng(seed).integers(0, x.size, (1000, x.size))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", stats.ConstantInputWarning)
        draws = np.array([stats.spearmanr(x[i], y[i]).statistic for i in rows])
    finite = draws[np.isfinite(draws)]
    low, high = np.percentile(finite, [2.5, 97.5])
    return float(low), float(high), int(draws.size - finite.size)


def test_the_interval_is_a_seeded_percentile_bootstrap_of_1000():
    """D-040's parameters, against a bootstrap that shares no code with the module.

    Mutations caught: seed 43; 200 resamples; the basic method; a seed argument ignored."""
    rng = np.random.default_rng(5)
    poisson = rng.poisson(rng.gamma(0.5, 40, 300))
    larger = (poisson.tolist(), np.where(rng.random(300) < 0.3, 0, poisson + 2).tolist())
    cases = {"tied": (TIED_A, TIED_B), "four": ([1, 5, 9, 3], [2, 4, 12, 3]), "larger": larger}
    for name, (a, b) in cases.items():
        result = concordance.agreement(a, b)
        low, high, dropped = independent_interval(a, b)
        assert (result.spearman_low, result.spearman_high) == pytest.approx(
            (low, high), abs=1e-12
        ), name
        assert result.resamples_dropped == dropped, name
    assert independent_interval(*cases["four"])[2] > 0
    reseeded = concordance.agreement(*larger, seed=43)
    expected = independent_interval(*larger, seed=43)
    assert expected[:2] != pytest.approx(independent_interval(*larger)[:2], abs=1e-12)
    assert (reseeded.spearman_low, reseeded.spearman_high) == pytest.approx(expected[:2], abs=1e-12)


def test_tied_resamples_are_dropped_and_counted():
    """A resample of tied values has no rank correlation.

    Mutation caught: scipy's own interval, which a single NaN resample turns into NaN."""
    result = concordance.agreement([0, 0, 0, 0, 0, 0, 1, 1, 1, 1], [0, 1, 2, 3, 0, 0, 1, 1, 1, 1])
    assert result.resamples_dropped > 0
    assert result.spearman_low is not None and math.isfinite(result.spearman_low)
    assert result.spearman_low <= result.spearman <= result.spearman_high


def test_bad_input_is_refused():
    """Mutation caught: no check on negative counts (log2 of a negative count + 1 is NaN)."""
    with pytest.raises(ValueError, match="one length"):
        concordance.agreement([1, 2], [1, 2, 3])
    with pytest.raises(ValueError, match="non-negative"):
        concordance.agreement([1, -2], [1, 2])
