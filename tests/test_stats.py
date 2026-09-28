"""The call rule and the four p-value adjustments: oracles, boundaries and the volcano cutoff."""

import subprocess
import sys

import numpy as np
import polars as pl
import pytest
from synth import design_counts

from ctx.dataset import Contrast
from ctx.de.deseq import fit
from ctx.de.stats import ADJUSTMENTS, ALPHA, LFC_MIN, adjust, classify, p_line

BASIC = Contrast("basic", "condition", "control", "case", (6, 6))
# statsmodels names for the same four methods: the oracle, installed but not a ctx dependency.
STATSMODELS = {
    "padj": "fdr_bh",
    "padj_by": "fdr_by",
    "padj_holm": "holm",
    "padj_bonferroni": "bonferroni",
}
EPS = 1e-9


@pytest.fixture(scope="module")
def filtered():
    """A fit where Cook's and independent filtering both remove genes, as in test_na_reasons.

    Without removed genes, BH over every p and BH over the tested set are the same numbers,
    and the test could not tell them apart.
    """
    synthetic = design_counts(condition_lfc=1.0)
    rng = np.random.default_rng(1)
    names = synthetic.samples["sample"].to_list()
    low = pl.DataFrame(
        {"gene_id": [f"low{i:03d}" for i in range(600)],
         **{c: rng.negative_binomial(10, 10 / 30, 600).astype(np.int32) for c in names}}
    )  # fmt: skip
    counts = pl.concat([synthetic.counts, low]).with_columns(
        pl.when(pl.col("gene_id") == "g00500").then(pl.col(names[0]) * 100 + 5000)
        .otherwise(pl.col(names[0])).alias(names[0])
    )  # fmt: skip
    return fit(counts, synthetic.samples, BASIC, {}, threads=1).results


def random_p(seed: int = 3) -> np.ndarray:
    """Unsorted p-values with ties, p = 1, and some small enough to pass Bonferroni."""
    rng = np.random.default_rng(seed)
    p = np.concatenate([rng.uniform(size=400) ** 3, rng.uniform(0, 1e-5, 20), np.ones(15)])
    p = np.concatenate([p, p[:40], np.full(10, p[7])])
    return rng.permutation(p)


def statsmodels_adjust(p: np.ndarray, method: str) -> np.ndarray:
    from statsmodels.stats.multitest import multipletests  # noqa: PLC0415 - the oracle only

    return multipletests(p, alpha=ALPHA, method=STATSMODELS[method])[1]


def test_oracle_covers_every_method():
    assert list(STATSMODELS) == list(ADJUSTMENTS)


# --- against pydeseq2 and statsmodels ---------------------------------------------------------


def test_bh_equals_pydeseq2_padj(filtered):
    reasons = set(filtered["na_reason"])
    assert {"cooks_outlier", "independent_filter"} <= reasons, reasons
    tested = filtered.filter(pl.col("padj").is_not_null())
    ours = adjust(tested["pvalue"].to_numpy(), "padj")
    assert np.max(np.abs(ours - tested["padj"].to_numpy())) <= 1e-12
    # The set matters: BH over every gene with a p-value gives other numbers.
    every = filtered.filter(pl.col("pvalue").is_not_null())
    over_all = adjust(every["pvalue"].to_numpy(), "padj")[every["padj"].is_not_null().to_numpy()]
    assert np.max(np.abs(over_all - tested["padj"].to_numpy())) > 1e-3


@pytest.mark.parametrize("method", list(ADJUSTMENTS))
def test_fit_columns_over_tested_set(filtered, method):
    column = filtered[method]
    assert column.is_null().to_list() == filtered["padj"].is_null().to_list()
    tested = filtered.filter(pl.col("padj").is_not_null())
    oracle = statsmodels_adjust(tested["pvalue"].to_numpy(), method)
    assert np.max(np.abs(tested[method].to_numpy() - oracle)) <= 1e-12


@pytest.mark.parametrize("method", list(ADJUSTMENTS))
def test_adjust_equals_statsmodels(method):
    p = random_p()
    assert np.max(np.abs(adjust(p, method) - statsmodels_adjust(p, method))) <= 1e-12


@pytest.mark.parametrize("method", list(ADJUSTMENTS))
def test_adjust_monotone_in_p(method):
    p = np.sort(random_p())
    adjusted = adjust(p, method)
    assert np.all(np.diff(adjusted) >= 0)
    assert np.all(adjusted >= p) and np.all(adjusted <= 1)


def test_adjust_refuses_bad_input():
    with pytest.raises(ValueError, match="no NaN"):
        adjust(np.array([0.1, np.nan]), "padj_holm")
    with pytest.raises(ValueError, match="unknown adjustment"):
        adjust(np.array([0.1]), "fdr_bh")


def test_imports_without_polars():
    # Pyodide has no polars (TODO.md D4): with polars unreachable, the module still loads.
    code = (
        "import sys; sys.modules['polars'] = None\n"
        "import numpy as np\n"
        "from ctx.de import stats\n"
        "print(stats.adjust(np.array([0.01, 0.5]), 'padj_holm')[0], "
        "stats.p_line('padj_holm', stats.ALPHA, 2, 1))\n"
    )
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "0.02 0.025"


# --- the call rule at its boundaries ----------------------------------------------------------


def test_owner_thresholds():
    assert (ALPHA, LFC_MIN) == (0.05, 1.0)


@pytest.mark.parametrize(
    ("padj", "lfc", "call"),
    [
        (ALPHA, 2.0, "ns"),
        (np.nextafter(ALPHA, 0), 2.0, "up"),
        (0.01, LFC_MIN, "up"),
        (0.01, LFC_MIN - EPS, "ns"),
        (0.01, -LFC_MIN, "down"),
        (0.01, -(LFC_MIN - EPS), "ns"),
        (0.01, 0.0, "ns"),
        (None, 3.0, "ns"),
        (float("nan"), 3.0, "ns"),
    ],
)
def test_classify_boundaries(padj, lfc, call):
    frame = pl.DataFrame(
        {"padj": [padj], "lfc": [lfc]}, schema={"padj": pl.Float64, "lfc": pl.Float64}
    )
    assert frame.select(classify(pl.col("padj"), pl.col("lfc"))).item() == call


def test_classify_takes_its_thresholds():
    frame = pl.DataFrame({"padj": [0.08, 0.08], "lfc": [0.5, -0.5]})
    called = frame.select(classify(pl.col("padj"), pl.col("lfc"), alpha=0.1, lfc_min=0.5))
    assert called.to_series().to_list() == ["up", "down"]


# --- the volcano's horizontal line ------------------------------------------------------------


BASES = (*ADJUSTMENTS, "pvalue")


def passing(p: np.ndarray, basis: str, alpha: float) -> np.ndarray:
    """Which genes pass the cutoff on one basis: the rule, from the adjustment itself."""
    return (p if basis == "pvalue" else adjust(p, basis)) < alpha


def line_for(p: np.ndarray, basis: str, alpha: float) -> float:
    line = p_line(basis, alpha, p.size, int(passing(p, basis, alpha).sum()))
    assert line is not None
    return line


@pytest.mark.parametrize("basis", BASES)
@pytest.mark.parametrize("alpha", [0.1, ALPHA, 0.01, 0.001])
def test_line_splits_passing_from_failing(basis, alpha):
    # Exact because every adjustment is monotone in p, over data with ties and p = 1.
    p = random_p()
    assert ((p < line_for(p, basis, alpha)) == passing(p, basis, alpha)).all()


@pytest.mark.parametrize("basis", BASES)
@pytest.mark.parametrize("alpha", [0.1, ALPHA, 0.001])
def test_line_is_the_bar_the_weakest_gene_cleared(basis, alpha):
    # Any value between the weakest passing p and the strongest failing one would split them;
    # the line is the bar itself: the weakest passing gene moved just under it still passes,
    # just over it fails. Ties removed, so moving one gene does not reorder its twin.
    p = np.unique(random_p())
    passes = passing(p, basis, alpha)
    line = line_for(p, basis, alpha)
    weakest = np.flatnonzero(passes)[np.argmax(p[passes])]
    for factor, passes_after in ((1 - 1e-9, True), (1 + 1e-9, False)):
        moved = p.copy()
        moved[weakest] = line * factor
        assert passing(moved, basis, alpha)[weakest] == passes_after, factor


def test_line_moves_with_the_cutoff():
    # The reported fault: no BH-adjusted p lies between 0.02 and 0.05, so both call the same
    # two genes, and a line at the weakest passing gene's p did not move between them.
    p = np.array([1e-6, 2e-6, 0.3, 0.6, 0.9])
    assert passing(p, "padj", 0.05).sum() == passing(p, "padj", 0.02).sum() == 2
    assert line_for(p, "padj", 0.05) == pytest.approx(2.5 * line_for(p, "padj", 0.02))


@pytest.mark.parametrize("basis", BASES)
def test_line_when_none_pass(basis):
    # The first bar: what the strongest gene would have had to clear. No gene is under it.
    p = np.array([0.2, 0.5, 0.5, 1.0])
    assert not passing(p, basis, ALPHA).any()
    line = line_for(p, basis, ALPHA)
    first = {"pvalue": ALPHA, "padj_by": ALPHA / (4 * (1 + 1 / 2 + 1 / 3 + 1 / 4))}
    assert line == pytest.approx(first.get(basis, ALPHA / 4))
    assert not (p < line).any()


def test_no_line_without_tested_genes():
    assert p_line("padj", ALPHA, 0, 0) is None
    with pytest.raises(ValueError, match="unknown basis"):
        p_line("fdr_bh", ALPHA, 10, 1)
