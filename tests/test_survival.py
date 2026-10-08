"""ctx.clinical.survival, score and cutpoint against R and statsmodels, and their page twin
ctxSurv (templates/_surv.js.j2) against them (GS plan, tests 3 to 6). The R oracle is
tests/fixtures/survival/fixture.json, written by make_fixture.R under R 4.4.2 and survival
3.7.0; each test names the mutation it catches."""

import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import SurvfuncRight, survdiff

from ctx.clinical import ClinicalError, cutpoint, score, survival

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = json.loads((ROOT / "tests/fixtures/survival/fixture.json").read_text())
SCRIPT = ROOT / "src/ctx/report/templates/_surv.js.j2"
Z95 = 1.959963984540054


def column(rows: list[dict[str, Any]], name: str, kind: type = float) -> np.ndarray:
    """One column of the fixture's cohort."""
    return np.array([r[name] for r in rows], dtype=kind)


def nan(values: list[Any]) -> np.ndarray:
    """R's NA, which jsonlite writes as the text NA, as NaN."""
    return np.array([np.nan if v == "NA" else v for v in values], float)


A = FIXTURE["a"]
T, E, X1, X2 = column(A, "time"), column(A, "event", int), column(A, "x1"), column(A, "x2")
G2, G3 = column(A, "g2", str), column(A, "g3", str)


def test_kaplan_meier_median_rmst_and_follow_up_equal_rs():
    """Catches: Greenwood bands on the plain scale (R's log-log), the median taken at the first
    survival below one half rather than at or below it, its interval from the wrong band,
    RMST without the first step from 1, and follow-up as the median time."""
    curve, r = survival.kaplan_meier(T, E), FIXTURE["km"]
    assert curve.time.tolist() == r["time"]
    assert curve.at_risk.tolist() == r["n_risk"] and curve.events.tolist() == r["n_event"]
    np.testing.assert_allclose(curve.survival, r["surv"], rtol=1e-12)
    np.testing.assert_allclose(curve.lower, nan(r["lower"]), rtol=1e-10)
    np.testing.assert_allclose(curve.upper, nan(r["upper"]), rtol=1e-10)
    assert survival.median(curve) == (r["median"], r["median_lower"], r["median_upper"])
    assert survival.rmst(curve, FIXTURE["horizon"]) == pytest.approx(r["rmean"], rel=1e-12)
    assert survival.follow_up(T, E) == r["reverse_median"]
    by = FIXTURE["km_by_g2"]
    for i, level in enumerate(("high", "low")):
        got = survival.median(survival.kaplan_meier(T[level == G2], E[level == G2]))
        assert got == (by["median"][i][0], by["median_lower"][i][0], by["median_upper"][i][0])


def test_a_median_on_a_flat_stretch_at_one_half_is_its_midpoint():
    """Six patients whose estimate sits at exactly 0.5 from time 3 to 5: R 4.4.2 with
    conf.type "log-log" gives 4, the midpoint, with interval 1 to NA (run 2026-10-08).
    Catches the median taken where the estimate first falls below one half (5), and the flat
    stretch read from a censored row."""
    curve = survival.kaplan_meier(np.array([1, 2, 3, 4, 5, 6.0]), np.array([1, 1, 1, 0, 1, 0]))
    found = survival.median(curve)
    assert found[:2] == (4.0, 1.0) and np.isnan(found[2])


def test_kaplan_meier_equals_statsmodels():
    """statsmodels' SurvfuncRight as the second oracle of the estimate itself."""
    fit = SurvfuncRight(T, E)
    curve = survival.kaplan_meier(T, E)
    events = curve.events > 0
    np.testing.assert_allclose(curve.survival[events], fit.surv_prob, rtol=1e-12)


@pytest.mark.parametrize(("groups", "key"), [(G2, "logrank_g2"), (G3, "logrank_g3")])
def test_logrank_equals_r_and_statsmodels(groups, key):
    """Catches: the variance without the (n - d) / (n - 1) factor for tied events."""
    chi, df, p = survival.logrank(T, E, groups)
    assert (chi, df) == (pytest.approx(FIXTURE[key]["chisq"], rel=1e-12), FIXTURE[key]["df"])
    assert p == pytest.approx(FIXTURE[key]["p"], rel=1e-10)
    assert chi == pytest.approx(survdiff(T, E, groups)[0], rel=1e-12)


def test_the_split_statistic_is_the_logrank_of_two_groups():
    """survival.Splits, the scan's vectorised form, equals logrank on each split. Catches a
    second definition drifting from the first."""
    cuts = np.percentile(X1, [20, 35, 50, 65, 80])
    high = X1[:, None] >= cuts[None, :]
    expected = [survival.logrank(T, E, h)[0] for h in high.T]
    np.testing.assert_allclose(survival.Splits(T, E)(high), expected, rtol=1e-12)


@pytest.mark.parametrize(("key", "design", "breslow"), [
    ("cox_x1_efron", (X1,), False), ("cox_x1x2_efron", (X1, X2), False),
    ("cox_x1x2_breslow", (X1, X2), True)])  # fmt: skip
def test_cox_ph_test_and_concordance_equal_rs(key, design, breslow):
    """Catches: Breslow in place of Efron, the PH test ranking event times only (the
    preview's definition) or weighting each patient by its own time, and concordance that
    counts tied event times."""
    x, r = np.column_stack(design), FIXTURE[key]
    fit = survival.cox(T, E, x, breslow=breslow)
    np.testing.assert_allclose(fit.beta[0], r["coef"], rtol=1e-8)
    np.testing.assert_allclose(np.sqrt(np.diag(fit.cov[0])), r["se"], rtol=1e-7)
    np.testing.assert_allclose(fit.loglik[0], r["loglik"], rtol=1e-12)
    assert not fit.refused[0]
    tests, total, p, p_total = survival.cox_zph(T, E, x, fit.beta[0], breslow=breslow)
    # R's table: one row per term, then GLOBAL, which equals the term for one covariate.
    # Both fits stop when the log likelihood moves by under 1e-9 of itself, so beta agrees to
    # about 1e-9 and a small statistic, here 0.0006, to about 2e-7.
    np.testing.assert_allclose([*tests, total], r["zph_chisq"], rtol=1e-6)
    np.testing.assert_allclose([*p, p_total], r["zph_p"], rtol=1e-6)
    assert survival.concordance(T, E, x @ fit.beta[0]) == pytest.approx(r["concordance"], rel=1e-12)


def test_cox_equals_statsmodels_phreg_with_efron():
    """statsmodels 0.15.0's PHReg, Efron ties, on the fixture and on a larger tied cohort."""
    rng = np.random.default_rng(3)
    big_t = rng.integers(1, 60, 600).astype(float)
    big_e = (rng.random(600) < 0.4).astype(int)
    big_x = rng.normal(size=(600, 3))
    for t, e, x in ((T, E, np.column_stack([X1, X2])), (big_t, big_e, big_x)):
        ref = PHReg(t, x, status=e, ties="efron").fit()
        fit = survival.cox(t, e, x)
        np.testing.assert_allclose(fit.beta[0], ref.params, rtol=1e-7, atol=1e-10)
        np.testing.assert_allclose(np.sqrt(np.diag(fit.cov[0])), ref.bse, rtol=1e-6)


def test_a_stratified_fit_and_its_ph_test_equal_rs():
    """Catches: one baseline over the strata, and ranks computed within a stratum."""
    fit = survival.cox(T, E, X1[:, None], strata=G2)
    r = FIXTURE["cox_x1_strata"]
    assert fit.beta[0, 0] == pytest.approx(r["coef"], rel=1e-8)
    tests, _, _, _ = survival.cox_zph(T, E, X1[:, None], fit.beta[0], strata=G2)
    assert tests[0] == pytest.approx(r["zph_chisq"][0], rel=1e-7)


def test_a_monotone_likelihood_is_refused():
    """One non-zero value, on a censored patient, shaped like ENSG00000248769 (GS plan): R
    stops at beta = -8.68 and warns the coefficient may be infinite. Catches the check
    removed, when the fit reports a finite, converged hazard ratio."""
    b, r = FIXTURE["b"], FIXTURE["monotone"]
    fit = survival.cox(column(b, "time"), column(b, "event", int), column(b, "x")[:, None])
    assert "may be infinite" in r["warning"]
    assert fit.beta[0, 0] == pytest.approx(r["coef"], rel=1e-6)
    assert fit.refused.tolist() == [True]
    assert np.isnan(survival.wald(fit).hr).all()


def test_a_batch_of_fits_equals_each_fit_alone():
    """The per-gene screen fits genes as one batch. Catches a batch that shares a step size,
    a halving or a stop between its members."""
    rng = np.random.default_rng(5)
    genes = rng.normal(size=(6, len(T)))
    genes[3] = 0.0
    genes[3, 1] = 1.5  # a monotone member must not stop or bend the others
    batch = survival.cox(T, E, np.stack([np.column_stack([g, X2]) for g in genes]))
    for i, g in enumerate(genes):
        alone = survival.cox(T, E, np.column_stack([g, X2]))
        # One stopping rule, so the two agree to the step it leaves: about 1e-8 here.
        np.testing.assert_allclose(batch.beta[i], alone.beta[0], rtol=1e-7)
        assert batch.refused[i] == alone.refused[0]


def test_power_follows_hsieh_and_lavori():
    """The plan's figures for 145 PFI events: 59.3 % at HR 1.2 per SD, 88.5 % at 1.3, alpha
    0.05. Catches a one-sided alpha or the log of the HR left out."""
    assert survival.power(145, 1.2, 0.05) == pytest.approx(59.3, abs=0.05)
    assert survival.power(145, 1.3, 0.05) == pytest.approx(88.5, abs=0.05)
    assert survival.power(145, 1 / 1.3, 0.05) == pytest.approx(88.5, abs=0.05)


def test_the_score_follows_its_order():
    """GS plan test 4, on numbers worked by hand: two genes over four patients, the second
    selected down. Catches z computed from the unrounded values, the raw log2 CPM averaged
    without z, the sign dropped, standardising before the mean, and the filter run on z."""
    values = np.array([[1.0004, 2.0, 3.0, 4.0], [8.0, 6.0, 6.0, 4.0], [0.1, 0.2, 0.1, 0.3]])
    kept, flat = score.expressed(score.stored(values), 1.0)
    assert kept.tolist() == [True, True, False] and flat.tolist() == [False, False, False]
    kept_values = score.stored(values)[kept]
    assert kept_values[0, 0] == 1000  # 1.0004 is stored as 1.000
    z = score.zscores(kept_values)
    got = score.score(z, np.array([1, -1]))
    a = np.array([1, 2, 3, 4]) - 2.5
    b = np.array([8, 6, 6, 4]) - 6.0
    mean = (a / a.std() - b / b.std()) / 2
    np.testing.assert_allclose(got, (mean - mean.mean()) / mean.std(), rtol=1e-12)
    with pytest.raises(ClinicalError):
        score.score(z[:0], np.array([]))


def test_splits_put_the_cut_in_the_higher_group():
    """Catches the median split by > rather than >=, and tertiles by numpy's nearest rank."""
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    assert score.split(values, "median").tolist() == [0, 0, 0, 1, 1, 1]
    assert score.split(np.array([1.0, 2.0, 3.0]), "median").tolist() == [0, 1, 1]
    assert score.split(values, "tertiles").tolist() == [0, 0, 1, 1, 2, 2]


HARNESS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const S = ctxSurv, z = input.z;
const out = {};
for (const [name, call] of Object.entries(input.calls)) {
  const [fn, ...args] = call;
  out[name] = ({
    km: (t, e) => S.kaplanMeier(t, e, z),
    median: (t, e) => S.median(S.kaplanMeier(t, e, z)),
    rmst: (t, e, h) => S.rmst(S.kaplanMeier(t, e, z), h),
    follow: (t, e) => S.followUp(t, e, z),
    logrank: (t, e, g, k) => S.logrank(t, e, g, k),
    cox: (t, e, x, s) => S.cox(t, e, x, s),
    zph: (t, e, x, s) => { const f = S.cox(t, e, x, s); return S.zph(t, e, x, f.beta, s); },
    c: (t, e, x) => {
      const f = S.cox(t, e, x);
      return S.concordance(t, e, x.map(r => r.reduce((a, v, i) => a + v * f.beta[i], 0)));
    },
    wald: (t, e, x) => S.wald(S.cox(t, e, x), 0, z),
    percentile: (v, qs) => qs.map(q => S.percentile(v, q)),
    draws: (seed, n) => {
      const d = S.mulberry32(seed);
      return Array.from({length: n}, () => d());
    },
    shuffle: (v, seed) => S.shuffle(v, S.mulberry32(seed)),
    scan: (t, e, v, p, seed) => S.scan(t, e, v, 20, 80, p, seed),
    score: (rows, signs) => S.score(rows, signs),
    split: (v, rule) => S.split(v, rule),
  })[fn](...args);
}
console.log(JSON.stringify(out, (k, v) => (typeof v === "number" && !isFinite(v) ? String(v) : v)));
"""


def node(calls: dict[str, list[Any]]) -> dict[str, Any]:
    """Run ctxSurv under node on the calls; NaN comes back as the text NaN."""
    binary = shutil.which("node")
    assert binary, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = re.sub(r"\{%-? (?:end)?raw -?%\}", "", SCRIPT.read_text())
    done = subprocess.run([binary, "-e", script + HARNESS],
                          input=json.dumps({"z": Z95, "calls": calls}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def numbers(values: Any) -> np.ndarray:
    """Node's numbers, its NaN text read back."""
    return np.array([np.nan if v == "NaN" else v for v in values], float)


def test_the_page_twin_equals_python():
    """ctxSurv computes what ctx.clinical computes, on the R fixture. Catches either language
    changing alone: a band, a tie rule, the PH test's ranks, the concordance's ties, the
    halving or the refusal."""
    t, e, x = T.tolist(), E.tolist(), np.column_stack([X1, X2])
    stage = np.where(G2 == "high", 1, 0)
    b = FIXTURE["b"]
    got = node({
        "km": ["km", t, e], "median": ["median", t, e], "rmst": ["rmst", t, e, 1000],
        "follow": ["follow", t, e], "lr2": ["logrank", t, e, stage.tolist(), 2],
        "lr3": ["logrank", t, e, np.unique(G3, return_inverse=True)[1].tolist(), 3],
        "cox": ["cox", t, e, x.tolist(), None], "zph": ["zph", t, e, x.tolist(), None],
        "strata": ["zph", t, e, X1[:, None].tolist(), G2.tolist()],
        "c": ["c", t, e, x.tolist()], "wald": ["wald", t, e, x.tolist()],
        "monotone": ["cox", column(b, "time").tolist(), column(b, "event", int).tolist(),
                     [[v] for v in column(b, "x")], None],
    })  # fmt: skip
    curve = survival.kaplan_meier(T, E)
    for name in ("survival", "lower", "upper"):
        np.testing.assert_allclose(numbers(got["km"][name]), getattr(curve, name), rtol=1e-12)
    assert tuple(got["median"]) == survival.median(curve)
    assert got["rmst"] == pytest.approx(survival.rmst(curve, 1000), rel=1e-12)
    assert got["follow"] == survival.follow_up(T, E)
    assert got["lr2"]["chi2"] == pytest.approx(survival.logrank(T, E, G2)[0], rel=1e-12)
    assert got["lr3"]["chi2"] == pytest.approx(survival.logrank(T, E, G3)[0], rel=1e-12)
    assert got["lr2"]["p"] == pytest.approx(survival.logrank(T, E, G2)[2], rel=1e-6)
    fit = survival.cox(T, E, x)
    # Each language stops when its log likelihood moves by under 1e-9 of itself; sums taken
    # in another order leave beta agreeing to about 3e-9 here.
    np.testing.assert_allclose(got["cox"]["beta"], fit.beta[0], rtol=1e-7)
    np.testing.assert_allclose(got["cox"]["cov"], fit.cov[0], rtol=1e-6)
    tests, total, _, _ = survival.cox_zph(T, E, x, fit.beta[0])
    np.testing.assert_allclose([*got["zph"]["chi2"], got["zph"]["global"]], [*tests, total],
                               rtol=1e-6)  # fmt: skip
    strata = survival.cox(T, E, X1[:, None], strata=G2)
    assert got["strata"]["chi2"][0] == pytest.approx(
        survival.cox_zph(T, E, X1[:, None], strata.beta[0], strata=G2)[0][0], rel=1e-6
    )
    assert got["c"] == pytest.approx(survival.concordance(T, E, x @ fit.beta[0]), rel=1e-12)
    w = survival.wald(fit)
    assert got["wald"]["hr"] == pytest.approx(w.hr[0], rel=1e-7)
    assert got["wald"]["p"] == pytest.approx(w.p[0], rel=1e-6)
    assert got["monotone"]["refused"] is True


def test_percentiles_draws_and_shuffles_agree_exactly():
    """GS plan test 6: the same generator, seed and shuffle order in both languages, and the
    same percentile, bit for bit. Catches another generator, another seed, a shuffle from
    the front, and a percentile by another interpolation."""
    rng = np.random.default_rng(11)
    values = rng.normal(size=257).round(6)
    qs = [20, 33.333333333333336, 50, 66.66666666666667, 80, *range(21, 80, 7)]
    draw = survival.mulberry32(42)
    got = node({"pct": ["percentile", values.tolist(), qs], "draws": ["draws", 42, 5],
                "shuffle": ["shuffle", list(range(10)), 42]})  # fmt: skip
    assert got["pct"] == np.percentile(values, qs).tolist()
    assert got["draws"] == [draw() for _ in range(5)]
    assert got["shuffle"] == survival.mulberry32(42).shuffle(np.arange(10)).tolist()


def test_the_scan_gives_the_same_corrected_p_in_both_languages():
    """Catches an off-by-one scan bound (19 or 81), and a permutation drawn differently."""
    got = node({"scan": ["scan", T.tolist(), E.tolist(), X1.tolist(), 200, 42]})["scan"]
    found = cutpoint.scan(T, E, X1, permutations=200, seed=42)
    assert got["cuts"] == list(range(20, 81)) == found.cuts.tolist()
    assert got["beaten"] == found.beaten and got["corrected"] == found.corrected
    assert got["best"] == found.best
    np.testing.assert_allclose(got["statistic"], found.statistic, rtol=1e-10)
    np.testing.assert_allclose(got["maxima"], found.maxima, rtol=1e-10)
    assert 0 < found.corrected < 1


def test_the_page_score_and_split_equal_pythons():
    """ctxSurv.score from the stored int16 rows and ctxSurv.split, against score.py."""
    rng = np.random.default_rng(2)
    rows = score.stored(rng.normal(5, 2, size=(5, 40)))
    signs = np.array([1, -1, 1, 1, -1])
    got = node({"score": ["score", rows.tolist(), signs.tolist()],
                "median": ["split", X1.tolist(), "median"],
                "tertiles": ["split", X1.tolist(), "tertiles"]})  # fmt: skip
    np.testing.assert_allclose(got["score"], score.score(score.zscores(rows), signs), rtol=1e-10)
    assert got["median"] == score.split(X1, "median").tolist()
    assert got["tertiles"] == score.split(X1, "tertiles").tolist()
