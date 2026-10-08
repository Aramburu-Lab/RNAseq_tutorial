"""Survival statistics in numpy, one definition each, mirrored by the page's ``ctxSurv``.

Kaplan-Meier with log(-log) Greenwood bands and R's median; the log-rank test with the
(n - d) / (n - 1) ties factor; restricted mean survival time; the Cox model with Efron ties,
R's ``coxph`` default, fitted by Newton's method with step halving as R fits it, and refused
where R would warn that a coefficient may be infinite; R's ``cox.zph(transform = "rank")``
test; Harrell's C as R's ``concordance``; Hsieh and Lavori's power. statsmodels 0.15.0 and
R 4.4.2 with survival 3.7.0 are the oracles (``tests/test_survival.py``,
``tests/fixtures/survival``). Times are in any unit; the page's are months.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.stats import chi2, norm, rankdata

LEVEL = 0.95
# R's coxph.control: eps on the relative change in log likelihood, 20 iterations, and the
# infinite-coefficient check's toler.inf = sqrt(eps).
EPS = 1e-9
ITERATIONS = 20
TOLER_INF = EPS**0.5
HALVINGS = 30
# R's quantile.survfit compares survival with p within sqrt(machine epsilon).
TOLERANCE = np.sqrt(np.finfo(float).eps)


@dataclass(frozen=True)
class Curve:
    """A Kaplan-Meier estimate at every time a patient has the event or is censored.

    Attributes:
        time: Distinct times, ascending.
        at_risk: Patients followed just before each time.
        events: Events at each time.
        censored: Patients censored at each time.
        survival: The estimate just after each time.
        lower: Lower bound of the pointwise interval, log(-log) Greenwood; NaN where the
            estimate is 0.
        upper: Upper bound, likewise.
    """

    time: np.ndarray
    at_risk: np.ndarray
    events: np.ndarray
    censored: np.ndarray
    survival: np.ndarray
    lower: np.ndarray
    upper: np.ndarray


def kaplan_meier(time: np.ndarray, event: np.ndarray, level: float = LEVEL) -> Curve:
    """The Kaplan-Meier estimate with R's ``conf.type = "log-log"`` interval.

    Args:
        time: Follow-up time per patient.
        event: 1 for an event, 0 for censored.
        level: The interval's coverage.

    Returns:
        The curve.
    """
    time, event = np.asarray(time, float), np.asarray(event, int)
    times, inverse = np.unique(time, return_inverse=True)
    events = np.bincount(inverse, weights=event, minlength=len(times))
    followed = np.bincount(inverse, minlength=len(times))
    at_risk = len(time) - np.concatenate([[0], np.cumsum(followed)[:-1]])
    survival = np.cumprod(1.0 - events / at_risk)
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(np.cumsum(events / (at_risk * (at_risk - events))))
        log_s = np.log(survival)
        spread = norm.ppf(0.5 + level / 2) * se / log_s
        upper = np.exp(-np.exp(np.log(-log_s) + spread))
        lower = np.exp(-np.exp(np.log(-log_s) - spread))
    edge = (survival == 1) | (survival == 0)
    fill = np.where(survival == 0, np.nan, 1.0)
    return Curve(times, at_risk, events.astype(int), (followed - events).astype(int), survival,
                 np.where(edge, fill, lower), np.where(edge, fill, upper))  # fmt: skip


def quantile(times: np.ndarray, survival: np.ndarray, p: float = 0.5) -> float:
    """The first time a curve falls to ``1 - p`` or below, as R's ``quantile.survfit``.

    Repeated values keep their first time; on a flat stretch exactly at the target the
    quantile is the midpoint to the next drop. NaN when the curve never reaches it.
    """
    target = 1 - p
    keep = ~np.isnan(survival)
    x, y = times[keep], survival[keep]
    first = np.concatenate([[True], y[1:] != y[:-1]])
    x, y = x[first], y[first]
    below = np.flatnonzero(y <= target + TOLERANCE)
    if not below.size:
        return float("nan")
    i = int(below[0])
    if i == len(y) - 1 or abs(y[i] - target) > TOLERANCE:
        return float(x[i])
    return float((x[i] + x[i + 1]) / 2)


def median(curve: Curve) -> tuple[float, float, float]:
    """Median survival and its interval, as R's ``quantile(fit, 0.5)``.

    Where the estimate falls to one half, then where its lower band does (the interval's low
    end), then its upper band (the high end).
    """
    return (quantile(curve.time, curve.survival), quantile(curve.time, curve.lower),
            quantile(curve.time, curve.upper))  # fmt: skip


def follow_up(time: np.ndarray, event: np.ndarray) -> float:
    """Median follow-up by reverse Kaplan-Meier: censoring as the event."""
    return median(kaplan_meier(time, 1 - np.asarray(event, int)))[0]


def rmst(curve: Curve, horizon: float) -> float:
    """Area under the step curve from 0 to the horizon, as R's ``rmean``."""
    kept = curve.time < horizon
    edges = np.concatenate([[0.0], curve.time[kept], [horizon]])
    heights = np.concatenate([[1.0], curve.survival[kept]])
    return float(np.sum(heights * np.diff(np.maximum(edges, 0.0))))


def logrank(time: np.ndarray, event: np.ndarray, groups: Sequence[str] | np.ndarray
            ) -> tuple[float, int, float]:  # fmt: skip
    """The log-rank test of k groups, as R's ``survdiff``: chi-squared, df and p.

    Each event time contributes observed minus expected events per group and the
    hypergeometric covariance with the (n - d) / (n - 1) factor for ties.
    """
    time, event = np.asarray(time, float), np.asarray(event, int)
    levels, member = np.unique(np.asarray(groups), return_inverse=True)
    one_hot = np.eye(len(levels))[member]
    times = np.unique(time[event == 1])
    n_at = (time[None, :] >= times[:, None]).astype(float) @ one_hot
    d_at = ((time[None, :] == times[:, None]) & (event[None, :] == 1)).astype(float) @ one_hot
    n, d = n_at.sum(axis=1), d_at.sum(axis=1)
    observed_minus = (d_at - d[:, None] * n_at / n[:, None]).sum(axis=0)
    ties = np.where(n > 1, (n - d) / np.maximum(n - 1, 1), 0.0)
    share = n_at / n[:, None]
    weight = (d * ties)[:, None, None]
    cov = (weight * (np.einsum("ti,ij->tij", share, np.eye(len(levels)))
                     - share[:, :, None] * share[:, None, :])).sum(axis=0)  # fmt: skip
    k = len(levels) - 1
    statistic = float(observed_minus[:k] @ np.linalg.solve(cov[:k, :k], observed_minus[:k]))
    return statistic, k, float(chi2.sf(statistic, k))


class Splits:
    """``logrank`` for two groups, vectorised over splits for the cut-point scan.

    The event times, risk sets and tied events are found once; each call takes patients by
    splits, True for high, and returns one chi-squared per split. The tests hold it equal to
    ``logrank``.
    """

    def __init__(self, time: np.ndarray, event: np.ndarray) -> None:
        """Sort the patients by time and find each event time's risk set and events."""
        self.order = np.argsort(time, kind="stable")
        t, e = np.asarray(time, float)[self.order], np.asarray(event, int)[self.order]
        times = np.unique(t[e == 1])
        self.at = np.searchsorted(t, times, side="left")
        self.risk = (len(t) - self.at).astype(float)[:, None]
        self.hits = ((t[None, :] == times[:, None]) & (e[None, :] == 1)).astype(float)
        self.deaths = self.hits.sum(axis=1)[:, None]
        self.ties = (self.risk - self.deaths) / np.maximum(self.risk - 1, 1)

    def __call__(self, high: np.ndarray) -> np.ndarray:
        """Chi-squared per column of ``high``, patients by splits in the patients' order."""
        g = np.asarray(high, float)[self.order]
        in_high = np.cumsum(g[::-1], axis=0)[::-1][self.at]
        share = in_high / self.risk
        excess = (self.hits @ g - self.deaths * share).sum(axis=0)
        variance = (self.deaths * share * (1 - share) * self.ties).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(np.where(variance > 0, excess**2 / variance, 0.0))


@dataclass(frozen=True)
class RiskSets:
    """Where each event time's risk set and tied events sit, patients sorted by stratum, time.

    Attributes:
        order: Patients in that order.
        start: Per event time, the first sorted position at risk.
        end: Per event time, one past its stratum's last position.
        tied: Per event time, the number of events.
        deaths: Sorted positions of the events, grouped by event time.
        groups: Per event time, where its events start in ``deaths``.
    """

    order: np.ndarray
    start: np.ndarray
    end: np.ndarray
    tied: np.ndarray
    deaths: np.ndarray
    groups: np.ndarray


def risk_sets(time: np.ndarray, event: np.ndarray, strata: np.ndarray | None = None) -> RiskSets:
    """The risk sets of every event time, within strata when given."""
    time, event = np.asarray(time, float), np.asarray(event, int)
    stratum = (
        np.zeros(len(time), int) if strata is None else np.unique(strata, return_inverse=True)[1]
    )
    order = np.lexsort((time, stratum))
    t, e, s = time[order], event[order], stratum[order]
    keys = np.stack([s, t], axis=1)[e == 1]
    distinct, first = np.unique(keys, axis=0, return_index=True)
    positions = np.flatnonzero(e == 1)
    start = np.array([np.flatnonzero((s == k[0]) & (t >= k[1]))[0] for k in distinct], dtype=int)
    end = np.array([np.flatnonzero(s == k[0])[-1] + 1 for k in distinct], dtype=int)
    groups = np.sort(first)
    tied = np.diff(np.concatenate([groups, [len(positions)]]))
    return RiskSets(order, start, end, tied, positions, groups)


def efron(sets: RiskSets, x: np.ndarray, eta: np.ndarray, breslow: bool = False
          ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:  # fmt: skip
    """Log partial likelihood, score and information at a linear predictor, Efron's ties.

    Args:
        sets: ``risk_sets`` of the patients.
        x: Covariates, batches by patients by covariates, in the patients' own order.
        eta: Linear predictor, batches by patients.
        breslow: Breslow's ties instead, for the oracle tests alone.

    Returns:
        Log likelihood (B), score (B, p), information (B, p, p).
    """
    loglik, score, info = efron_terms(sets, x, eta, breslow)
    return loglik, score.sum(axis=1), info.sum(axis=1)


def efron_terms(sets: RiskSets, x: np.ndarray, eta: np.ndarray, breslow: bool = False
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:  # fmt: skip
    """``efron`` with each event time's score and information apart.

    Shapes (B, J, p) and (B, J, p, p): the proportional-hazards test weights them by time.
    """
    xs, es = x[:, sets.order], eta[:, sets.order]
    w = np.exp(es - es.max(axis=1, keepdims=True))
    shift = es.max(axis=1)

    def tail(values: np.ndarray) -> np.ndarray:
        """Sums from each position to the end, a zero appended."""
        rev = np.cumsum(values[:, ::-1], axis=1)[:, ::-1]
        return np.concatenate([rev, np.zeros_like(rev[:, :1])], axis=1)

    wx = w[..., None] * xs
    wxx = wx[..., :, None] * xs[..., None, :]
    s0 = tail(w)[:, sets.start] - tail(w)[:, sets.end]
    s1 = tail(wx)[:, sets.start] - tail(wx)[:, sets.end]
    s2 = tail(wxx)[:, sets.start] - tail(wxx)[:, sets.end]
    d0 = np.add.reduceat(w[:, sets.deaths], sets.groups, axis=1)
    d1 = np.add.reduceat(wx[:, sets.deaths], sets.groups, axis=1)
    d2 = np.add.reduceat(wxx[:, sets.deaths], sets.groups, axis=1)
    loglik = es[:, sets.deaths].sum(axis=1) - (shift[:, None] * sets.tied[None, :]).sum(axis=1)
    score = np.add.reduceat(xs[:, sets.deaths], sets.groups, axis=1)
    info = np.zeros((*s2.shape,))
    for step in range(int(sets.tied.max(initial=0))):
        on = (step < sets.tied).astype(float)
        f = on * (0.0 if breslow else step / sets.tied)
        a0 = np.where(on[None, :] > 0, s0 - f * d0, 1.0)
        a1 = s1 - f[None, :, None] * d1
        a2 = s2 - f[None, :, None, None] * d2
        mean = a1 / a0[..., None]
        loglik -= (on[None, :] * np.log(a0)).sum(axis=1)
        score -= on[None, :, None] * mean
        cov = a2 / a0[..., None, None] - mean[..., :, None] * mean[..., None, :]
        info += on[None, :, None, None] * cov
    return loglik, score, info


@dataclass(frozen=True)
class CoxFit:
    """Cox fits of a batch of designs over the same patients.

    Attributes:
        beta: Coefficients, batches by covariates.
        cov: Inverse information at ``beta``.
        loglik: Log partial likelihood at 0 and at ``beta``, batches by 2.
        iterations: Newton steps taken per batch.
        infinite: Per coefficient, R's flag "coefficient may be infinite": a monotone
            likelihood, which ``wald`` refuses to report.
    """

    beta: np.ndarray
    cov: np.ndarray
    loglik: np.ndarray
    iterations: np.ndarray
    infinite: np.ndarray

    @property
    def refused(self) -> np.ndarray:
        """Per batch, whether any coefficient may be infinite."""
        return np.asarray(self.infinite.any(axis=1))


def cox(time: np.ndarray, event: np.ndarray, x: np.ndarray, strata: np.ndarray | None = None,
        breslow: bool = False) -> CoxFit:  # fmt: skip
    """Fit Cox models by Newton-Raphson from 0, as R's ``coxph.fit``.

    A step that lowers the likelihood is halved; the fit stops when the log likelihood
    changes by less than ``EPS`` relative to it, or after ``ITERATIONS`` steps. A coefficient
    is flagged infinite as R flags it: the Newton step still exceeds ``EPS`` and
    ``TOLER_INF`` times the coefficient.

    Args:
        time: Follow-up time per patient.
        event: 1 for an event, 0 for censored.
        x: Covariates, patients by covariates, or batches by patients by covariates.
        strata: A label per patient; each stratum has its own baseline hazard.
        breslow: Breslow's ties, for the oracle tests alone.

    Returns:
        The fits.
    """
    x = np.asarray(x, float)
    x = x[None] if x.ndim == 2 else x  # noqa: PLR2004 - one design or a batch of them
    batches, _, p = x.shape
    sets = risk_sets(time, event, strata)
    beta = np.zeros((batches, p))
    loglik, score, info = efron(sets, x, np.zeros(x.shape[:2]), breslow)
    start = loglik.copy()
    steps = np.zeros(batches, int)
    done = np.zeros(batches, bool)
    for _ in range(ITERATIONS):
        trial = beta + np.where(done[:, None], 0.0, _solve(info, score))
        new = efron(sets, x, np.einsum("bnp,bp->bn", x, trial), breslow)
        for _ in range(HALVINGS):
            worse = (new[0] < loglik) & ~done
            if not worse.any():
                break
            trial = np.where(worse[:, None], (beta + trial) / 2, trial)
            new = efron(sets, x, np.einsum("bnp,bp->bn", x, trial), breslow)
        settled = np.abs(1 - loglik / new[0]) <= EPS
        steps += ~done
        beta = np.where(done[:, None], beta, trial)
        loglik, score, info = (np.where(done, loglik, new[0]),
                               np.where(done[:, None], score, new[1]),
                               np.where(done[:, None, None], info, new[2]))  # fmt: skip
        done |= settled
        if done.all():
            break
    cov = np.linalg.pinv(info)
    step = np.abs(np.einsum("bp,bpq->bq", score, cov))
    infinite = (step > EPS) & (step > TOLER_INF * np.abs(beta))
    return CoxFit(beta, cov, np.stack([start, loglik], axis=1), steps, infinite)


def _solve(info: np.ndarray, score: np.ndarray) -> np.ndarray:
    """Newton steps of a batch; a singular information takes its pseudo-inverse."""
    return np.asarray(np.einsum("bpq,bq->bp", np.linalg.pinv(info), score))


@dataclass(frozen=True)
class Wald:
    """One coefficient's hazard ratio, interval and two-sided Wald p; NaN where refused."""

    hr: np.ndarray
    low: np.ndarray
    high: np.ndarray
    p: np.ndarray
    se: np.ndarray


def wald(fit: CoxFit, column: int = 0, level: float = LEVEL) -> Wald:
    """The hazard ratio of one covariate per batch, NaN where the fit is refused."""
    beta = fit.beta[:, column]
    se = np.sqrt(fit.cov[:, column, column])
    z = norm.ppf(0.5 + level / 2)
    with np.errstate(over="ignore", invalid="ignore"):  # a refused fit's bounds; NaN below
        values = np.stack([np.exp(beta), np.exp(beta - z * se), np.exp(beta + z * se),
                           2 * norm.sf(np.abs(beta / se)), se])  # fmt: skip
    hr, low, high, p, se = np.where(fit.refused[None, :], np.nan, values)
    return Wald(hr, low, high, p, se)


def cox_zph(  # noqa: PLR0913 - the fit's own data, then its options
    time: np.ndarray, event: np.ndarray, x: np.ndarray, beta: np.ndarray,
    strata: np.ndarray | None = None, breslow: bool = False,
) -> tuple[np.ndarray, float, np.ndarray, float]:  # fmt: skip
    """R's ``cox.zph(fit, transform = "rank")``: the score test of beta(t) = beta + theta g(t).

    g is the rank of every patient's time, ties averaged, centred on its mean over the events,
    as R computes it; ranking the event times alone was the preview's definition (GS plan).
    At an event time every patient at risk carries x g(t), with g of that time: the score
    for theta weights each event time's score by g, and its information by g and g squared.

    Returns:
        Chi-squared and p per covariate, then the global chi-squared and p.
    """
    time, x = np.asarray(time, float), np.asarray(x, float)
    x = x[:, None] if x.ndim == 1 else x
    ranks = rankdata(time)
    g = ranks - ranks[np.asarray(event) == 1].mean()
    sets = risk_sets(time, event, strata)
    _, score, info = efron_terms(sets, x[None], (x @ beta)[None], breslow)
    at = g[sets.order][sets.deaths][sets.groups]
    p = x.shape[1]
    u = np.concatenate([np.zeros(p), np.einsum("j,jp->p", at, score[0])])
    cross = np.einsum("j,jpq->pq", at, info[0])
    full = np.block([[info[0].sum(axis=0), cross],
                     [cross, np.einsum("j,jpq->pq", at**2, info[0])]])  # fmt: skip
    tests = np.empty(p)
    for k in range(p):
        keep = [*range(p), p + k]
        sub = u[keep]
        tests[k] = sub @ np.linalg.solve(full[np.ix_(keep, keep)], sub)
    total = float(u @ np.linalg.solve(full, u))
    return tests, total, chi2.sf(tests, 1), float(chi2.sf(total, p))


def concordance(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    """Harrell's C as R's ``concordance``: of comparable pairs, the share ordered by risk.

    A pair is comparable when the shorter time is an event, or both times are equal and only
    one is an event; tied risks count one half; pairs of tied event times are left out.
    """
    time, event, risk = (np.asarray(v, float) for v in (time, event, risk))
    earlier = (time[:, None] < time[None, :]) | ((time[:, None] == time[None, :])
                                                 & (event[None, :] == 0))  # fmt: skip
    comparable = earlier & (event[:, None] == 1)
    np.fill_diagonal(comparable, False)
    higher = (risk[:, None] > risk[None, :]) + 0.5 * (risk[:, None] == risk[None, :])
    return float(higher[comparable].sum() / comparable.sum())


def power(events: np.ndarray | float, hr: np.ndarray | float, alpha: float) -> np.ndarray:
    """Power in % of the two-sided Wald test of a per-SD hazard ratio (Hsieh and Lavori 2000)."""
    shift = np.sqrt(events) * np.abs(np.log(hr)) - norm.ppf(1 - alpha / 2)
    return np.asarray(100.0 * norm.cdf(shift))


def mulberry32(seed: int) -> "Mulberry":
    """The page's random generator, so a permutation p is the same in Python and the page."""
    return Mulberry(seed)


class Mulberry:
    """mulberry32 in 32-bit unsigned arithmetic, as ``ctxSurv.mulberry32`` draws it."""

    MASK = 0xFFFFFFFF

    def __init__(self, seed: int) -> None:
        """Start from a seed, taken modulo 2 ** 32."""
        self.state = seed & self.MASK

    def __call__(self) -> float:
        """The next draw in [0, 1)."""
        m = self.MASK
        self.state = (self.state + 0x6D2B79F5) & m
        a = self.state
        t = ((a ^ (a >> 15)) * (1 | a)) & m
        t = ((t + (((t ^ (t >> 7)) * (61 | t)) & m)) & m) ^ t
        return ((t ^ (t >> 14)) & m) / 4294967296

    def shuffle(self, values: np.ndarray) -> np.ndarray:
        """A Fisher-Yates permutation of a copy, from the last element down."""
        out = values.copy()
        for i in range(len(out) - 1, 0, -1):
            j = int(self() * (i + 1))
            out[i], out[j] = out[j], out[i]
        return out
