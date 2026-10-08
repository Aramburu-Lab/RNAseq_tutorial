"""The M10 page's analyses as tables.

Survival by group, Cox models, the gene screen, expression by clinical variable and
confounding. Page words are ``page.py``'s; the statistics are ``survival``'s, one each.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu, norm, spearmanr
from scipy.stats import false_discovery_control as bh

from ctx.clinical import harmonise, survival

# Fewest patients a level needs to enter a Kruskal-Wallis test or a Cramér's V table.
MIN_GROUP = 10
# Age groups of the confounding map, in years: below the first edge, between, at or above.
AGE_EDGES = (50, 60, 70)
ALPHA = 0.05
TARGET_POWER = 80.0


@dataclass(frozen=True)
class Outcome:
    """One endpoint's analysis rows: the patients with a flag and a time, in table order."""

    patients: pl.DataFrame
    months: np.ndarray
    event: np.ndarray


def outcome(patients: pl.DataFrame, endpoint: str) -> Outcome:
    """The endpoint's analysis rows, times in months."""
    rows = harmonise.analysed(patients, endpoint)
    return Outcome(rows, rows[f"{endpoint}_months"].to_numpy().astype(float),
                   rows[endpoint].to_numpy().astype(int))  # fmt: skip


def km_frame(months: np.ndarray, event: np.ndarray, labels: np.ndarray,
             names: Sequence[str]) -> pl.DataFrame:  # fmt: skip
    """Each group's Kaplan-Meier estimate, long, in ``KmData``'s columns; names in order."""
    frames = []
    for name in names:
        member = labels == name
        if not member.any():
            continue
        c = survival.kaplan_meier(months[member], event[member])
        frames.append(pl.DataFrame({
            "group": name, "time": c.time, "at_risk": c.at_risk, "events": c.events,
            "censored": c.censored, "survival": c.survival, "lower": c.lower, "upper": c.upper,
        }))  # fmt: skip
    return pl.concat(frames).with_columns(pl.col("lower", "upper").fill_nan(None))


def groups(months: np.ndarray, event: np.ndarray, labels: np.ndarray, names: Sequence[str],
           horizon: float) -> pl.DataFrame:  # fmt: skip
    """Per group: patients, events, median with its interval, and RMST to the horizon."""
    rows = []
    for name in names:
        member = labels == name
        if not member.any():
            continue
        c = survival.kaplan_meier(months[member], event[member])
        mid, low, high = survival.median(c)
        rows.append({"group": name, "patients": int(member.sum()),
                     "events": int(event[member].sum()), "median": mid, "low": low, "high": high,
                     "rmst": survival.rmst(c, horizon)})  # fmt: skip
    return pl.DataFrame(rows).with_columns(pl.col("median", "low", "high").fill_nan(None))


@dataclass(frozen=True)
class Model:
    """One Cox model's first term: hazard ratio, interval, p, PH test p and C; None if refused."""

    patients: int
    events: int
    hr: float | None
    low: float | None
    high: float | None
    p: float | None
    ph_p: float | None
    c_index: float | None


def cox_model(months: np.ndarray, event: np.ndarray, x: np.ndarray,
              strata: np.ndarray | None = None) -> Model:  # fmt: skip
    """Fit one model; its first column's numbers, or None for each where it is refused."""
    fit = survival.cox(months, event, x, strata)
    n, events = len(months), int(event.sum())
    if fit.refused[0]:
        return Model(n, events, None, None, None, None, None, None)
    w = survival.wald(fit)
    _, _, p, _ = survival.cox_zph(months, event, x, fit.beta[0], strata)
    c = survival.concordance(months, event, x @ fit.beta[0])
    return Model(n, events, float(w.hr[0]), float(w.low[0]), float(w.high[0]), float(w.p[0]),
                 float(p[0]), c)  # fmt: skip


def covariates(rows: pl.DataFrame, names: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    """The covariates as columns, and which rows have every one."""
    if not names:
        return np.empty((rows.height, 0)), np.ones(rows.height, bool)
    values = rows.select(names).to_numpy().astype(float)
    return values, ~np.isnan(values).any(axis=1)


def screen(out: Outcome, z: np.ndarray, extra: np.ndarray, complete: np.ndarray,
           chunk: int = 500) -> pl.DataFrame:  # fmt: skip
    """Cox per gene, alone and with the covariates, z per gene over each model's patients.

    Args:
        out: The endpoint's analysis rows.
        z: Stored log2 CPM of the genes over those rows, genes by patients, as floats.
        extra: The covariates over the rows.
        complete: Rows with every covariate.
        chunk: Genes fitted per batch.

    Returns:
        Per gene: hr, low, high and p of each model, its q by Benjamini-Hochberg over the
        genes, and whether it was refused.
    """

    def fit(values: np.ndarray, rows: np.ndarray, others: np.ndarray) -> dict[str, np.ndarray]:
        """One model over every gene, in chunks."""
        sd = values[:, rows].std(axis=1, keepdims=True)
        scaled = (values[:, rows] - values[:, rows].mean(axis=1, keepdims=True)) / sd
        parts = []
        for start in range(0, len(scaled), chunk):
            block = scaled[start : start + chunk]
            design = np.concatenate(
                [block[:, :, None], np.broadcast_to(others, (len(block), *others.shape))], axis=2
            )
            parts.append(survival.wald(survival.cox(out.months[rows], out.event[rows], design)))
        return {
            k: np.concatenate([getattr(w, k) for w in parts]) for k in ("hr", "low", "high", "p")
        }

    alone = fit(z, np.ones(len(out.months), bool), np.empty((len(out.months), 0)))
    adjusted = fit(z, complete, extra[complete])
    frame = {f"{k}_{m}": v for m, model in (("alone", alone), ("adjusted", adjusted))
             for k, v in model.items()}  # fmt: skip
    table = pl.DataFrame(frame).with_columns(pl.all().fill_nan(None))
    for model in ("alone", "adjusted"):
        p = table[f"p_{model}"].to_numpy()
        q = np.full(len(p), np.nan)
        ok = ~np.isnan(p.astype(float))
        q[ok] = bh(p[ok].astype(float))
        table = table.with_columns(pl.Series(f"q_{model}", q).fill_nan(None))
    return table


def _test(values: np.ndarray, labels: np.ndarray, kind: str) -> tuple[int, float, float]:
    """(n, effect, p) of one test of expression against a variable.

    ``mann_whitney``: labels positive and negative, effect AUC of positive over negative.
    ``kruskal``: levels with at least ``MIN_GROUP`` patients, effect epsilon squared.
    ``spearman``: a numeric variable, effect rho.
    """
    if kind == "spearman":
        keep = ~np.isnan(labels.astype(float))
        rho, p = spearmanr(values[keep], labels[keep].astype(float))
        return int(keep.sum()), float(rho), float(p)
    if kind == "mann_whitney":
        a, b = values[labels == "positive"], values[labels == "negative"]
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        return len(a) + len(b), float(u / (len(a) * len(b))), float(p)
    present = np.array([v is not None for v in labels])
    levels = [lv for lv in np.unique(labels[present]) if (labels == lv).sum() >= MIN_GROUP]
    parts = [values[labels == lv] for lv in levels]
    h, p = kruskal(*parts)
    n = sum(len(x) for x in parts)
    return n, float(h / (n - 1)), float(p)


def expression_tests(patients: pl.DataFrame, values: Mapping[str, np.ndarray]) -> pl.DataFrame:
    """Each build gene against ER, stage, histology, menopause and age; q over all of them.

    Args:
        patients: The patient table.
        values: log2 CPM per build gene symbol, over the patients in table order.

    Returns:
        gene, variable, test, n, effect, p, q.
    """
    plans = (("er", "mann_whitney"), ("stage", "kruskal"), ("stage", "spearman"),
             ("histology", "kruskal"), ("menopause", "kruskal"), ("age", "spearman"))  # fmt: skip
    rows = []
    for gene, v in values.items():
        for variable, kind in plans:
            labels = patients[variable].to_numpy()
            if kind == "spearman":
                labels = np.array([np.nan if x is None else x for x in labels], float)
            n, effect, p = _test(np.asarray(v, float), labels, kind)
            rows.append({"gene": gene, "variable": variable, "test": kind, "n": n,
                         "effect": effect, "p": p})  # fmt: skip
    table = pl.DataFrame(rows)
    return table.with_columns(q=pl.Series(bh(table["p"].to_numpy())))


def cramers_v(sheet: pl.DataFrame, a: str, b: str) -> tuple[float | None, int]:
    """Cramér's V of two columns over rows with both, uncorrected; None if one is constant."""
    both = sheet.select(a, b).drop_nulls().filter((pl.col(a) != "") & (pl.col(b) != ""))
    counts = both.group_by(a, b).len().pivot(on=b, index=a, values="len").fill_null(0)
    table = counts.drop(a).to_numpy()
    if min(table.shape) < 2:  # noqa: PLR2004 - a constant column has no association
        return None, both.height
    statistic = chi2_contingency(table, correction=False).statistic
    return float(np.sqrt(statistic / (both.height * (min(table.shape) - 1)))), both.height


def age_group(age: pl.Expr) -> pl.Expr:
    """Age in ``AGE_EDGES`` groups, as text: "<50", "50-59", "60-69", ">=70"."""
    edges = AGE_EDGES
    expr: Any = pl.when(age < edges[0]).then(pl.lit(f"<{edges[0]}"))
    for low, high in zip(edges, edges[1:], strict=False):
        expr = expr.when(age < high).then(pl.lit(f"{low}-{high - 1}"))
    grouped: pl.Expr = expr.when(age.is_not_null()).then(pl.lit(f">={edges[-1]}"))
    return grouped


def detectable(events: int, alpha: float = ALPHA, target: float = TARGET_POWER) -> float:
    """The per-SD hazard ratio a Wald test detects with ``target`` % power at these events."""
    shift = norm.ppf(1 - alpha / 2) + norm.ppf(target / 100)
    return float(np.exp(shift / np.sqrt(events)))
