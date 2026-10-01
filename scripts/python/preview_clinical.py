#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side preview builder; the M10 page on made-up data until ctx clinical exists
# ///
"""Course preview of M10 Clinical and survival, drawn on made-up data until the module is built.

The survival figures come from a made-up tumour cohort drawn here from a seeded generator:
expression of one gene, age, stage and a progression-free interval with censoring. lifelines is
not in the analysis environment, so Kaplan-Meier, the log-rank test and a Breslow Cox fit are
computed with numpy below. The confounding map is computed on the made-up dataset's own sample
sheet. The page text is ``content/clinical.toml``, which the built module will reuse. Run
through ``build_preview.py --only clinical``.
"""

from __future__ import annotations

import time
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import chi2, chi2_contingency, norm

from ctx import runtime
from ctx.de.run import Writer, newest
from ctx.io.naming import MODULES, stem
from ctx.report import text
from ctx.report.builder import (
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    write,
)
from ctx.viz.contour import Mark, SurfaceData
from ctx.viz.curve import CurveData
from ctx.viz.heatmap import HeatmapData, Swatch
from ctx.viz.scatter import Line
from ctx.viz.style import HEAT_ALPHAS, ROLES

# The made-up dataset's name on the page, one word for every preview (content/ui.toml).
SHOWN = ("ui", "preview.dataset")

PAGE = MODULE = "clinical"
SEED = 10
PATIENTS = 420
ENDPOINT = "PFI"
HORIZON = 60  # months: the RMST horizon
RISK_MONTHS = (0, 24, 48, 72, 96, 120)
COVARIATES = ("age", "stage")
SCAN = (20, 80)  # percentiles scanned, as the Human Protein Atlas and bc-GenExMiner scan
PERMUTATIONS = 1000
ALPHA = 0.05
TARGET_POWER = 80.0
LEVEL = 0.95
DAYS_PER_MONTH = 365.25 / 12
TICK = 0.015  # half-height of a censor tick, in survival fraction
MEDIAN, HALF = 50, 0.5
CONFOUNDERS = ("condition", "cell_type", "batch", "sex")
V_EDGES = (0.1, 0.3, 0.5)  # Cohen's small, medium and large for one degree of freedom


def fill(key: str, **values: Any) -> str:
    """A template of this page's content file, filled; ``text.fill`` on the page."""
    return text.fill(PAGE, key, **values)


def lines(key: str, **values: Any) -> tuple[str, ...]:
    """Every template of a list in the page file, filled as strictly as ``text.fill`` fills one."""
    filled = []
    for template in text.words(PAGE, key):
        if missing := text.placeholders(template) - values.keys():
            raise text.TextError(f"content/{PAGE}.toml {key!r} needs {sorted(missing)}")
        filled.append(template.format_map(values))
    return tuple(filled)


def listed(names: tuple[str, ...]) -> str:
    """Names as a sentence lists them: "age and stage"."""
    if len(names) == 1:
        return names[0]
    return fill("results.and", rest=", ".join(names[:-1]), last=names[-1])


def columns(what: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) per column of one table or figure, from the page file."""
    return {name: (unit, meaning) for name, (unit, meaning) in text.raw(PAGE, what).items()}


def split(values: np.ndarray, percentile: float) -> np.ndarray:
    """High is at or above the percentile: the one split rule of the curves and the scan."""
    return np.asarray(values >= np.percentile(values, percentile))


def made_up_cohort(rng: np.random.Generator) -> pl.DataFrame:
    """Patients with expression, age, stage and a censored time to progression, all made up.

    Expression rises with stage, so part of its hazard ratio is stage's, which the adjusted
    model shows. Follow-up runs to twelve years, as in TCGA, with some patients lost earlier.
    """
    stage = rng.choice(np.arange(1, 5), PATIENTS, p=(0.2, 0.45, 0.25, 0.1))
    age = np.clip(rng.normal(58.0, 12.0, PATIENTS), 26.0, 89.0).round()
    log2cpm = 6.0 + 0.4 * (stage - 2.2) + rng.normal(0.0, 1.0, PATIENTS)
    z = (log2cpm - log2cpm.mean()) / log2cpm.std()
    hazard = 0.0045 * np.exp(0.25 * z + 0.45 * (stage - 1) + 0.015 * (age - 58.0))
    event_days = rng.exponential(1.0 / hazard) * DAYS_PER_MONTH
    lost = rng.exponential(9000.0, PATIENTS)
    censor_days = np.minimum(rng.uniform(120.0, 4400.0, PATIENTS), lost)
    days = np.ceil(np.minimum(event_days, censor_days)).astype(int)
    group = np.where(split(z, MEDIAN), fill("split.high"), fill("split.low"))
    return pl.DataFrame({
        "patient": [f"P{i:03d}" for i in range(PATIENTS)], "age": age, "stage": stage,
        "log2cpm": log2cpm, "z": z, "group": group, "days": days,
        "months": days / DAYS_PER_MONTH, "event": (event_days <= censor_days).astype(int),
    })  # fmt: skip


def kaplan_meier(months: np.ndarray, event: np.ndarray) -> pl.DataFrame:
    """The estimate at every time a patient has the event or is censored, with its interval.

    The interval is Greenwood's variance on the log(-log) scale, as lifelines draws it.
    """
    times, first = np.unique(np.sort(months), return_index=True)
    at_risk = len(months) - first
    events = np.array([event[months == t].sum() for t in times])
    survival = np.cumprod(1.0 - events / at_risk)
    with np.errstate(divide="ignore", invalid="ignore"):
        greenwood = np.cumsum(events / (at_risk * (at_risk - events)))
        spread = np.where(greenwood > 0, norm.ppf(0.5 + LEVEL / 2) * np.sqrt(greenwood)
                          / np.log(survival), 0.0)  # fmt: skip
        bounds = [survival ** np.exp(sign * spread) for sign in (-1.0, 1.0)]
    lower, upper = (np.where(np.isfinite(b), b, survival) for b in bounds)
    followed = np.bincount(np.searchsorted(times, months), minlength=len(times))
    return pl.DataFrame({"months": times, "at_risk": at_risk, "events": events,
                         "censored": followed - events, "survival": survival, "lower": lower,
                         "upper": upper})  # fmt: skip


def first_below_half(km: pl.DataFrame, column: str) -> float:
    """The first time a curve reaches one half or below; NaN when it never does."""
    reached = km.filter(pl.col(column) <= HALF)["months"]
    return float(reached[0]) if reached.len() else float("nan")


def rmst(km: pl.DataFrame, horizon: float) -> float:
    """Area under the step curve from 0 to the horizon, in months."""
    kept = km.filter(pl.col("months") < horizon)
    edges = np.concatenate([[0.0], kept["months"].to_numpy(), [horizon]])
    heights = np.concatenate([[1.0], kept["survival"].to_numpy()])
    return float(np.sum(heights * np.diff(edges)))


def logrank(months: np.ndarray, event: np.ndarray, high: np.ndarray) -> np.ndarray:
    """Log-rank chi-squared of each column of ``high`` (patients by splits), ties handled."""
    order = np.argsort(months, kind="stable")
    t, e, g = months[order], event[order].astype(float), high[order].astype(float)
    _, first = np.unique(t, return_index=True)
    at_risk = (len(t) - first).astype(float)[:, None]
    share = np.cumsum(g[::-1], axis=0)[::-1][first] / at_risk
    deaths = np.add.reduceat(e, first)[:, None]
    deaths_high = np.add.reduceat(e[:, None] * g, first, axis=0)
    ties = (at_risk - deaths) / np.maximum(at_risk - 1.0, 1.0)
    excess = (deaths_high - deaths * share).sum(axis=0)
    return np.asarray(excess**2 / (deaths * share * (1.0 - share) * ties).sum(axis=0))


def cox(months: np.ndarray, event: np.ndarray, x: np.ndarray) -> dict[str, np.ndarray]:
    """Cox fit by Newton-Raphson on Breslow's partial likelihood.

    Returns:
        ``beta``, its covariance ``cov``, the Schoenfeld residuals ``resid`` (one row per
        event) and the event ``times``, both in time order.
    """
    order = np.argsort(months, kind="stable")
    t, e, xs = months[order], event[order].astype(bool), x[order]
    start = np.searchsorted(t, t, side="left")  # tied times share one risk set

    def at_risk(values: np.ndarray) -> np.ndarray:
        """Sum over each patient's risk set: everyone whose time is not earlier."""
        return np.cumsum(values[::-1], axis=0)[::-1][start]

    beta = np.zeros(xs.shape[1])
    for _ in range(50):
        w = np.exp(xs @ beta)
        s0 = at_risk(w)[:, None]
        mean = at_risk(w[:, None] * xs) / s0
        second = at_risk(w[:, None, None] * xs[:, :, None] * xs[:, None, :]) / s0[:, :, None]
        info = (second - mean[:, :, None] * mean[:, None, :])[e].sum(axis=0)
        step = np.linalg.solve(info, (xs - mean)[e].sum(axis=0))
        beta = beta + step
        if np.abs(step).max() < 1e-10:  # noqa: PLR2004 - Newton converges in a few steps
            break
    return {"beta": beta, "cov": np.linalg.inv(info), "resid": (xs - mean)[e], "times": t[e]}


def ph_test(fit: dict[str, np.ndarray]) -> np.ndarray:
    """Grambsch-Therneau p per covariate: scaled Schoenfeld residuals against ranked time."""
    rank = np.argsort(np.argsort(fit["times"], kind="stable")).astype(float)
    centred = rank - rank.mean()
    n_events = len(rank)
    scaled = fit["beta"] + n_events * fit["resid"] @ fit["cov"]
    statistic = (centred @ scaled) ** 2 / (n_events * np.diag(fit["cov"]) * (centred @ centred))
    return np.asarray(chi2.sf(statistic, 1))


def c_index(months: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    """Harrell's C: of pairs where the earlier time is an event, the share ordered by risk."""
    earlier = (months[:, None] < months[None, :]) & (event[:, None] == 1)
    higher = (risk[:, None] > risk[None, :]) + 0.5 * (risk[:, None] == risk[None, :])
    return float(higher[earlier].sum() / earlier.sum())


def power(events: np.ndarray, hr: np.ndarray) -> np.ndarray:
    """Power in % of the two-sided Wald test of a per-SD hazard ratio (Hsieh and Lavori 2000)."""
    shift = np.sqrt(events) * np.abs(np.log(hr)) - norm.ppf(1 - ALPHA / 2)
    return np.asarray(100.0 * norm.cdf(shift))


def km_path(km: pl.DataFrame, group: str) -> list[tuple[Any, ...]]:
    """Vertices of one group's step line, a censor tick drawn into it at each censored time.

    The curve contract draws lines and no markers, so a tick is a short vertical run of the line.
    """
    step, tick = fill("figures.km.step"), fill("figures.km.censored")
    rows: list[tuple[Any, ...]] = [(group, 0.0, 1.0, km["at_risk"][0], step)]
    previous = 1.0
    for t, n, d, c, s in km.select("months", "at_risk", "events", "censored", "survival").rows():
        rows.append((group, t, previous, n, step))
        if d:
            rows.append((group, t, s, n, step))
        if c:
            rows += [(group, t, s + TICK, n, tick), (group, t, s - TICK, n, tick),
                     (group, t, s, n, tick)]  # fmt: skip
        previous = s
    return rows


def group_summary(estimates: dict[str, pl.DataFrame], cohort: pl.DataFrame) -> ResultTable:
    """Patients, events, median with its interval, and RMST, per group."""
    rows = []
    for name, km in estimates.items():
        mid = [first_below_half(km, c) for c in ("survival", "lower", "upper")]
        shown = [fill("summary.not_reached") if np.isnan(v) else f"{v:.1f}" for v in mid]
        members = cohort.filter(pl.col("group") == name)
        rows.append((name, str(members.height), str(members["event"].sum()), shown[0],
                     fill("summary.interval", low=shown[1], high=shown[2]),
                     f"{rmst(km, HORIZON):.1f}"))  # fmt: skip
    heads = tuple(fill(f"summary.columns.{k}", level=100 * LEVEL)
                  for k in text.raw(PAGE, "summary.columns"))  # fmt: skip
    title = fill("summary.title", endpoint=ENDPOINT, horizon=HORIZON)
    return ResultTable(title, heads, tuple(rows), numeric=(1, 2, 3, 5))


def risk_table(cohort: pl.DataFrame) -> ResultTable:
    """Patients still followed at each of ``RISK_MONTHS``, per group."""
    rows = tuple(
        (name, *(str(members.filter(pl.col("months") >= t).height) for t in RISK_MONTHS))
        for name in (fill("split.low"), fill("split.high"))
        for members in (cohort.filter(pl.col("group") == name),)
    )
    heads = (fill("risk.group"), *(str(t) for t in RISK_MONTHS))
    return ResultTable(fill("risk.title"), heads, rows, numeric=tuple(range(1, len(heads))))


def survival_group(writer: Writer, cohort: pl.DataFrame, gene: str) -> tuple[Group, dict[str, Any]]:
    """Kaplan-Meier figure, the per-group summary and the risk table; the numbers they print."""
    months, event = cohort["months"].to_numpy(), cohort["event"].to_numpy()
    low, high = fill("split.low"), fill("split.high")
    estimates = {name: kaplan_meier(months[m], event[m]) for name in (low, high)
                 for m in (cohort["group"].to_numpy() == name,)}  # fmt: skip
    long = pl.concat([k.insert_column(0, pl.lit(n).alias("group")) for n, k in estimates.items()])
    writer.table(long, "km_estimates", columns("columns.km_estimates"))
    statistic = float(logrank(months, event, (cohort["group"] == high).to_numpy()[:, None])[0])
    medians = [first_below_half(estimates[n], "survival") for n in (low, high)]
    numbers = {
        "source": fill("preview.source"), "n": cohort.height, "events": int(event.sum()),
        "endpoint": ENDPOINT, "gene": gene, "chi2": statistic, "p": f"{chi2.sf(statistic, 1):.2g}",
        "horizon": HORIZON, "rmst": rmst(estimates[high], HORIZON) - rmst(estimates[low], HORIZON),
        **{f"median_{g}": fill("summary.not_reached") if np.isnan(m) else f"{m:.1f} months"
           for g, m in zip(("low", "high"), medians, strict=True)},
    }  # fmt: skip
    path = pl.DataFrame([r for n, k in estimates.items() for r in km_path(k, n)],
                        schema=list(columns("columns.km")), orient="row")  # fmt: skip
    data = CurveData(
        frame=path, columns=columns("columns.km"), series="group", x="months", y="survival",
        x_label=fill("figures.km.x"), y_label=fill("figures.km.y", **numbers), colour="group",
        colours={low: ROLES["reference"], high: ROLES["comparison"]},
        caption=fill("figures.km.caption", **numbers), hover=("at_risk", "mark"),
    )  # fmt: skip
    words = FigureText(title=fill("figures.km.title", **numbers),
                       lead=fill("figures.km.lead", **numbers),
                       reading=lines("figures.km.reading", **numbers))  # fmt: skip
    items = (writer.figure(data, "km", words), group_summary(estimates, cohort), risk_table(cohort))
    return Group(fill("results.survival", **numbers), items), numbers


def fit_row(name: str, cohort: pl.DataFrame, used: list[str]) -> dict[str, Any]:
    """One Cox fit: the expression term's hazard ratio, interval, p and PH test, and C."""
    months, event = cohort["months"].to_numpy(), cohort["event"].to_numpy()
    x = cohort.select(used).to_numpy().astype(float)
    fit = cox(months, event, x)
    beta, se = float(fit["beta"][0]), float(np.sqrt(fit["cov"][0, 0]))
    wald = norm.ppf(0.5 + LEVEL / 2) * se
    return {
        "model": name, "patients": cohort.height, "events": int(event.sum()),
        "hr": float(np.exp(beta)), "low": float(np.exp(beta - wald)),
        "high": float(np.exp(beta + wald)), "p": float(2 * norm.sf(abs(beta) / se)),
        "ph_p": float(ph_test(fit)[0]), "c_index": c_index(months, event, x @ fit["beta"]),
    }  # fmt: skip


def cox_group(writer: Writer, cohort: pl.DataFrame, gene: str) -> tuple[Group, dict[str, Any]]:
    """Univariable and adjusted Cox fits, as a table; the adjusted fit for the page's numbers."""
    covariates = listed(COVARIATES)
    rows = [fit_row(fill("cox.alone", gene=gene), cohort, ["z"]),
            fit_row(fill("cox.adjusted", gene=gene, covariates=covariates), cohort,
                    ["z", *COVARIATES])]  # fmt: skip
    writer.table(pl.DataFrame(rows), "cox", columns("columns.cox"))
    shown = tuple(
        (r["model"], str(r["patients"]), str(r["events"]), f"{r['hr']:.2f}",
         fill("cox.ci", low=r["low"], high=r["high"]), f"{r['p']:.2g}", f"{r['ph_p']:.2g}",
         f"{r['c_index']:.3f}")
        for r in rows
    )  # fmt: skip
    heads = text.words(PAGE, "cox.columns")
    table = ResultTable(fill("cox.title", endpoint=ENDPOINT, gene=gene), heads, shown,
                        numeric=(1, 2, 3, 5, 6, 7))  # fmt: skip
    lead = fill("results.cox_lead", gene=gene, hr_uni=rows[0]["hr"], hr_adj=rows[1]["hr"],
                covariates=covariates)  # fmt: skip
    return Group(fill("results.cox"), (table,), lead=lead), rows[1]


def scan_group(
    writer: Writer, cohort: pl.DataFrame, numbers: dict[str, Any], rng: np.random.Generator
) -> Group:
    """Log-rank p at every cut-off, beside p = alpha before and after correcting for the scan.

    The correction is the maximally selected statistic's permutation p (Lausen and Schumacher
    1992): expression is shuffled against the outcomes and the scan repeated each time.
    """
    months, event, z = (cohort[c].to_numpy() for c in ("months", "event", "z"))
    cuts = np.arange(SCAN[0], SCAN[1] + 1)

    def scan(values: np.ndarray) -> np.ndarray:
        """Log-rank chi-squared at every cut-off of these values."""
        return logrank(months, event, np.stack([split(values, float(q)) for q in cuts], axis=1))

    observed = scan(z)
    maxima = np.array([scan(rng.permutation(z)).max() for _ in range(PERMUTATIONS)])
    beaten = int(np.sum(maxima >= observed.max()))
    corrected = (1 + beaten) / (PERMUTATIONS + 1)
    p = chi2.sf(observed, 1)
    values = {**numbers, "low": SCAN[0], "high": SCAN[1], "alpha": ALPHA, "cuts": len(cuts),
              "perms": PERMUTATIONS, "min_p": f"{p.min():.2g}", "best": int(cuts[p.argmin()]),
              "corrected": f"{corrected:.2g}", "floor": f"{1 / (PERMUTATIONS + 1):.2g}",
              "median_p": f"{p[cuts == MEDIAN][0]:.2g}"}  # fmt: skip
    line = fill("figures.scan.line")
    frame = pl.DataFrame({"line": line, "percentile": cuts, "chi2": observed, "p": p,
                          "cutoff": np.percentile(cohort["log2cpm"].to_numpy(), cuts),
                          "neglog10p": -np.log10(p)})  # fmt: skip
    ends = (float(SCAN[0]), float(SCAN[1]))
    needed = float(-np.log10(chi2.sf(np.quantile(maxima, 1 - ALPHA), 1)))
    data = CurveData(
        frame=frame, columns=columns("columns.scan"), series="line", x="percentile",
        y="neglog10p", x_label=fill("figures.scan.x", **values), y_label=fill("figures.scan.y"),
        colour="line", colours={line: ROLES["comparison"]}, hover=("p", "cutoff"),
        caption=fill("figures.scan.caption", **values),
        thresholds=(
            Line(fill("figures.scan.uncorrected", **values), ends, (-np.log10(ALPHA),) * 2, True),
            Line(fill("figures.scan.corrected", **values), ends, (needed, needed), True),
        ),
    )  # fmt: skip
    words = FigureText(title=fill("figures.scan.title", **values),
                       lead=fill("figures.scan.lead" if beaten else "figures.scan.lead_floor",
                                 **values),
                       reading=lines("figures.scan.reading", **values))  # fmt: skip
    return Group(fill("results.scan"), (writer.figure(data, "cut_scan", words),))


def power_group(writer: Writer, numbers: dict[str, Any], adjusted: dict[str, Any]) -> Group:
    """Power over events and hazard ratio per SD, this cohort's events and HR marked."""
    hrs = np.round(np.arange(1.05, 2.01, 0.05), 2)
    grid = pl.DataFrame([(d, float(h)) for d in range(10, 401, 10) for h in hrs],
                        schema=["events", "hr"], orient="row")  # fmt: skip
    grid = grid.with_columns(
        pl.Series("power", power(grid["events"].to_numpy(), grid["hr"].to_numpy()))
    )
    hr = max(adjusted["hr"], 1 / adjusted["hr"])  # power is the same either side of 1
    values = {**numbers, "alpha": ALPHA, "target": TARGET_POWER, "hr": hr,
              "power": float(power(np.array(numbers["events"]), np.array(hr)))}  # fmt: skip
    data = SurfaceData(
        frame=grid, columns=columns("columns.power"), x="events", y="hr", z="power",
        x_label=fill("figures.power.x"), y_label=fill("figures.power.y"),
        z_label=fill("figures.power.z"), caption=fill("figures.power.caption", **values),
        isolines=(TARGET_POWER,), isoline_format=".0f",
        mark=Mark(float(numbers["events"]), hr, fill("figures.power.mark")),
    )  # fmt: skip
    words = FigureText(title=fill("figures.power.title"),
                       lead=fill("figures.power.lead", **values),
                       reading=lines("figures.power.reading", **values))  # fmt: skip
    return Group(fill("results.power"), (writer.figure(data, "power", words),))


def cramers_v(sheet: pl.DataFrame, a: str, b: str) -> tuple[float | None, int]:
    """Cramér's V of two columns over samples with both, uncorrected; None if one is constant."""
    both = sheet.select(a, b).drop_nulls().filter((pl.col(a) != "") & (pl.col(b) != ""))
    counts = both.group_by(a, b).len().pivot(on=b, index=a, values="len").fill_null(0)
    table = counts.drop(a).to_numpy()
    if min(table.shape) < 2:  # noqa: PLR2004 - a constant column has no association
        return None, both.height
    statistic = chi2_contingency(table, correction=False).statistic
    return float(np.sqrt(statistic / (both.height * (min(table.shape) - 1)))), both.height


def confounding_group(writer: Writer, sheet: pl.DataFrame, dataset: str) -> Group:
    """Cramér's V between every pair of the sheet's design columns, binned by Cohen's steps."""
    names = [fill("figures.confounding.below", a=V_EDGES[0]),
             *(fill("figures.confounding.between", a=a, b=b)
               for a, b in zip(V_EDGES, V_EDGES[1:], strict=False)),
             fill("figures.confounding.above", a=V_EDGES[-1])]  # fmt: skip
    rows = []
    for a in CONFOUNDERS:
        for b in CONFOUNDERS:
            v, n = (None, sheet.height) if a == b else cramers_v(sheet, a, b)
            level = None if v is None else names[int(np.searchsorted(V_EDGES, v, "right"))]
            rows.append((a, b, v, level, n))
    frame = pl.DataFrame(rows, schema=list(columns("columns.confounding")), orient="row")
    top = frame.drop_nulls("v").sort("v", descending=True).row(0, named=True)
    values = {"source": fill("preview.sheet", dataset=dataset), "n": sheet.height,
              "a": top["row_variable"], "b": top["column_variable"], "v": top["v"]}  # fmt: skip
    shades = zip(names[1:], HEAT_ALPHAS, strict=True)
    data = HeatmapData(
        frame=frame, columns=columns("columns.confounding"), row="row_variable",
        column="column_variable", value="v", level="level",
        row_label=fill("figures.confounding.row"), column_label=fill("figures.confounding.column"),
        swatches=(Swatch(names[0], None), *(Swatch(n, ROLES["comparison"], a) for n, a in shades)),
        caption=fill("figures.confounding.caption", **values),
        row_order=CONFOUNDERS, column_order=CONFOUNDERS,
    )  # fmt: skip
    words = FigureText(title=fill("figures.confounding.title"),
                       lead=fill("figures.confounding.lead", **values),
                       reading=lines("figures.confounding.reading"))  # fmt: skip
    return Group(fill("results.confounding"), (writer.figure(data, "confounding", words),))


def key_numbers(cohort: pl.DataFrame, adjusted: dict[str, Any]) -> tuple[KeyNumber, ...]:
    """Patients, events, median follow-up by reverse Kaplan-Meier, and the adjusted HR."""
    events = int(cohort["event"].sum())
    reverse = kaplan_meier(cohort["months"].to_numpy(), 1 - cohort["event"].to_numpy())
    follow_up = first_below_half(reverse, "survival")
    return (
        KeyNumber(fill("keys.patients"), f"{cohort.height:,}"),
        KeyNumber(fill("keys.events", endpoint=ENDPOINT), f"{events:,}",
                  fill("keys.events_detail", pct=100 * events / cohort.height)),
        KeyNumber(fill("keys.follow_up"), fill("keys.follow_up_value", months=follow_up),
                  fill("keys.follow_up_detail")),
        KeyNumber(fill("keys.hr"), f"{adjusted['hr']:.2f}",
                  fill("keys.hr_detail", level=100 * LEVEL, low=adjusted["low"],
                       high=adjusted["high"])),
    )  # fmt: skip


def cards() -> tuple[ParameterCard, ...]:
    """The options ``ctx clinical`` will take, defaults filled from this script's constants."""
    defaults = {"endpoint": ENDPOINT, "horizon": HORIZON, "covariates": ",".join(COVARIATES)}
    return tuple(
        ParameterCard(flag, fill(f"planned.{flag}.meaning"),
                      fill(f"planned.{flag}.default", **defaults), fill(f"planned.{flag}.effect"))
        for flag in text.raw(PAGE, "planned")
    )  # fmt: skip


def page(root: Path, dataset: str, day: date) -> Path:
    """Write the M10 preview page of one dataset under ``outputs/12_clinical/``.

    Args:
        root: The ctx root holding the made-up dataset's gather outputs.
        dataset: The dataset's name, as the gather outputs carry it.
        day: The date in every file name.

    Returns:
        The page's path.
    """
    started = time.perf_counter()
    gather = root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    sheet_file = newest(gather, dataset, "samplesheet")
    sheet = pl.read_parquet(sheet_file)
    gene = str(pl.read_parquet(newest(gather, dataset, "genes"))["gene_name"][0])
    out = root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    release = {fill("preview.release_name"): fill("preview.release", seed=SEED)}
    writer = Writer(out, dataset, day, release)
    rng = np.random.default_rng(SEED)
    cohort = made_up_cohort(rng)
    writer.table(cohort, "cohort", columns("columns.cohort"))
    survival, numbers = survival_group(writer, cohort, gene)
    cox_results, adjusted = cox_group(writer, cohort, gene)
    results: list[Result] = [
        text.word("ui", "preview.note"), survival, cox_results,
        scan_group(writer, cohort, numbers, rng), power_group(writer, numbers, adjusted),
        confounding_group(writer, sheet, dataset),
    ]  # fmt: skip
    command = (f"ctx clinical --dataset {dataset} --endpoint {ENDPOINT} --split median "
               f"--horizon {HORIZON} --covariates {','.join(COVARIATES)}")  # fmt: skip
    usage, budget = runtime.measure(started), runtime.budget(MODULE)
    built = ModulePage(
        title=fill("title", dataset=text.word(*SHOWN)), **text.skeleton(PAGE), parameters=cards(),
        results=tuple(results), key_numbers=key_numbers(cohort, adjusted),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children, profile=fill("preview.profile"),
            budget=fill("preview.budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=(sheet_file.name,),
        ),
    )  # fmt: skip
    return write(built, writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html")
