#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side preview builder; the M10 page on made-up data until ctx clinical exists
# ///
"""Course preview of M10 Clinical and survival, drawn on made-up data for the made-up dataset.

The cohort is made up here from a seeded generator: expression of one gene, age, stage and a
progression-free interval with censoring. Every statistic and figure is ``ctx.clinical``'s, the
built module's one definition: Kaplan-Meier, log-rank, Efron Cox with its PH test and C, the
scan with its permutation p, power and Cramér's V, so the preview cannot disagree with the
page it previews. The confounding map is computed on the made-up dataset's own sample sheet.
The page text is ``content/clinical.toml``. Run through ``build_preview.py --only clinical``.
"""

from __future__ import annotations

import time
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

from ctx import runtime
from ctx.clinical import analyses, cutpoint, figures, score, survival
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
from ctx.viz.heatmap import HeatmapData, Swatch
from ctx.viz.style import HEAT_ALPHAS, ROLES

# The made-up dataset's name on the page, one word for every preview (content/ui.toml).
SHOWN = ("ui", "preview.dataset")

PAGE = MODULE = "clinical"
SEED = 10
PATIENTS = 420
ENDPOINT = "PFI"
HORIZON = 60  # months: the RMST horizon
COVARIATES = ("age", "stage")
DAYS_PER_MONTH = 365.25 / 12
MEDIAN = 50  # the scan's percentile of the pre-specified median split
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
    group = np.array(groups())[score.split(z, "median")]
    return pl.DataFrame({
        "patient": [f"P{i:03d}" for i in range(PATIENTS)], "age": age, "stage": stage,
        "log2cpm": log2cpm, "z": z, "group": group, "days": days,
        "months": days / DAYS_PER_MONTH, "event": (event_days <= censor_days).astype(int),
    })  # fmt: skip


def groups() -> tuple[str, str]:
    """The median split's two groups, low first, in the preview's words."""
    return fill("split.low"), fill("split.high")


def outcome(cohort: pl.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Months and event flags, in cohort order."""
    return cohort["months"].to_numpy().astype(float), cohort["event"].to_numpy().astype(int)


def group_summary(table: pl.DataFrame) -> ResultTable:
    """Patients, events, median with its interval, and RMST, per group (``analyses.groups``)."""

    def months(value: float | None) -> str:
        return fill("summary.not_reached") if value is None else f"{value:.1f}"

    rows = tuple((r["group"], str(r["patients"]), str(r["events"]), months(r["median"]),
                  fill("summary.interval", low=months(r["low"]), high=months(r["high"])),
                  f"{r['rmst']:.1f}") for r in table.iter_rows(named=True))  # fmt: skip
    heads = tuple(fill(f"summary.columns.{k}", level=100 * survival.LEVEL)
                  for k in text.raw(PAGE, "summary.columns"))  # fmt: skip
    title = fill("summary.title", endpoint=ENDPOINT, horizon=HORIZON)
    return ResultTable(title, heads, rows, numeric=(1, 2, 3, 5))


def survival_group(writer: Writer, cohort: pl.DataFrame, gene: str) -> tuple[Group, dict[str, Any]]:
    """Kaplan-Meier figure with its risk rows, and the per-group summary; the numbers they print.

    Every number is ``ctx.clinical``'s: the curves, medians and RMST by ``analyses``, the
    log-rank test by ``survival.logrank``, the figure by ``figures.km``, as the built page.
    """
    months, event = outcome(cohort)
    labels, names = cohort["group"].to_numpy(), groups()
    estimates = analyses.km_frame(months, event, labels, names)
    writer.table(estimates, "km_estimates", figures.columns("columns.km_steps"))
    summary = analyses.groups(months, event, labels, names, HORIZON)
    chi, _, p = survival.logrank(months, event, labels)
    medians = dict(summary.select("group", "median").rows())
    rmst = dict(summary.select("group", "rmst").rows())
    numbers = {
        "source": fill("preview.source"), "n": cohort.height, "events": int(event.sum()),
        "endpoint": ENDPOINT, "gene": gene, "chi2": chi, "p": f"{p:.2g}", "horizon": HORIZON,
        "rmst": rmst[names[1]] - rmst[names[0]],
        **{f"median_{g}": fill("summary.not_reached") if medians[n] is None
           else f"{medians[n]:.1f} months" for g, n in zip(("low", "high"), names, strict=True)},
    }  # fmt: skip
    data = figures.km(estimates, names, fill("figures.km.caption", **numbers), HORIZON)
    words = FigureText(title=fill("figures.km.title", **numbers),
                       lead=fill("figures.km.lead", **numbers),
                       reading=lines("figures.km.reading", **numbers))  # fmt: skip
    items = (writer.figure(data, "km", words), group_summary(summary))
    return Group(fill("results.survival", **numbers), items), numbers


def fit_row(name: str, cohort: pl.DataFrame, used: list[str]) -> dict[str, Any]:
    """One Cox fit by ``analyses.cox_model``: the expression term's HR, interval, p, PH p and C."""
    months, event = outcome(cohort)
    model = analyses.cox_model(months, event, cohort.select(used).to_numpy().astype(float))
    return {"model": name, **vars(model)}


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


def scan_group(writer: Writer, cohort: pl.DataFrame, numbers: dict[str, Any]) -> Group:
    """Log-rank p at every cut-off with the scan's corrected p: ``cutpoint.scan``, as built."""
    months, event = outcome(cohort)
    found = cutpoint.scan(months, event, cohort["log2cpm"].to_numpy())
    values = {**numbers, "low": cutpoint.LOW, "high": cutpoint.HIGH, "alpha": analyses.ALPHA,
              "cuts": len(found.cuts), "perms": cutpoint.PERMUTATIONS,
              "min_p": f"{found.p.min():.2g}", "best": found.best,
              "corrected": f"{found.corrected:.2g}",
              "floor": f"{1 / (cutpoint.PERMUTATIONS + 1):.2g}",
              "median_p": f"{found.p[found.cuts == MEDIAN][0]:.2g}"}  # fmt: skip
    data = figures.scan(found, fill("figures.scan.caption", **values))
    words = FigureText(title=fill("figures.scan.title", **values),
                       lead=fill("figures.scan.lead" if found.beaten else "figures.scan.lead_floor",
                                 **values),
                       reading=lines("figures.scan.reading", **values))  # fmt: skip
    return Group(fill("results.scan"), (writer.figure(data, "cut_scan", words),))


def power_group(writer: Writer, numbers: dict[str, Any]) -> Group:
    """Power over events and hazard ratio per SD by ``figures.power``.

    The point marks this cohort's events and the hazard ratio they detect, as the built page
    marks it, never an observed one.
    """
    events = int(numbers["events"])
    values = {**numbers, "alpha": analyses.ALPHA, "target": analyses.TARGET_POWER,
              "hr": analyses.detectable(events)}  # fmt: skip
    data = figures.power(events, fill("figures.power.caption", **values))
    words = FigureText(title=fill("figures.power.title"),
                       lead=fill("figures.power.lead", **values),
                       reading=lines("figures.power.reading", **values))  # fmt: skip
    return Group(fill("results.power"), (writer.figure(data, "power", words),))


def confounding_group(writer: Writer, sheet: pl.DataFrame, dataset: str) -> Group:
    """Cramér's V between every pair of the sheet's design columns (``analyses.cramers_v``)."""
    names = [fill("figures.confounding.below", a=V_EDGES[0]),
             *(fill("figures.confounding.between", a=a, b=b)
               for a, b in zip(V_EDGES, V_EDGES[1:], strict=False)),
             fill("figures.confounding.above", a=V_EDGES[-1])]  # fmt: skip
    rows = []
    for a in CONFOUNDERS:
        for b in CONFOUNDERS:
            v, n = (None, sheet.height) if a == b else analyses.cramers_v(sheet, a, b)
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
    months, event = outcome(cohort)
    return (
        KeyNumber(fill("keys.patients"), f"{cohort.height:,}"),
        KeyNumber(fill("keys.events", endpoint=ENDPOINT), f"{int(event.sum()):,}",
                  fill("keys.events_detail", pct=100 * event.sum() / cohort.height)),
        KeyNumber(fill("keys.follow_up"),
                  fill("keys.follow_up_value", months=survival.follow_up(months, event)),
                  fill("keys.follow_up_detail")),
        KeyNumber(fill("keys.hr"), f"{adjusted['hr']:.2f}",
                  fill("keys.hr_detail", level=100 * survival.LEVEL, low=adjusted["low"],
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
        scan_group(writer, cohort, numbers), power_group(writer, numbers),
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
