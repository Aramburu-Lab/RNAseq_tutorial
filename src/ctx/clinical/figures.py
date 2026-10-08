"""The M10 figures' contracts and words.

Kaplan-Meier curves, expression by ER, the scan, the power map and the confounding map, each
worded from ``content/clinical.toml``.
"""

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import chi2

from ctx.clinical import analyses, survival
from ctx.clinical.cutpoint import Scan
from ctx.report import text
from ctx.report.builder import FigureText
from ctx.viz.box import BoxData, box_data
from ctx.viz.contour import Mark, SurfaceData
from ctx.viz.curve import CurveData
from ctx.viz.heatmap import HeatmapData, Swatch
from ctx.viz.km import KmData
from ctx.viz.scatter import Line
from ctx.viz.style import DATA_COLOURS, HEAT_ALPHAS, KM, ROLES

PAGE = "clinical"
MONTHS_PER_YEAR = 12
V_EDGES = (0.1, 0.3, 0.5)  # Cohen's small, medium and large for one degree of freedom
HRS = tuple(round(1.05 + 0.05 * i, 2) for i in range(20))
EVENTS = tuple(range(10, 401, 10))


def columns(key: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) per column, from the page file."""
    return {name: (unit, meaning) for name, (unit, meaning) in text.raw(PAGE, key).items()}


def word(key: str) -> str:
    """A fixed text of the page file."""
    return text.word(PAGE, key)


def words(name: str, **values: Any) -> FigureText:
    """Title, lead and reading guide of one figure from ``figures.<name>``, filled."""
    reading = tuple(template.format_map(values) for template in
                    text.words(PAGE, f"figures.{name}.reading"))  # fmt: skip
    return FigureText(title=text.fill(PAGE, f"figures.{name}.title", **values),
                      lead=text.fill(PAGE, f"figures.{name}.lead", **values),
                      reading=reading)  # fmt: skip


def colours(groups: Sequence[str]) -> dict[str, str]:
    """Two groups take the reference and comparison roles; three put a data colour between."""
    if len(groups) == 2:  # noqa: PLR2004 - a pair takes the two roles
        return {groups[0]: ROLES["reference"], groups[1]: ROLES["comparison"]}
    return dict(zip(groups, (ROLES["reference"], DATA_COLOURS[2], ROLES["comparison"]),
                    strict=False))  # fmt: skip


def risk_step(last: float) -> float:
    """Months between at-risk counts: whole years, at most ``KM["risk_columns"]`` columns.

    ``ctxLive.riskStep`` is its page twin.
    """
    years = max(1, math.ceil(last / (KM["risk_columns"] - 1) / MONTHS_PER_YEAR))
    return float(MONTHS_PER_YEAR * years)


def km(frame: pl.DataFrame, groups: Sequence[str], caption: str, months: float,
       live: Mapping[str, Any] | None = None) -> KmData:  # fmt: skip
    """A Kaplan-Meier figure in months, at-risk counts every ``risk_step`` months."""
    last = float(frame["time"].max()) if frame.height else months  # type: ignore[arg-type]
    step = risk_step(last)
    times = tuple(float(t) for t in np.arange(0, last + step, step) if t <= last)
    return KmData(
        frame=frame, columns=columns("columns.km_steps"), groups=tuple(groups),
        colours=colours(groups), x_label=word("figures.axes.time"),
        y_label=word("figures.axes.free"), risk_label=word("figures.axes.risk"),
        risk_times=times, caption=caption, live=live,
    )  # fmt: skip


def er_box(patients: pl.DataFrame, values: Mapping[str, np.ndarray]) -> BoxData:
    """log2 CPM of the build genes, ER-negative against ER-positive (IHC)."""
    called = patients["er"].to_numpy()
    named = np.flatnonzero(np.isin(called, ("positive", "negative")))
    names = {v: word(f"er.{v}") for v in ("positive", "negative")}
    ids = patients["patient"].to_list()
    frame = pl.DataFrame([{"feature": gene, "group": names[called[i]], "sample": ids[i],
                           "value": float(v[i])}
                          for gene, v in values.items() for i in named])  # fmt: skip
    return box_data(frame, columns=columns("columns.box"),
                    feature_label=word("figures.expression.feature"),
                    value_label=word("figures.expression.value"),
                    groups=(names["negative"], names["positive"]))  # fmt: skip


def scan(found: Scan, caption: str, live: Mapping[str, Any] | None = None) -> CurveData:
    """-log10 p of every cut, with p = alpha uncorrected and the scan's corrected line."""
    line = word("figures.scan_score.line")
    alpha = analyses.ALPHA
    ends = (float(found.cuts[0]), float(found.cuts[-1]))
    needed = float(-np.log10(chi2.sf(np.quantile(found.maxima, 1 - alpha), 1)))
    frame = pl.DataFrame({"line": line, "percentile": found.cuts, "chi2": found.statistic,
                          "p": found.p, "cutoff": found.thresholds,
                          "neglog10p": -np.log10(found.p)})  # fmt: skip
    return CurveData(
        frame=frame, columns=columns("columns.scan_score"), series="line", x="percentile",
        y="neglog10p", x_label=word("figures.scan_score.x"),
        y_label=word("figures.scan_score.y"), colour="line",
        colours={line: ROLES["comparison"]}, hover=("p", "cutoff"), caption=caption,
        thresholds=(
            Line(text.fill(PAGE, "figures.scan_score.uncorrected", alpha=alpha), ends,
                 (float(-np.log10(alpha)),) * 2, True),
            Line(text.fill(PAGE, "figures.scan_score.corrected", alpha=alpha), ends,
                 (needed, needed), True),
        ),
        live=live,
    )  # fmt: skip


def power(events: int, caption: str) -> SurfaceData:
    """Power over events and per-SD hazard ratio.

    The mark sits at this endpoint's events and the hazard ratio they detect with 80 % power,
    never at an observed one (GS plan: power at the observed HR is a function of its p).
    """
    grid = pl.DataFrame([(d, h) for d in EVENTS for h in HRS], schema=["events", "hr"],
                        orient="row")  # fmt: skip
    grid = grid.with_columns(pl.Series("power", survival.power(
        grid["events"].to_numpy(), grid["hr"].to_numpy(), analyses.ALPHA)))  # fmt: skip
    mark = Mark(float(events), analyses.detectable(events), word("figures.power_map.mark"))
    return SurfaceData(
        frame=grid, columns=columns("columns.power_map"), x="events", y="hr", z="power",
        x_label=word("figures.power_map.x"), y_label=word("figures.power_map.y"),
        z_label=word("figures.power_map.z"), caption=caption,
        isolines=(analyses.TARGET_POWER,), isoline_format=".0f", mark=mark,
    )  # fmt: skip


def confounding(patients: pl.DataFrame, names: Sequence[str]
                ) -> tuple[HeatmapData, dict[str, Any]]:  # fmt: skip
    """Cramér's V between every pair of patient variables, binned by Cohen's steps.

    The caption is ``figures.confounding_map.caption``, filled with the patients and the
    strongest pair, which the call also returns.
    """
    table = patients.with_columns(analyses.age_group(pl.col("age")).alias("age_group"),
                                  pl.col("stage").cast(pl.String))  # fmt: skip
    bins = [text.fill(PAGE, "figures.confounding_map.below", a=V_EDGES[0]),
            *(text.fill(PAGE, "figures.confounding_map.between", a=a, b=b)
              for a, b in zip(V_EDGES, V_EDGES[1:], strict=False)),
            text.fill(PAGE, "figures.confounding_map.above", a=V_EDGES[-1])]  # fmt: skip
    label = {n: word(f"patients.variables.{n}") for n in names}
    rows = []
    for a in names:
        for b in names:
            v, n = (None, table.height) if a == b else analyses.cramers_v(table, a, b)
            level = None if v is None else bins[int(np.searchsorted(V_EDGES, v, "right"))]
            rows.append((label[a], label[b], v, level, n))
    frame = pl.DataFrame(rows, schema=list(columns("columns.confounding_map")), orient="row")
    top = frame.drop_nulls("v").sort("v", descending=True).row(0, named=True)
    caption = text.fill(PAGE, "figures.confounding_map.caption", n=table.height)
    shades = zip(bins[1:], HEAT_ALPHAS, strict=True)
    data = HeatmapData(
        frame=frame, columns=columns("columns.confounding_map"), row="row_variable",
        column="column_variable", value="v", level="level",
        row_label=word("figures.confounding_map.row"),
        column_label=word("figures.confounding_map.column"),
        swatches=(Swatch(bins[0], None), *(Swatch(b, ROLES["comparison"], a) for b, a in shades)),
        caption=caption, row_order=tuple(label[n] for n in names),
        column_order=tuple(label[n] for n in names),
    )  # fmt: skip
    return data, top
