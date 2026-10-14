"""The M10 figures' contracts and words.

Kaplan-Meier curves, expression by ER, the scan, the power map and the confounding map, each
worded from ``content/clinical.toml``.
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
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
# The competing causes of death: the CDR's DSS event, and an OS death that is not one.
CANCER, OTHER = 1, 2
# Stage groups of the box: AJCC I to II against III to IV.
STAGE_GROUPS = ((1, 2), (3, 4))


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


def two_levels(patients: pl.DataFrame, variable: str) -> tuple[np.ndarray, tuple[str, str]]:
    """Each patient's level of a two-level variable, None where it has none; the two names.

    ``er``: ER-negative, the reference, against ER-positive, by IHC. ``stage_group``: stage
    I to II, the reference, against III to IV; no stage is no level.

    Returns:
        The level per patient in table order, and (reference, comparison) in page words.
    """
    if variable == "er":
        names = {v: word(f"er.{v}") for v in ("negative", "positive")}
        levels = [names.get(v) if v else None for v in patients["er"].to_list()]
    else:
        names = {k: word(f"box.{k}") for k in ("early", "late")}
        early, late = STAGE_GROUPS
        levels = [None if s is None else names["early" if s in early else "late"]
                  for s in patients["stage"].to_list()]  # fmt: skip
        assert all(s is None or s in early + late for s in patients["stage"].to_list())  # noqa: S101 - harmonise gives stage 1 to 4
    first, second = names.values()
    return np.array(levels, dtype=object), (first, second)


def box(patients: pl.DataFrame, values: Mapping[str, np.ndarray], variable: str) -> BoxData:
    """log2 CPM of the build genes, by a two-level variable's levels (``two_levels``)."""
    levels, groups = two_levels(patients, variable)
    named = np.flatnonzero(levels != None)  # noqa: E711 - elementwise over an object array
    ids = patients["patient"].to_list()
    frame = pl.DataFrame([{"feature": gene, "group": levels[i], "sample": ids[i],
                           "value": float(v[i])}
                          for gene, v in values.items() for i in named])  # fmt: skip
    return box_data(frame, columns=columns(f"columns.box.{variable}"),
                    feature_label=word("figures.expression.feature"),
                    value_label=word("figures.expression.value"), groups=groups)  # fmt: skip


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


def steps(times: np.ndarray, values: np.ndarray, start: float) -> tuple[list[float], list[float]]:
    """A step function's vertices from time 0: each value held until the next time."""
    xs, ys, last = [0.0], [start], start
    for t, v in zip(times.tolist(), values.tolist(), strict=True):
        xs += [t, t]
        ys += [last, v]
        last = v
    return xs, ys


def competing(competed: "Competing", caption: str) -> CurveData:
    """Cumulative incidence of death from the cancer and from other causes, and 1 - KM of DSS.

    1 - KM of DSS counts other deaths as censored, and so overstates the cancer's incidence.
    """
    lines = {k: word(f"figures.competing.lines.{k}") for k in ("cancer", "other", "naive")}
    parts = []
    for key, (times, values) in competed.curves().items():
        xs, ys = steps(times, values, 0.0)
        parts.append(pl.DataFrame({"line": lines[key], "time": xs, "incidence": ys}))
    colours = dict(zip(lines.values(), (ROLES["comparison"], DATA_COLOURS[2], ROLES["reference"]),
                       strict=True))  # fmt: skip
    return CurveData(
        frame=pl.concat(parts), columns=columns("columns.competing"), series="line", x="time",
        y="incidence", x_label=word("figures.axes.time"), y_label=word("figures.competing.y"),
        colour="line", colours=colours, caption=caption,
    )  # fmt: skip


@dataclass(frozen=True)
class Competing:
    """Deaths by cause over the patients with both OS and DSS.

    The CDR's DSS flag names a death from the cancer; an OS death with DSS 0 is one from other
    causes.

    Attributes:
        months: OS time per patient.
        cause: 0 censored, 1 death from the cancer, 2 death from other causes.
        unknown: OS deaths left out for want of a DSS flag, the cause not known.
    """

    months: np.ndarray
    cause: np.ndarray
    unknown: int

    @classmethod
    def of(cls, patients: pl.DataFrame) -> "Competing":
        """From the patient table; the CDR gives DSS the OS time, which is checked."""
        both = patients.filter(pl.col("OS").is_not_null() & pl.col("DSS").is_not_null())
        if (both["OS_days"] != both["DSS_days"]).any():
            raise ValueError("a patient's DSS time differs from its OS time")
        if ((both["DSS"] == 1) & (both["OS"] == 0)).any():
            raise ValueError("a patient died of the cancer (DSS) without dying (OS)")
        unknown = patients.filter((pl.col("OS") == 1) & pl.col("DSS").is_null()).height
        os_, dss = both["OS"].to_numpy(), both["DSS"].to_numpy()
        cause = np.where(dss == 1, CANCER, np.where(os_ == 1, OTHER, 0))
        return cls(both["OS_months"].to_numpy().astype(float), cause, unknown)

    def curves(self) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        """Incidence of each cause, and 1 - KM of the cancer's with other deaths censored."""
        times, (cancer, other) = survival.cumulative_incidence(
            self.months, self.cause, (CANCER, OTHER)
        )
        km = survival.kaplan_meier(self.months, (self.cause == CANCER).astype(int))
        return {"cancer": (times, cancer), "other": (times, other),
                "naive": (km.time, 1 - km.survival)}  # fmt: skip

    def at(self, horizon: float) -> dict[str, float]:
        """Each curve's value at the horizon, in %."""
        out = {}
        for key, (times, values) in self.curves().items():
            held = values[times <= horizon]
            out[key] = 100 * float(held[-1]) if held.size else 0.0
        return out


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
    numeric = [n for n in ("stage", "grade") if n in names]
    table = patients.with_columns(analyses.age_group(pl.col("age")).alias("age_group"),
                                  *(pl.col(n).cast(pl.String) for n in numeric))  # fmt: skip
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
