"""M2's band map flag view (V2): every sample on every metric the run judged, as a view.

The figure controls plan's question 7 recommended the flag map as the band map's default; the
build keeps the absolute bands and offers the flags as a view, looked up at the run's k and
strata (D-073). A cell is flagged where ``checkpoint.judge`` flagged that sample on that metric,
by its stratum's median ± k MAD or by its absolute band; not flagged where a band judged it;
not judged where its stratum held too few samples for a band and the metric has no absolute
band; no value where the sample has none. Metrics the registry keeps are judged and never
flagged, so their rows are left out and the caption counts them. Every word is
``content/qc.toml``'s ``figures.flag_map``.
"""

import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import polars as pl

from ctx.qc.checkpoint import STRATUM, _stratum_label, strata_of
from ctx.qc.registry import Registry
from ctx.report import text
from ctx.report.builder import MAX_BYTES, FigureBlock, Group, ModulePage, rendered
from ctx.viz.heatmap import HeatmapData, Swatch
from ctx.viz.style import BAND_COLOURS

PAGE = "qc"
LEVELS = ("flagged", "passed", "unjudged", "none")
# "No value" is the band map's grey, "not judged" the same grey lighter: neither is a band.
UNJUDGED_ALPHA = 0.45


def level(value: object, flagged: bool, judged: bool) -> str:
    """A cell's level: no value, flagged, passed where a band judged it, else not judged.

    The one rule the map's cells follow; a flag needs a value, as ``checkpoint.judge`` gives
    one only on a value.
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "none"
    if flagged:
        return "flagged"
    return "passed" if judged else "unjudged"


def flag_map(values: pl.DataFrame, judged: tuple[pl.DataFrame, pl.DataFrame],
             units: pl.DataFrame, registry: Registry,
             strata: tuple[str, ...]) -> HeatmapData | None:  # fmt: skip
    """Every sample's judgement on every metric the run flags by, one row per metric.

    Args:
        values: ``sample`` then one column per metric, as ``checkpoint.judge`` read them.
        judged: The run's bands and flags, from ``checkpoint.judge``.
        units: One row per sample with every strata column, ``LENGTH_GROUP`` among them.
        registry: The loaded registry.
        strata: The run's ``--strata``.

    Returns:
        The map; None where no metric can flag.
    """
    bands, flags = judged
    metrics = [m for m in registry.metrics.values()
               if m.column in values.columns and m.action != "keep"]  # fmt: skip
    if not metrics:
        return None
    kept = sum(m.column in values.columns and m.action == "keep"
               for m in registry.metrics.values())  # fmt: skip
    flagged = set(flags.select("sample", "metric").iter_rows())
    banded = {(m, s) for m, s, low, high in bands.select("metric", STRATUM, "low", "high").rows()
              if low is not None or high is not None}  # fmt: skip
    words = {k: text.word(PAGE, f"figures.flag_map.levels.{k}") for k in LEVELS}
    rows = []
    for row in values.join(units, on="sample", how="left").iter_rows(named=True):
        for m in metrics:
            stratum = _stratum_label(row, strata_of(m, strata))
            judged_here = m.kind != "numeric" or m.banded or (m.id, stratum) in banded
            cell = level(row[m.column], (row["sample"], m.id) in flagged, judged_here)
            rows.append((row["sample"], m.title, words[cell]))
    frame = pl.DataFrame(rows, schema=["sample", "metric", "level"], orient="row").with_columns(
        value=pl.lit(None, pl.Float64))  # fmt: skip
    order = tuple(sorted(str(u) for u in values["sample"].unique()))
    colours = {"flagged": BAND_COLOURS["bad"], "passed": BAND_COLOURS["good"]}
    swatches = tuple(Swatch(words[k], colours.get(k), UNJUDGED_ALPHA if k == "unjudged" else 1.0)
                     for k in LEVELS)  # fmt: skip
    counts = {k: int((frame["level"] == words[k]).sum()) for k in LEVELS}
    caption = text.fill(PAGE, "figures.flag_map.caption", n=len(order), m=len(metrics),
                        k=f"{registry.mad_multiplier:g}", kept=kept, **counts)  # fmt: skip
    return HeatmapData(
        frame=frame, row="metric", column="sample", value="value", level="level",
        columns={"sample": ("none", text.word(PAGE, "columns.sample")),
                 "metric": ("none", text.word(PAGE, "figures.flag_map.metric")),
                 "level": ("none", text.word(PAGE, "figures.flag_map.level")),
                 "value": ("none", text.word(PAGE, "figures.band_map.value"))},
        row_label=text.word(PAGE, "figures.band_map.row_label"),
        column_label=text.word(PAGE, "figures.band_map.hover.sample"),
        column_title=text.fill(PAGE, "figures.band_map.x",
                               unit=text.word(PAGE, "figures.band_map.unit.sample"), n=len(order)),
        column_order=order, swatches=swatches, caption=caption, width="double",
    )  # fmt: skip


def fitted(page: ModulePage, path: Path) -> ModulePage:
    """The page, or, where it would pass ``builder.MAX_BYTES``, its band map without the view.

    A note in the view's place says so. On Step-1 the view added 709,561 B to a 14,193,047 B
    page (2026-10-19), 97,392 B under the cap: a larger cohort can cross it, and a page over
    the cap stops the run.
    """
    if len(rendered(page, path).encode()) <= MAX_BYTES:
        return page
    return replace(page, results=tuple(_without_view(r) for r in page.results))


def _without_view(result: Any) -> Any:
    """A group with its band map's view dropped and the note after it; anything else as is."""
    if not isinstance(result, Group):
        return result
    items: list[Any] = []
    for item in result.items:
        if isinstance(item, FigureBlock) and item.stem.endswith("_band_map") and (
                item.controls or {}).get("kind") == "view":  # fmt: skip
            items += [
                replace(item, controls=None),
                text.fill(PAGE, "figures.flag_map.left_out", cap=MAX_BYTES),
            ]
        else:
            items.append(item)
    return replace(result, items=tuple(items))


def counts(data: HeatmapData) -> dict[str, int]:
    """The map's cells per level, by level key, for its lead."""
    words = {text.word(PAGE, f"figures.flag_map.levels.{k}"): k for k in LEVELS}
    found = dict(data.frame["level"].value_counts().rows())
    return {words[w]: int(found.get(w, 0)) for w in words} | {"cells": data.frame.height}
