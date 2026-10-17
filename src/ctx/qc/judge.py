"""M2's judge panel (work-order item 5): what the page needs to judge a strip again.

A strip's metric is judged by ``checkpoint.judge`` at the run's k and strata; the page judges
it again at another k, within other strata, and for agreement and distance with other peers,
through ``ctxJudge`` (templates ``_qc.js.j2``), the twin ``tests/test_qc_judge.py`` holds to
``checkpoint.judge`` under node (D-012). What does not depend on k or the strata is decided
here once: a value's absolute band and a strand call, whose exact arithmetic stays in Python.
The status table, the key numbers and the confounding tests stay at the run's own judgement.
"""

import math
from collections.abc import Sequence
from typing import Any

import polars as pl

from ctx.qc import checkpoint
from ctx.qc.registry import Metric, Registry
from ctx.report import text
from ctx.viz.strip import SPREAD, StripData
from ctx.viz.style import BAND_ALPHA, LINE_WIDTHS, MARKERS, STRIP, THEMES, tint

PAGE = "qc"
# The multipliers the panel offers around the registry's 3 (figure controls plan).
K_STEPS = (2.0, 2.5, 3.0, 3.5, 4.0, 5.0)
# The values the page takes per peer grouping, which ``Measured.by_peers`` holds.
PEERED = ("replicate_spearman", "pca_distance")


def grouping_label(columns: Sequence[str]) -> str:
    """A grouping as the selects name it: its columns joined, or the word for none."""
    return " and ".join(columns) if columns else text.word(PAGE, "judge.no_strata")


def unit_table(values: pl.DataFrame, units: pl.DataFrame, key: str,
               columns: Sequence[str]) -> pl.DataFrame:  # fmt: skip
    """Every unit of the values, in their order, with each grouping column as the judge names it.

    Text, a missing value read "None", as ``checkpoint._stratum_label`` prints it.
    """
    kept = [c for c in dict.fromkeys(columns) if c in units.columns]
    return values.select(key).join(units.select(key, *kept), on=key, how="left",
                                   maintain_order="left").with_columns(
        [pl.col(c).cast(pl.String).fill_null("None") for c in kept])  # fmt: skip


def decided(frame: pl.DataFrame, metric: Metric) -> tuple[dict[int, str], list[str | None]]:
    """What k and the strata never change: each unit's bad-band reason and strand call.

    The reasons come from ``checkpoint._absolute`` for a banded metric; a strand call from
    ``Metric.strandedness``, in its exact arithmetic. ``frame`` holds ``sample`` and the
    metric's column, one row per unit in the page's order.
    """
    index = {s: i for i, s in enumerate(frame["sample"])}
    absolute: dict[int, str] = {}
    if metric.banded and metric.kind == "numeric":
        absolute = {index[s]: r for s, r in checkpoint._absolute(frame, metric, "")}  # noqa: SLF001
    calls: list[str | None] = [None] * frame.height
    if metric.kind == "strandedness":
        calls = [None if v is None else metric.strandedness(float(v))
                 for v in frame[metric.column]]  # fmt: skip
    return absolute, calls


def _metric(metric: Metric) -> dict[str, Any]:
    """The metric's part of the judge, as ``ctxJudge.judge`` reads it."""
    return {"id": metric.id, "title": metric.title, "kind": metric.kind, "scale": metric.scale,
            "direction": metric.direction, "keep": metric.action == "keep",
            "banded": metric.banded,
            "within": metric.judged_within == checkpoint.LENGTH_GROUP}  # fmt: skip


def _words(key: str) -> dict[str, str]:
    """The reason and note templates the twin fills, and the panel's own lines."""
    raw = {f: text.raw(PAGE, f"reasons.{f}") for f in ("low", "high", "strand")}
    raw |= {"small": text.raw(PAGE, "notes.small_stratum"),
            "spread": text.word(PAGE, "notes.no_spread"),
            "none": text.word(PAGE, "strata.none"),
            "lead": text.raw(PAGE, "figures.strip.lead"),
            "line": text.raw(PAGE, "judge.line"), "missing": text.raw(PAGE, "judge.missing"),
            "unit": text.word(PAGE, f"figures.band_map.unit.{key}"),
            **{f: text.word(PAGE, f"judge.{f}") for f in ("label", "k", "strata", "peers",
                                                           "noscript")}}  # fmt: skip
    return raw


def _draw(data: StripData) -> dict[str, Any]:
    """Where the strip's traces sit and how a stratum's bar is drawn (``StripData``)."""
    light = THEMES["white"]
    fixed = len(data.spans) + len(data.references)
    return {"fixed": fixed, "points": fixed + 2 * len(data.ranges), "jitter": STRIP["jitter"],
            "half": STRIP["range"], "spread": SPREAD,
            "symbols": [MARKERS[0][1], MARKERS[1][1]],
            "range": {"name": data.range_label, "fill": tint(light.reference, BAND_ALPHA),
                      "colour": light.reference, "width": LINE_WIDTHS["line"]}}  # fmt: skip


def settings(  # noqa: PLR0913 - one strip's panel: its metric, values, units and the run's rule
    metric: Metric,
    data: StripData,
    *,
    values: pl.DataFrame,
    units: pl.DataFrame,
    key: str,
    registry: Registry,
    options: tuple[list[tuple[str, ...]], list[tuple[str, ...]]],
    by_peers: Sequence[pl.DataFrame] = (),
) -> dict[str, Any]:
    """The page's judge panel for one strip, ``kind`` "judge".

    Args:
        metric: The strip's metric.
        values: ``key`` and the metric's column, every unit the run judged, in order.
        units: ``key`` and every grouping column, one row per unit.
        key: ``sample`` or ``run_accession``.
        registry: The registry at the run's k.
        data: The strip as built, whose traces the page redraws.
        options: The strata and the peer groupings (``checkpoint.groupings``), the run's first.
        by_peers: Each peer grouping's agreement and distance, for a peered metric.

    Returns:
        The settings ``ctxJudge`` and the page's panel read.
    """
    strata, peers = options
    columns = [c for o in (*strata, *peers) for c in o] + [checkpoint.LENGTH_GROUP]
    table = unit_table(values, units, key, columns)
    frame = values.select(pl.col(key).alias("sample"), metric.column)
    absolute, calls = decided(frame, metric)
    peered = metric.column in PEERED and len(by_peers) > 1
    cfg: dict[str, Any] = {
        "kind": "judge", "scope": key, "metric": _metric(metric),
        "units": {c: table[c].to_list() for c in table.columns},
        "values": [None if v is None or math.isnan(v) else v
                   for v in frame[metric.column].cast(pl.Float64)],
        "absolute": {str(i): r for i, r in absolute.items()}, "calls": calls,
        "ks": list(K_STEPS), "k": registry.mad_multiplier,
        "strata": [[grouping_label(o), list(o)] for o in strata],
        "min_samples": registry.min_samples, "scaled": registry.mad_scaled,
        "words": _words(key), "draw": _draw(data),
    }  # fmt: skip
    if peered:
        cfg["peers"] = [[grouping_label(o), list(o)] for o in peers]
        cfg["peered"] = [values.select(key).join(
            p.select(pl.col("sample").alias(key), metric.column), on=key, how="left",
            maintain_order="left")[metric.column].to_list() for p in by_peers]  # fmt: skip
    return cfg
