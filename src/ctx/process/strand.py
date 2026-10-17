"""M3's forward share per sample, its strand calls redrawn at the reader's thresholds (WO8).

Each sample's call is the registry's ``Metric.strandedness`` on its two featureCounts totals,
in exact arithmetic; the page's twin ``ctxStrand.call`` (``templates/_strand.js.j2``) decides
it the same way, on the same two integers, and tests/test_strand_page.py holds the two
together, the edges where floats fail among them. The figure is a dot per sample, coloured
by its call, with the thresholds as lines; a reader moves the two thresholds and every call,
colour, line and count follows. The figure files draw nf-core/rnaseq's thresholds.
"""

from collections.abc import Mapping
from dataclasses import replace
from typing import Any

import polars as pl

from ctx.qc.registry import Metric
from ctx.report import text
from ctx.report.figures import columnar, finite
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import BAND_COLOURS, ROLES

PAGE = "process"
CALLS = ("forward", "reverse", "unstranded", "undetermined")
# By role: the two stranded calls the comparison and reference colours, an unstranded library
# the good band's green, a sample between the bands the usable band's orange.
CALL_COLOURS = {"forward": ROLES["comparison"], "reverse": ROLES["reference"],
                "unstranded": BAND_COLOURS["good"],
                "undetermined": BAND_COLOURS["usable"]}  # fmt: skip
# The thresholds a reader may type: a stranded share past a coin's half, an unstranded width
# within the whole range.
STRANDED_RANGE, UNSTRANDED_RANGE = (0.5, 1.0), (0.0, 1.0)


def _lines(n: int, stranded: float, unstranded: float) -> tuple[Line, ...]:
    """The thresholds across every sample: forward, reverse and the unstranded band's edges."""
    span = (0.5, n + 0.5)
    half = unstranded / 2
    named = (("forward", stranded), ("reverse", 1 - stranded),
             ("unstranded", 0.5 - half), ("unstranded", 0.5 + half))  # fmt: skip
    return tuple(Line(text.word(PAGE, f"figures.share.line_{k}"), span, (y, y), dash=True)
                 for k, y in named)  # fmt: skip


def share_plot(strand: pl.DataFrame, metric: Metric, caption: str) -> ScatterData:
    """The forward share per sample, coloured by its call, with the live thresholds' settings.

    Args:
        strand: ``tables.strandedness``: sample, assigned_s1, assigned_s2, forward_share, call.
        metric: The registry's ``gene_strand_bias``: its stranded and unstranded thresholds.
        caption: The filled caption.
    """
    assert metric.stranded is not None and metric.unstranded is not None  # noqa: S101 - load() checked
    shown = strand.filter(pl.col("forward_share").is_not_null()).select(
        "sample", s1="assigned_s1", s2="assigned_s2", share="forward_share", call="call",
    ).with_row_index("rank", offset=1).with_columns(pl.col("rank").cast(pl.Int64))  # fmt: skip
    names = {c: text.word(PAGE, f"figures.share.calls.{c}") for c in CALLS}
    columns = {c: (u, text.word(PAGE, f"figures.share.columns.{c}"))
               for c, u in (("rank", "none"), ("sample", "none"), ("s1", "reads"),
                            ("s2", "reads"), ("share", "fraction"), ("call", "none"))}  # fmt: skip
    data = ScatterData(
        frame=shown, columns=columns, x="rank", y="share",
        x_label=text.word(PAGE, "figures.share.x"), y_label=text.word(PAGE, "figures.share.y"),
        colour="call", colours={c: CALL_COLOURS[c] for c in CALLS}, names=names,
        caption=caption, hover=("sample", "s1", "s2"),
        lines=_lines(shown.height, metric.stranded, metric.unstranded),
        labels={c: text.word(PAGE, f"figures.share.labels.{c}") for c in ("sample", "s1", "s2")},
    )  # fmt: skip
    return replace(data, live=live(data, metric))


def live(data: ScatterData, metric: Metric) -> dict[str, Any]:
    """The page script's settings: every sample's two totals, the thresholds and the words."""
    assert metric.stranded is not None and metric.unstranded is not None  # noqa: S101 - load() checked
    rows = data.frame.select("rank", "sample", "s1", "s2", "share").to_dicts()
    return {
        "kind": "strand", "rows": columnar([finite(r) for r in rows]),
        "stranded": f"{metric.stranded:g}", "unstranded": f"{metric.unstranded:g}",
        "ranges": {"stranded": STRANDED_RANGE, "unstranded": UNSTRANDED_RANGE},
        "calls": list(CALLS), "traces": {c: i for i, c in enumerate(data.colours)},
        "lines": [len(data.colours) + i for i in range(4)], "hover": list(data.hover),
        "summary": text.raw(PAGE, "figures.share.summary"),
        "words": {k: text.word(PAGE, f"figures.share.form.{k}")
                  for k in ("label", "stranded", "unstranded", "reset")},
        "moved": text.raw(PAGE, "figures.share.moved"),
        "default_summary": summary(data.frame, metric.stranded, metric.unstranded),
    }  # fmt: skip


def counted(calls: list[str]) -> Mapping[str, int]:
    """Samples per call, every call named."""
    return {c: calls.count(c) for c in CALLS}


def summary(frame: pl.DataFrame, stranded: float, unstranded: float) -> str:
    """The live line at the build's thresholds, as the page fills it again."""
    found = counted(frame["call"].to_list())
    return text.fill(PAGE, "figures.share.summary", n=frame.height, stranded_at=f"{stranded:g}",
                     unstranded_at=f"{unstranded:g}", **found)  # fmt: skip
