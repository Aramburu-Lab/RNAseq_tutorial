"""FastQC's own plots behind its status verdicts, read through MultiQC, one line per unit.

The status tabs of M2 judge FastQC's pass, warn and fail; these figures show what each verdict
was made on (D-052): mean quality, adapter and N content along the read, the per-read GC
histogram over MultiQC's theoretical curve, and the top overrepresented sequence's share,
each over the display bands the registry cites. MultiQC 1.35 writes no table of per-base
sequence content, so that verdict keeps no figure; the page says so. Every number comes from
MultiQC's plot data (D-028), never from a FastQC zip.
"""

import re
from collections.abc import Collection, Mapping
from dataclasses import dataclass, replace

import polars as pl

from ctx.qc import figures, multiqc_io
from ctx.qc.registry import Metric
from ctx.report import text
from ctx.viz.curve import CurveData
from ctx.viz.scatter import Line
from ctx.viz.strip import Span, StripData
from ctx.viz.style import MAX_DISTINCT, category_colours

PAGE = "qc"
# The registry metric each MultiQC plot shows, in the order the page draws them.
PLOTS = {
    "fastqc_per_base_quality": "fastqc_per_base_sequence_quality_plot",
    "fastqc_adapter": "fastqc_adapter_content_plot",
    "fastqc_n_content": "fastqc_per_base_n_content_plot",
    "fastqc_gc_distribution": "fastqc_per_sequence_gc_content_plot_Percentages",
}
OVERREPRESENTED = ("fastqc_overrepresented", "fastqc_overrepresented_sequences_plot",
                   "Top overrepresented sequence")  # fmt: skip
# The status whose plot MultiQC 1.35 does not export.
NO_TABLE = "fastqc_sequence_content"
# MultiQC's own row for the theoretical GC curve its config names (registry.multiqc_config).
THEORETICAL = "Theoretical GC Content"
# Adapter rows are "<sample> - <adapter>"; FastQC judges the largest adapter at each position.
ADAPTER_SEPARATOR = " - "
# A cell of MultiQC's FastQC plot data: "(x, y)".
POINT = re.compile(r"\(\s*([-0-9.eE]+)\s*,\s*([-0-9.eE]+)\s*\)")


@dataclass(frozen=True)
class Plots:
    """FastQC's plot data for one set of units, runs or samples.

    Attributes:
        key: The unit column: ``run_accession`` or ``sample``.
        curves: Per registry metric, ``key``, ``x`` and ``y``, one row per point of read 1.
        theoretical: MultiQC's theoretical GC curve, ``x`` and ``y``; None without one.
        overrepresented: ``key`` and ``top_pct``, the top sequence's share of read 1.
    """

    key: str
    curves: Mapping[str, pl.DataFrame]
    theoretical: pl.DataFrame | None
    overrepresented: pl.DataFrame


def _points(row: tuple[str | None, ...]) -> list[tuple[float, float]]:
    """The (x, y) pairs of one row of MultiQC's plot data."""
    return [(float(x), float(y)) for cell in row if cell for x, y in POINT.findall(cell)]


def read(run: multiqc_io.MultiqcRun, names: Mapping[str, str], first: Collection[str],
         key: str) -> Plots:  # fmt: skip
    """Every FastQC plot MultiQC wrote, for the read-1 rows of the units ``names`` maps.

    Args:
        run: The MultiQC run over the FastQC reports.
        names: Every name MultiQC gives a FastQC report, both mates, to its unit.
        first: The names of read 1, or of a single-end run's one file.
        key: The unit column of the tables returned.

    Returns:
        The curves, the theoretical GC curve and the overrepresented shares.

    Raises:
        QcError: If a MultiQC row names no known report.
    """
    curves, theoretical = {}, None
    for metric, plot in PLOTS.items():
        path = run.data / f"{plot}.txt"
        if not path.is_file():
            continue
        table = pl.read_csv(path, separator="\t", infer_schema=False)
        labels = [str(s).split(ADAPTER_SEPARATOR)[0] for s in table["Sample"]]
        multiqc_io.map_names([n for n in labels if n != THEORETICAL], names)
        rows = []
        for label, row in zip(labels, table.rows(), strict=True):
            if label == THEORETICAL:
                theoretical = pl.DataFrame(_points(row[1:]), schema=["x", "y"], orient="row")
            elif label in first:
                rows += [(names[label], x, y) for x, y in _points(row[1:])]
        frame = pl.DataFrame(rows, schema={key: pl.String, "x": pl.Float64, "y": pl.Float64},
                             orient="row")  # fmt: skip
        curves[metric] = frame.group_by(key, "x", maintain_order=True).agg(pl.col("y").max())
    return Plots(key, curves, theoretical, _overrepresented(run, names, first, key))


def _overrepresented(run: multiqc_io.MultiqcRun, names: Mapping[str, str],
                     first: Collection[str], key: str) -> pl.DataFrame:  # fmt: skip
    """The top overrepresented sequence's share of each unit's read 1, in %."""
    path = run.data / f"{OVERREPRESENTED[1]}.txt"
    empty = pl.DataFrame(schema={key: pl.String, "top_pct": pl.Float64})
    if not path.is_file():
        return empty
    table = pl.read_csv(path, separator="\t", infer_schema=False)
    multiqc_io.map_names(table["Sample"], names)
    table = multiqc_io.numeric(table.filter(pl.col("Sample").is_in(list(first))),
                               [OVERREPRESENTED[2]])  # fmt: skip
    return table.select(pl.col("Sample").replace_strict(names).alias(key),
                        pl.col(OVERREPRESENTED[2]).alias("top_pct"))  # fmt: skip


def _spans(metric: Metric) -> tuple[Span, ...]:
    """The display bands FastQC's limits give, labelled as FastQC's verdicts read."""
    word = {b: text.word(PAGE, f"figures.fastqc.bands.{b}") for b in ("good", "usable", "bad")}
    return tuple(Span(b, word[b], low, high) for b, low, high in metric.plot_spans())


def curve_figure(plots: Plots, metric: Metric, units: pl.DataFrame,
                 colour: str | None) -> CurveData | None:  # fmt: skip
    """One FastQC plot, a line per unit over the metric's display bands; None without data.

    Args:
        plots: ``read``'s plot data.
        metric: The status metric the plot shows.
        units: ``plots.key``, ``sample`` and the colour column, one row per unit.
        colour: The sheet column colouring the lines, or None for one colour.

    Returns:
        The curves; each unit its own colour and dash while they fit, else its group's.
    """
    frame = plots.curves.get(metric.id)
    if frame is None or frame.is_empty():
        return None
    key, words = plots.key, f"figures.fastqc.{metric.id}"
    extra = [c for c in ("sample", colour) if c and c != key and c in units.columns]
    frame = frame.join(
        units.select(key, *extra).unique(key, keep="first", maintain_order=True),
        on=key,
        how="left",
        maintain_order="left",
    )
    if colour is None:
        colour = "group"
        frame = frame.with_columns(group=pl.lit(text.word(PAGE, "figures.fates.all")))
    frame = frame.with_columns(pl.col(colour).cast(pl.String))
    distinct = frame[key].n_unique() <= MAX_DISTINCT
    units_of = text.raw(PAGE, f"{words}.units")
    columns = {c: (units_of.get(c, "none"), text.fill(PAGE, "figures.fates.colour", column=c)
                   if c == colour else text.word(PAGE, f"figures.fastqc.columns.{c}"))
               for c in frame.columns}  # fmt: skip
    theoretical = plots.theoretical if metric.id == "fastqc_gc_distribution" else None
    lines = () if theoretical is None else (Line(
        text.word(PAGE, "figures.fastqc.theoretical"), tuple(theoretical["x"].to_list()),
        tuple(theoretical["y"].to_list())),)  # fmt: skip
    levels = [str(v) for v in frame[colour].unique(maintain_order=True)]
    hover = tuple(c for c in ("sample", colour) if c in frame.columns and c not in (key, "group"))
    return CurveData(
        frame=frame, columns=columns, series=key, x="x", y="y",
        x_label=text.word(PAGE, f"{words}.x"), y_label=text.word(PAGE, f"{words}.y"),
        colour=key if distinct else colour, colours=category_colours(levels),
        distinct=distinct, hover=hover if distinct else hover[:1], spans=_spans(metric),
        thresholds=lines, width="double",
        caption=text.fill(PAGE, f"{words}.caption", n=frame[key].n_unique(),
                          unit=text.word(PAGE, f"figures.band_map.unit.{key}")),
    )  # fmt: skip


def overrepresented_strip(plots: Plots, metric: Metric, units: pl.DataFrame,
                          strata: list[str]) -> StripData | None:  # fmt: skip
    """The top overrepresented sequence's share per unit, by stratum, over FastQC's limits.

    Drawn by ``figures.strip``, the one strip builder, as a numeric metric of its own whose
    bands are the status metric's display bands; nothing judges the share.
    """
    frame = plots.overrepresented
    if frame.is_empty():
        return None
    words = f"figures.fastqc.{metric.id}"
    shown = replace(metric, column="top_pct", kind="numeric", unit="%",
                    computed=text.word(PAGE, f"{words}.computed"),
                    meaning=text.word(PAGE, f"{words}.meaning"))  # fmt: skip
    judged = (pl.DataFrame(schema={"metric": pl.String, "stratum": pl.String,
                                   "median": pl.Float64, "low": pl.Float64,
                                   "high": pl.Float64}),
              pl.DataFrame(schema={"sample": pl.String, "metric": pl.String,
                                   "reason": pl.String}))  # fmt: skip
    return figures.strip(shown, frame, judged, units.select(plots.key, *strata), strata,
                         plots.key)  # fmt: skip
