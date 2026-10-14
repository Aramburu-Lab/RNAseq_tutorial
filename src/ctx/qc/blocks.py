"""The M2 page's blocks: each group of figures, its words from ``content/qc.toml``.

``ctx.qc.run`` measures and judges; this module turns a judgement into the page's groups.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import polars as pl

from ctx.de.run import Writer
from ctx.qc import checkpoint, fastqc, figures
from ctx.qc import registry as reg
from ctx.report import text
from ctx.report.builder import FigureBlock, FigureText, Group
from ctx.viz.contour import DENSITY_MIN_POINTS
from ctx.viz.style import MAX_CATEGORIES, MAX_DISTINCT

if TYPE_CHECKING:
    from ctx.qc.run import Options

PAGE = "qc"
# The metric kinds drawn as a strip: a value per sample or run on one axis.
STRIPPED = ("numeric", "strandedness")


@dataclass(frozen=True)
class Stage:
    """One stage's judgement: its values, bands and flags, all keyed by sample."""

    values: pl.DataFrame
    bands: pl.DataFrame
    flags: pl.DataFrame
    notes: tuple[str, ...] = ()
    run_values: pl.DataFrame | None = None
    run_flags: pl.DataFrame | None = None


@dataclass(frozen=True)
class Reads:
    """What the reads showed: values per run, the runs left unread, FastQC's plots, the grid."""

    values: pl.DataFrame
    skipped: pl.DataFrame
    plots: fastqc.Plots
    grid: pl.DataFrame
    grid_run: str


def fastqc_figures(
    writer: Writer, plots: fastqc.Plots, units: pl.DataFrame, registry: reg.Registry,
    options: "Options",
) -> list[FigureBlock | str]:  # fmt: skip
    """FastQC's own plots behind its verdicts, one tab each, then the note on what is missing.

    Args:
        writer: The step's writer.
        plots: FastQC's plot data, per run or per sample.
        units: One row per unit: ``plots.key``, ``sample`` and the sheet's columns.
        registry: The registry, whose status metrics name the plots and their bands.
        options: The run's options; the strata and peers.

    Returns:
        The figures, then the note on per-base sequence content.
    """
    colour = figures.colour_column(units, options.peers)
    unit = text.word(PAGE, f"figures.band_map.unit.{plots.key}")
    items: list[FigureBlock | str] = []
    for metric_id, what in zip(fastqc.PLOTS, ("base_quality", "fastqc_adapter", "fastqc_n_content",
                                               "gc_distribution"), strict=True):  # fmt: skip
        data = fastqc.curve_figure(plots, registry.metrics[metric_id], units, colour)
        if data is None:
            continue
        styles = (text.fill(PAGE, "figures.fastqc.styles_distinct", unit=unit) if data.distinct
                  else text.fill(PAGE, "figures.fastqc.styles_shared", most=MAX_DISTINCT,
                                 unit=unit))  # fmt: skip
        source = text.fill(PAGE, "figures.fastqc.source",
                           source=registry.metrics[metric_id].plot_source)  # fmt: skip
        items.append(figure(writer, data, what, f"fastqc.{metric_id}", styles, source,
                            n=data.frame[plots.key].n_unique(), unit=unit))  # fmt: skip
    metric = registry.metrics[fastqc.OVERREPRESENTED[0]]
    over = fastqc.overrepresented_strip(plots, metric, units, list(options.strata))
    if over is not None:
        items.append(figure(writer, over, "overrepresented", "fastqc.fastqc_overrepresented",
                            text.fill(PAGE, "figures.fastqc.source", source=metric.plot_source),
                            n=over.frame.height, unit=unit))  # fmt: skip
    if items:
        items.append(text.fill(PAGE, "notes.no_sequence_content",
                               title=registry.metrics[fastqc.NO_TABLE].title))  # fmt: skip
    return items


def read_figures(writer: Writer, seen: Reads, sheet: pl.DataFrame, registry: reg.Registry,
                 options: "Options") -> Group:  # fmt: skip
    """FastQC's plots and the fastp grid, with the runs left unread said first."""
    items: list[FigureBlock | str] = []
    if seen.skipped.height:
        items.append(text.fill(PAGE, "notes.skipped", n=seen.skipped.height,
                               runs=", ".join(seen.skipped["run_accession"].head(10))))  # fmt: skip
    items += fastqc_figures(writer, seen.plots, sheet, registry, options)
    grid = figures.grid_figure(seen.grid, options.fastp, seen.grid_run)
    if grid is None:
        items.append(figures.flat_grid(seen.grid, seen.grid_run))
    else:
        counts, fastp = figures.grid_counts(seen.grid), options.fastp
        floor = (text.fill(PAGE, "figures.grid.reading_floor", **counts),) if counts["none"] else ()
        items.append(figure(writer, grid, "fastp_grid", "grid", *floor,
                            lead="lead_floor" if counts["none"] else "lead", run=seen.grid_run,
                            q=fastp.q, l=fastp.l, bands=len(grid.edges()) - 1,
                            isoline=figures.GRID_ISOLINE, **counts))  # fmt: skip
    run_bands = figures.band_map(seen.values, registry, "run_accession")
    if run_bands is not None:
        items.append(band_map(writer, run_bands, "run_band_map", "run_accession"))
    return Group(text.word(PAGE, "groups.reads"), tuple(items))


def strips(writer: Writer, stage: Stage, units: pl.DataFrame, registry: reg.Registry,
           options: "Options", key: str) -> list[FigureBlock]:  # fmt: skip
    """A strip per numeric or strandedness metric the values hold, one tab each, in order."""
    values = stage.values if key == "sample" else stage.run_values
    flags = stage.flags if key == "sample" else stage.run_flags
    if values is None or flags is None:
        return []
    blocks = []
    for metric in registry.metrics.values():
        if metric.kind not in STRIPPED or metric.column not in values.columns:
            continue
        within = checkpoint.strata_of(metric, options.strata)
        data = figures.strip(metric, values, (stage.bands, flags), units, within, key)
        if data is None:
            continue
        lead = text.fill(PAGE, "figures.strip.lead", flagged=int(data.frame["flagged"].sum()),
                         n=data.frame.height, k=f"{registry.mad_multiplier:g}")  # fmt: skip
        reading = tuple(line for line in (metric.meaning, metric.rationale) if line)
        reading += _equal_values(data.frame["value"], metric) + _mixed_lengths(units, metric)
        if figures.shading_kind(metric) != "none":
            reading += (text.fill(PAGE, "figures.strip.source",
                                  source=metric.plot_source or metric.source),)  # fmt: skip
        text_around = FigureText(title=metric.title, tab=metric.title, lead=lead, reading=reading)
        blocks.append(writer.figure(data, f"strip_{metric.id}", text_around))
    return blocks


def _equal_values(values: pl.Series, metric: reg.Metric) -> tuple[str, ...]:
    """The line saying every unit holds one value, where they do: a cut, not a measurement.

    The course cut each teaching run to the same number of reads, which puts all 8 in
    ENCODE's bad depth band by construction (D-052, proposal decision 5).
    """
    if values.n_unique() != 1 or values.len() < 2:  # noqa: PLR2004 - one unit is not a cut
        return ()
    return (text.fill(PAGE, "figures.strip.equal", n=values.len(), value=figures.cut(values[0]),
                      unit=metric.unit),)  # fmt: skip


def _mixed_lengths(units: pl.DataFrame, metric: reg.Metric) -> tuple[str, ...]:
    """The line saying the cohort mixes read lengths by design, for a metric judged within them."""
    groups = units[checkpoint.LENGTH_GROUP].unique(maintain_order=True).to_list()
    if metric.judged_within != checkpoint.LENGTH_GROUP or len(groups) < 2:  # noqa: PLR2004
        return ()
    return (text.fill(PAGE, "figures.strip.mixed_lengths", n=len(groups),
                      groups="; ".join(str(g) for g in groups)),)  # fmt: skip


def lines(key: str, **values: Any) -> tuple[str, ...]:
    """The list of text at a key, each line's placeholders filled from ``values``.

    Raises:
        TextError: If a line needs a value the code did not give.
    """
    found = text.words(PAGE, key)
    if missing := {n for line in found for n in text.placeholders(line)} - values.keys():
        raise text.TextError(f"content/{PAGE}.toml {key!r} needs {sorted(missing)}")
    return tuple(line.format_map(values) for line in found)


def colour_line(column: str | None) -> str:
    """The reading-guide line naming the sheet column that colours a figure's marks."""
    if column is None:
        return text.fill(PAGE, "colour.none", most=MAX_CATEGORIES)
    return text.fill(PAGE, "colour.by", column=column)


def figure(writer: Writer, data: Any, what: str, key: str, *extra: str, lead: str = "lead",
           **values: Any) -> FigureBlock:  # fmt: skip
    """A page figure whose title, lead and reading guide are ``figures.<key>`` in the content.

    Args:
        writer: The step's writer.
        data: The figure's data contract.
        what: The figure's file description.
        key: Its table under ``figures`` in ``content/qc.toml``.
        extra: Lines added to the reading guide after the content's own.
        lead: The lead's key in that table, where a figure has more than one.
        values: The placeholders of the lead and the reading guide.

    Returns:
        The figure, written.
    """
    return writer.figure(data, what, FigureText(
        title=text.word(PAGE, f"figures.{key}.title"),
        lead=text.fill(PAGE, f"figures.{key}.{lead}", **values),
        reading=(*lines(f"figures.{key}.reading", **values), *extra)))  # fmt: skip


def band_map(writer: Writer, data: Any, what: str, key: str) -> FigureBlock:
    """The band map, its lead counting the bad and usable values, one column per sample or run."""
    bands = data.frame["band"]
    return figure(writer, data, what, "band_map", cells=bands.len(),
                  bad=int((bands == "bad").sum()), usable=int((bands == "usable").sum()),
                  unit=text.word(PAGE, f"figures.band_map.unit.{key}"))  # fmt: skip


def identity_figures(writer: Writer, stage: Stage, measured: checkpoint.Measured,
                     sheet: pl.DataFrame, registry: reg.Registry,
                     options: "Options") -> list[FigureBlock]:  # fmt: skip
    """Library size against genes, the composition scatter and the sex check, where measured."""
    items = []
    colour = figures.colour_column(sheet, options.peers)
    scatter = figures.density(stage.values, sheet, stage.flags, colour, registry)
    if scatter is not None:
        drawn = scatter.frame
        elsewhere = (drawn["flagged_by"] != "") & ~drawn["flagged"]
        items.append(figure(writer, scatter, "size_genes", "density", colour_line(colour),
                            flagged=int(drawn["flagged"].sum()), other=int(elsewhere.sum()),
                            n=drawn.height, minimum=DENSITY_MIN_POINTS, **figures.HDR))  # fmt: skip
    mix = figures.composition(stage.values, sheet, stage.flags, colour, registry)
    if mix is not None:
        items.append(figure(writer, mix, "composition", "composition", colour_line(colour),
                            n=mix.frame.height,
                            flagged=int((mix.frame["flagged_by"] != "").sum())))  # fmt: skip
    if measured.sexes is not None:
        sexes = figures.sex_scatter(measured.sexes)
        disagree = sexes.frame.filter(pl.col("agrees").not_()).height
        items.append(figure(writer, sexes, "sex_check", "sex", n=sexes.frame.height,
                            disagree=disagree))  # fmt: skip
    return items


def overview_groups(writer: Writer, stage: Stage, measured: checkpoint.Measured | None,
                    context: tuple[pl.DataFrame, pl.DataFrame], registry: reg.Registry,
                    options: "Options") -> list[Group]:  # fmt: skip
    """The groups every run draws: the glance, every metric, composition and identity."""
    sheet, status = context
    groups, glance, identity = [], [], []
    bands = figures.band_map(stage.values, registry, "sample")
    if bands is not None:
        glance.append(band_map(writer, bands, "band_map", "sample"))
    if measured is not None:
        dims = (*dict.fromkeys([*options.strata, *options.peers, "batch"]),)
        linked = figures.linked_pca(measured.scores, measured.share, sheet, status, dims)
        glance.append(figure(writer, linked, "linked_pca", "linked"))
        identity = identity_figures(writer, stage, measured, sheet, registry, options)
    units = sheet.unique("sample", keep="first").select(
        "sample", *options.strata, checkpoint.LENGTH_GROUP
    )
    metrics = strips(writer, stage, units, registry, options, "sample")
    for title, items, tabs in (("glance", glance, False), ("metrics", metrics, True),
                               ("identity", identity, False)):  # fmt: skip
        if items:
            groups.append(Group(text.word(PAGE, f"groups.{title}"), tuple(items),
                                tabs=tabs and len(items) > 1))  # fmt: skip
    return groups


def fate_figures(writer: Writer, fate: pl.DataFrame, values: pl.DataFrame, sheet: pl.DataFrame,
                 registry: reg.Registry, options: "Options") -> Group:  # fmt: skip
    """The read-fates bubble chart and the Sankey, as one group of the page."""
    colour = figures.colour_column(sheet, options.peers)
    bubble = figures.fates_bubble(fate, values, sheet, colour, registry)
    median = float(bubble.frame["useful_share"].median())  # type: ignore[arg-type]
    total, useful = int(fate["total_reads"].sum()), int(fate["num_mapped"].sum())
    blocks = (
        figure(writer, bubble, "read_fates", "fates", colour_line(colour), median=median,
               k=f"{registry.mad_multiplier:g}"),
        figure(writer, figures.fates_flow(fate), "read_flow", "flow", total=total, useful=useful,
               share=useful / total if total else 0.0),
    )  # fmt: skip
    return Group(text.word(PAGE, "groups.fates"), blocks)
