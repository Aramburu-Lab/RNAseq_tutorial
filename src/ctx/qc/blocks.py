"""The M2 page's blocks: each group of figures, its words from ``content/qc.toml``.

``ctx.qc.run`` measures and judges; this module turns a judgement into the page's groups.
"""

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

import polars as pl

from ctx.de.run import Writer
from ctx.qc import checkpoint, fastqc, figures, flag_map, identity, judge
from ctx.qc import registry as reg
from ctx.report import export, text
from ctx.report.builder import FigureBlock, FigureText, Group
from ctx.report.views import viewed
from ctx.viz.contour import DENSITY_MIN_POINTS
from ctx.viz.style import MAX_CATEGORIES, MAX_DISTINCT, figure_inches

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
    by_peers: tuple[pl.DataFrame, ...] = ()


@dataclass(frozen=True)
class Reads:
    """What the reads showed: values per run, the runs left unread, FastQC's plots, the grid.

    ``source`` says where the tools ran when it was not here: empty for the local run.
    """

    values: pl.DataFrame
    skipped: pl.DataFrame
    plots: fastqc.Plots
    grid: pl.DataFrame
    grid_run: str
    source: str = ""


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
    items: list[FigureBlock | str] = [seen.source] if seen.source else []
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
    """A strip per numeric or strandedness metric the values hold, one tab each, in order.

    Each carries the judge panel (``judge.settings``): ``units`` holds every grouping column
    the panel offers, and a moved panel's PDF and SVG draw the strip as rejudged.
    """
    options_offered = checkpoint.groupings(options.strata, options.peers, units.columns)
    values = stage.values if key == "sample" else stage.run_values
    flags = stage.flags if key == "sample" else stage.run_flags
    if values is None or flags is None:
        return []
    blocks: list[FigureBlock] = []
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
        block = writer.figure(data, f"strip_{metric.id}", text_around)
        panel = judge.settings(metric, data, values=values.select(key, metric.column),
                               units=units, key=key, registry=registry, options=options_offered,
                               by_peers=stage.by_peers)  # fmt: skip
        if any(b.controls and b.controls.get("units") for b in blocks):
            # One scope's strips share their units, which the page reads from the first:
            # 20 copies were 407,060 of Step-1's 612,692 bytes of panels.
            panel = {k: v for k, v in panel.items() if k != "units"}
        ready = export.settings(figure_inches(data.width, "strip"), data.caption, block.stem)
        blocks.append(replace(block, controls={**panel, "export": ready}))
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


def grouping_columns(sheet: pl.DataFrame, options: "Options") -> list[str]:
    """Every sheet column a strata or peer grouping of the judge panel names, each once."""
    strata, peers = checkpoint.groupings(options.strata, options.peers, sheet.columns)
    return [c for c in dict.fromkeys(c for o in (*strata, *peers) for c in o)
            if c != checkpoint.LENGTH_GROUP]  # fmt: skip


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


def band_map(writer: Writer, data: Any, what: str, key: str,
             flags: Any = None) -> FigureBlock:  # fmt: skip
    """The band map, its lead counting the bad and usable values, one column per sample or run.

    With ``flags``, the sample map's flag view (``flag_map``), each view with its own PDF and
    SVG (``exported``).
    """
    bands = data.frame["band"]
    values = {"cells": bands.len(), "bad": int((bands == "bad").sum()),
              "usable": int((bands == "usable").sum()),
              "unit": text.word(PAGE, f"figures.band_map.unit.{key}")}  # fmt: skip
    block = figure(writer, data, what, "band_map", **values)
    if flags is None:
        return block
    found = flag_map.counts(flags)
    options = [
        (text.word(PAGE, "figures.flag_map.view_bands"), writer.marked(data), FigureText(
            title=text.word(PAGE, "figures.band_map.title"),
            lead=text.fill(PAGE, "figures.band_map.lead", **values),
            reading=lines("figures.band_map.reading", **values))),
        (text.word(PAGE, "figures.flag_map.view_flags"), writer.marked(flags), FigureText(
            title=text.word(PAGE, "figures.flag_map.title"),
            lead=text.fill(PAGE, "figures.flag_map.lead", **found),
            reading=lines("figures.flag_map.reading"))),
    ]  # fmt: skip
    return exported(block, options, text.word(PAGE, "figures.flag_map.view"))


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
        items.append(sex_figure(writer, measured.sexes, sheet))
    return items


def own_sex(sexes: pl.DataFrame, sheet: pl.DataFrame) -> pl.DataFrame | None:
    """The sex check against each sample's own recorded sex; None where no donor sex differs.

    The method compares a transplanted sample with its donor (``identity.sex_check``); this is
    the rule it avoids, through the same agreement rule (``identity.against``).
    """
    if "donor_sex" not in sheet.columns or not (sheet["donor_sex"].fill_null("") != "").any():
        return None
    recorded = sexes.select("sample").join(
        sheet.unique("sample", keep="first", maintain_order=True).select("sample", "sex"),
        on="sample", how="left", maintain_order="left")  # fmt: skip
    return identity.against(sexes, recorded["sex"].fill_null("").to_list())


def sex_figure(writer: Writer, sexes: pl.DataFrame, sheet: pl.DataFrame) -> FigureBlock:
    """Recorded sex against expression; with donor sexes, a view against each sample's own.

    Work-order item 5 (D-051): by the method 0 of Step-1's 208 compared samples disagree,
    against their own sex 31 (figure controls plan), so the select shows why the method
    compares a transplanted sample with its donor. Both views' calls are one table.
    """
    shown = [(text.word(PAGE, "figures.sex.view_donor"), sexes, "caption")]
    if (own := own_sex(sexes, sheet)) is not None:
        shown.append((text.word(PAGE, "figures.sex.view_own"), own, "caption_own"))
    options = []
    for label, frame, caption in shown:
        data = figures.sex_scatter(frame, caption)
        words = FigureText(
            title=text.word(PAGE, "figures.sex.title"),
            lead=text.fill(PAGE, "figures.sex.lead", n=data.frame["agrees"].count(),
                           disagree=data.frame.filter(pl.col("agrees").not_()).height),
            reading=lines("figures.sex.reading"))  # fmt: skip
        options.append((label, data, words))
    block = writer.figure(options[0][1], "sex_check", options[0][2])
    if len(options) < 2:  # noqa: PLR2004 - a view needs a second option
        return block
    views = pl.concat([f.select("sample", "expected", "expression_sex", "agrees")
                       .with_columns(view=pl.lit(label)) for label, f, _ in shown])  # fmt: skip
    table(writer, views.select("view", "sample", "expected", "expression_sex", "agrees"),
          "sex_views")  # fmt: skip
    return exported(block, options, text.word(PAGE, "figures.sex.view"))


def exported(block: FigureBlock, options: list[tuple[str, Any, FigureText]], label: str,
             inches: tuple[float, float] | None = None) -> FigureBlock:  # fmt: skip
    """``viewed``, each option carrying the settings of its own PDF and SVG.

    The view kind exports the shown figure with the build's caption and labels (D-073): M2's
    own-sex view printed "A transplanted sample is compared with its donor" under 31 samples
    compared with their own sex. Each option here takes its own caption, size and layers;
    ``inches`` sizes a contract that names no size of its own, as FlowData names none.
    """
    shown = viewed(block, options, label)
    controls = dict(shown.controls or {})

    def settings(data: Any) -> dict[str, Any]:
        """One view's PDF and SVG settings: its own caption, size and layers."""
        layers = data.export_layers() if hasattr(data, "export_layers") else None
        return export.settings(inches or data.inches(), data.caption, block.stem, layers)

    each = [settings(data) for _, data, _ in options]
    controls["options"] = [{**o, "export": e} for o, e in zip(controls["options"], each,
                                                               strict=True)]  # fmt: skip
    controls["export"] = each[0]
    return replace(shown, controls=controls)


def table(writer: Writer, frame: pl.DataFrame, what: str) -> None:
    """A table the page's views read, each column described by ``tables.<what>`` in the content."""
    units = text.raw(PAGE, f"tables.{what}.units")
    words = {c: text.word(PAGE, f"tables.{what}.columns.{c}") for c in frame.columns}
    writer.table(frame, what, {c: (units.get(c, "none"), words[c]) for c in frame.columns})


def overview_groups(writer: Writer, stage: Stage, measured: checkpoint.Measured | None,
                    context: tuple[pl.DataFrame, pl.DataFrame], registry: reg.Registry,
                    options: "Options") -> list[Group]:  # fmt: skip
    """The groups every run draws: the glance, every metric, composition and identity."""
    sheet, status = context
    groups, glance, identity = [], [], []
    units = sheet.unique("sample", keep="first", maintain_order=True).select(
        "sample", *grouping_columns(sheet, options), checkpoint.LENGTH_GROUP
    )
    bands = figures.band_map(stage.values, registry, "sample")
    if bands is not None:
        flags = flag_map.flag_map(stage.values, (stage.bands, stage.flags), units, registry,
                                  tuple(options.strata))  # fmt: skip
        glance.append(band_map(writer, bands, "band_map", "sample", flags))
    if measured is not None:
        dims = (*dict.fromkeys([*options.strata, *options.peers, "batch"]),)
        linked = figures.linked_pca(measured.scores, measured.share, sheet, status, dims)
        glance.append(figure(writer, linked, "linked_pca", "linked"))
        identity = identity_figures(writer, stage, measured, sheet, registry, options)
    peered = replace(stage, by_peers=measured.by_peers) if measured is not None else stage
    metrics = strips(writer, peered, units, registry, options, "sample")
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
    blocks = (
        figure(writer, bubble, "read_fates", "fates", colour_line(colour), median=median,
               k=f"{registry.mad_multiplier:g}"),
        flow_figure(writer, fate, sheet, (*dict.fromkeys([*options.strata, *options.peers,
                                                         "batch"]),)),
    )  # fmt: skip
    return Group(text.word(PAGE, "groups.fates"), blocks)


def flow_words(fate: pl.DataFrame, view: str = "") -> FigureText:
    """The Sankey's title, lead and reading guide, its numbers summed over ``fate``'s samples.

    A view other than all samples is named in the title, so a switch of view shows on the figure.
    """
    total, useful = int(fate["total_reads"].sum()), int(fate["num_mapped"].sum())
    title = (
        text.fill(PAGE, "figures.flow.title_view", view=view)
        if view
        else text.word(PAGE, "figures.flow.title")
    )
    return FigureText(
        title=title,
        lead=text.fill(PAGE, "figures.flow.lead", total=total, useful=useful,
                       share=useful / total if total else 0.0),
        reading=lines("figures.flow.reading"))  # fmt: skip


def flow_figure(writer: Writer, fate: pl.DataFrame, sheet: pl.DataFrame,
                dims: tuple[str, ...]) -> FigureBlock:  # fmt: skip
    """Where every read goes, over all samples; a view per level of each of ``dims``.

    Work-order item 5 (D-051): the per-sample median decoy share is 0.260 in BM_4 against
    0.150 in LDN on Step-1 (figure controls plan), so where reads go depends on the cells. A
    level's links are its samples' fates summed by the same ``fates_flow``; every view's links
    are one table. A column with more levels than a page colours, or only one, offers none.
    """
    samples = sheet.unique("sample", keep="first", maintain_order=True)
    shown = [(text.word(PAGE, "figures.flow.view_all"), fate)]
    for column in (d for d in dims if d in samples.columns):
        levels = samples[column].cast(pl.String).fill_null("").unique(maintain_order=True)
        if not 1 < levels.len() <= MAX_CATEGORIES:
            continue
        for level in levels.sort().filter(levels.sort() != ""):
            kept = samples.filter(pl.col(column).cast(pl.String).fill_null("") == level)["sample"]
            part = fate.filter(pl.col("sample").is_in(kept.implode()))
            if part.height:
                label = text.fill(PAGE, "figures.flow.view_level", column=column, level=level)
                shown.append((label, part))
    options = [
        (label, figures.fates_flow(part), flow_words(part, label if i else ""))
        for i, (label, part) in enumerate(shown)
    ]
    block = writer.figure(options[0][1], "read_flow", options[0][2])
    if len(options) < 2:  # noqa: PLR2004 - a view needs a second option
        return block
    links = pl.concat([data.links.with_columns(view=pl.lit(label)) for label, data, _ in options])
    table(writer, links.select("view", "source", "target", "value"), "flow_views")
    return exported(block, options, text.word(PAGE, "figures.flow.view"),
                    figure_inches(options[0][1].width, "flow"))  # fmt: skip
