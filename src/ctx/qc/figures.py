"""The QC page's figures, built from the run's tables; every word from ``content/qc.toml``."""

from collections.abc import Sequence
from dataclasses import replace
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import median_abs_deviation

from ctx.qc import Fastp
from ctx.qc.identity import SEX_X, SEX_Y
from ctx.qc.registry import STATUS_BANDS as FASTQC_BANDS
from ctx.qc.registry import Metric, Registry
from ctx.report import text
from ctx.viz.contour import DENSITY_MIN_POINTS, HDR_MASSES, Mark, SurfaceData, above_floor
from ctx.viz.density import DensityData
from ctx.viz.flow import FlowData, Node
from ctx.viz.heatmap import HeatmapData, Swatch
from ctx.viz.linked import LinkedData
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.strip import Range, Span, StripData
from ctx.viz.style import BAND_COLOURS, MAX_CATEGORIES, ROLES, category_colours

PAGE = "qc"
# The composition shares a bubble's area adds up, where the counts had gene types.
SHARES = ("mt_fraction", "mt_rrna_fraction", "rrna_fraction")
# Trim Galore's fates and Salmon's, in the order the Sankey stacks them.
TRIM_FATES = ("reads_too_short", "reads_too_long", "reads_too_many_n", "reads_discarded_untrimmed")
SALMON_FATES = ("num_mapped", "num_decoy_fragments", "not_assigned")
# The two highest-density regions a density draws, as its caption and reading guide name them.
HDR = {"inner": f"{HDR_MASSES[0] * 100:.0f} %", "outer": f"{HDR_MASSES[1] * 100:.0f} %"}


def _band(useful: np.ndarray, k: float) -> tuple[float, float, float]:
    """The useful fraction's median and median ± k scaled MAD, clipped to 0 and 1."""
    median = float(np.median(useful))
    mad = float(median_abs_deviation(useful, scale="normal"))
    return median, max(median - k * mad, 0.0), min(median + k * mad, 1.0)


def fates_bubble(
    fates: pl.DataFrame, values: pl.DataFrame, sheet: pl.DataFrame, colour: str | None,
    registry: Registry,
) -> ScatterData:  # fmt: skip
    """Raw against useful fragments, log axes, one bubble per sample.

    The area is the MT + rRNA share of counted reads where the counts had gene types; the
    lines are y = x and the useful fraction's median ± k scaled MAD across the samples.
    """
    shares = [c for c in SHARES if c in values.columns]
    samples = sheet.unique("sample", keep="first")
    if colour is None:
        colour, samples = "group", samples.with_columns(
            group=pl.lit(text.word(PAGE, "figures.fates.all")))  # fmt: skip
    frame = fates.join(samples.select("sample", colour), on="sample", how="left",
                       maintain_order="left").join(
        values.select("sample", *shares), on="sample", how="left",
        maintain_order="left")  # fmt: skip
    frame = frame.with_columns(useful_share=pl.col("num_mapped") / pl.col("total_reads"))
    if shares:
        frame = frame.with_columns(mt_rrna_share=pl.sum_horizontal(shares))
    levels = [str(v) for v in frame[colour].unique(maintain_order=True)]
    _, low, high = _band(frame["useful_share"].to_numpy(), registry.mad_multiplier)
    raw = frame["total_reads"].to_numpy().astype(float)
    ends = (float(raw.min()) * 0.9, float(raw.max()) * 1.1)
    lines = (Line(text.word(PAGE, "figures.fates.identity"), ends, ends),
             *(Line(text.fill(PAGE, "figures.fates.band", k=f"{registry.mad_multiplier:g}"), ends,
                    (f * ends[0], f * ends[1]), dash=True) for f in (low, high)))  # fmt: skip
    units = text.raw(PAGE, "figures.fates.units")
    columns = {c: (units.get(c, "none"), text.fill(PAGE, "figures.fates.colour", column=c)
                   if c == colour else text.word(PAGE, f"figures.fates.columns.{c}"))
               for c in frame.columns}  # fmt: skip
    area = text.word(PAGE, "figures.fates.area" if shares else "figures.fates.no_area")
    caption = text.fill(PAGE, "figures.fates.caption", n=frame.height,
                        k=f"{registry.mad_multiplier:g}", area=area)  # fmt: skip
    return ScatterData(
        frame=frame.with_columns(pl.col(colour).cast(pl.String)), columns=columns,
        x="total_reads", y="num_mapped", x_label=text.word(PAGE, "figures.fates.x"),
        y_label=text.word(PAGE, "figures.fates.y"), colour=colour,
        colours=category_colours(levels), caption=caption, lines=lines, log_x=True, log_y=True,
        hover=("sample", "useful_share", *TRIM_FATES, "num_decoy_fragments", "not_assigned"),
        width="double", size="mt_rrna_share" if shares else None,
        size_label=text.word(PAGE, "figures.fates.size"),
    )  # fmt: skip


def fates_flow(fates: pl.DataFrame) -> FlowData:
    """Every read of the samples, summed: sequenced, Trim Galore's fates, then Salmon's."""
    total = {c: float(fates[c].sum()) for c in (*fates.columns[1:],)}
    word = {c: text.word(PAGE, f"figures.flow.nodes.{c}") for c in
            ("total_reads", "reads_written", *TRIM_FATES, *SALMON_FATES)}  # fmt: skip
    kept, lost = BAND_COLOURS["good"], None
    nodes = [Node("total_reads", word["total_reads"], 0, ROLES["reference"]),
             Node("reads_written", word["reads_written"], 1, kept),
             *(Node(c, word[c], 1, lost) for c in TRIM_FATES if total[c] > 0),
             Node("num_mapped", word["num_mapped"], 2, kept),
             *(Node(c, word[c], 2, lost) for c in SALMON_FATES[1:] if total[c] > 0)]  # fmt: skip
    links = [("total_reads", c, total[c]) for c in ("reads_written", *TRIM_FATES) if total[c] > 0]
    links += [("reads_written", c, total[c]) for c in SALMON_FATES if total[c] > 0]
    return FlowData(
        nodes=tuple(nodes),
        links=pl.DataFrame(links, schema=["source", "target", "value"], orient="row"),
        columns={
            c: ("reads" if c == "value" else "none", text.word(PAGE, f"figures.flow.columns.{c}"))
            for c in ("source", "target", "value")
        },  # fmt: skip
        stage_labels=tuple(text.words(PAGE, "figures.flow.stages")),
        caption=text.fill(PAGE, "figures.flow.caption", n=fates.height),
        width="double",
    )


def colour_column(sheet: pl.DataFrame, candidates: Sequence[str]) -> str | None:
    """The first candidate column with few enough levels to colour by; None if none has."""
    fits = [c for c in candidates if c in sheet.columns and sheet[c].n_unique() <= MAX_CATEGORIES]
    return fits[0] if fits else None


# The isoline a grid map draws: the settings that still keep this share of the reads.
GRID_ISOLINE = 90.0
# A grid cell keeping this share keeps no read, as -l past the reads' length does
# (``reads.GRID_L``); it is drawn as the floor, and the bands span the cells above it.
GRID_FLOOR = 0.0


def grid_counts(grid: pl.DataFrame) -> dict[str, Any]:
    """The grid's cells, those that keep reads and those on the floor, and what they keep."""
    kept = grid["passed_pct"].filter(above_floor(grid["passed_pct"].to_numpy(), GRID_FLOOR))
    return {"cells": grid.height, "kept": kept.len(), "none": grid.height - kept.len(),
            "low": kept.min(), "high": kept.max()}  # fmt: skip


def grid_figure(grid: pl.DataFrame, fastp: Fastp, run: str) -> SurfaceData | None:
    """Reads fastp keeps across -q and -l on one run, settings that keep none as the floor.

    None when the cells that keep reads all keep the same share: there is no band to draw.
    """
    counts = grid_counts(grid)
    if counts["low"] is None or counts["low"] == counts["high"]:
        return None
    units = text.raw(PAGE, "figures.grid.units")
    surface = SurfaceData(
        frame=grid, x="q", y="l", z="passed_pct",
        columns={c: (units.get(c, "none"), text.word(PAGE, f"figures.grid.columns.{c}"))
                 for c in grid.columns},  # fmt: skip
        x_label=text.word(PAGE, "figures.grid.x"), y_label=text.word(PAGE, "figures.grid.y"),
        z_label=text.word(PAGE, "figures.grid.z"), caption="",
        isolines=(GRID_ISOLINE,), isoline_format=".0f",
        mark=Mark(float(fastp.q), float(fastp.l), text.fill(PAGE, "figures.grid.mark", q=fastp.q,
                                                              l=fastp.l)),
        floor=GRID_FLOOR, floor_label=text.word(PAGE, "figures.grid.floor"),
    )  # fmt: skip
    # The band count is the one the surface draws, so the caption cannot name another.
    key = "figures.grid.caption_floor" if counts["none"] else "figures.grid.caption"
    return replace(surface, caption=text.fill(
        PAGE, key, run=run, u=fastp.u, nq=grid["q"].n_unique(), nl=grid["l"].n_unique(),
        isoline=GRID_ISOLINE, bands=len(surface.edges()) - 1, **counts))  # fmt: skip


def flat_grid(grid: pl.DataFrame, run: str) -> str:
    """The note in place of a grid with nothing to map.

    Every cell keeps the same share, or every cell above the floor does and the rest none.
    """
    counts = grid_counts(grid)
    if not counts["none"] or not counts["kept"]:
        return text.fill(PAGE, "notes.flat_grid", run=run, pct=float(grid["passed_pct"][0]))
    return text.fill(PAGE, "notes.flat_floor", run=run, pct=float(counts["low"]), **counts)


def stratum_labels(strata: Sequence[str]) -> pl.Expr:
    """Each row's stratum as the judge names it: its strata values joined by " / "."""
    parts = [pl.col(c).cast(pl.String).fill_null("None") for c in strata]
    return pl.concat_str(parts, separator=" / ").alias("stratum")


def strip(
    metric: Metric, values: pl.DataFrame, judged: tuple[pl.DataFrame, pl.DataFrame],
    units: pl.DataFrame, strata: Sequence[str], key: str,
) -> StripData | None:  # fmt: skip
    """One metric's values by stratum: flagged ones as triangles, its bands behind them.

    Args:
        metric: The registry metric.
        values: ``key`` and the metric's column.
        judged: The bands and the flags tables of the same judgement.
        units: ``key`` and the strata columns of every sample or run.
        strata: The strata columns.
        key: ``sample`` or ``run_accession``.

    Returns:
        The strip, or None where no unit has a value.
    """
    bands, flags = judged
    frame = (
        values.select(key, metric.column)
        .drop_nulls()
        .join(units, on=key, how="left", maintain_order="left")
    )
    if frame.is_empty():
        return None
    frame = frame.with_columns(stratum_labels(strata))
    reasons = dict(flags.filter(pl.col("metric") == metric.id).select("sample", "reason").rows())
    frame = frame.select(key, "stratum", pl.col(metric.column).cast(pl.Float64).alias("value"),
                         flagged=pl.col(key).is_in(list(reasons)),
                         reason=pl.col(key).replace_strict(reasons, default=""))  # fmt: skip
    present = set(frame["stratum"])
    own = bands.filter((pl.col("metric") == metric.id) & pl.col("median").is_not_null())
    rows = own.select("stratum", "median", "low", "high").rows()
    ranges = tuple(Range(*row) for row in rows if row[0] in present)
    spans, references = strip_shading(metric)
    log = metric.scale == "log10" and bool((frame["value"] > 0).all())
    return StripData(
        frame=frame.with_columns(group=pl.lit(text.word(PAGE, "figures.strip.all"))),
        columns={key: ("none", text.word(PAGE, f"columns.{key}")),
                 "stratum": ("none", text.word(PAGE, "figures.strip.stratum")),
                 "value": (metric.unit, metric.meaning),
                 "flagged": ("none", text.word(PAGE, "figures.strip.flagged")),
                 "reason": ("none", text.word(PAGE, "figures.strip.reason")),
                 "group": ("none", text.word(PAGE, "figures.strip.group"))},
        category="stratum", value="value", category_label=text.word(PAGE, "figures.strip.x"),
        value_label=metric.title, colour="group",
        colours={text.word(PAGE, "figures.strip.all"): ROLES["reference"]},
        caption=text.fill(PAGE, "figures.strip.caption", computed=metric.computed,
                          bands=text.word(PAGE, f"figures.strip.shading.{shading_kind(metric)}")),
        flagged="flagged", flag_label=text.word(PAGE, "figures.strip.flagged_label"),
        hover=(key, "reason"), spans=spans, ranges=ranges,
        range_label=text.word(PAGE, "figures.strip.range"), log=log, references=references,
    )  # fmt: skip


def shading_kind(metric: Metric) -> str:
    """What a metric's figures shade, as the caption names it.

    judged: the bands ``band()`` reads; display: ``plot_bands``, which judge nothing; line: a
    reference line alone; calls: a strand call's regions; none: nothing.
    """
    if metric.kind == "strandedness":
        return "calls"
    if metric.banded:
        return "judged"
    if metric.plot_bands is not None:
        return "display"
    return "line" if metric.plot_line is not None else "none"


def strip_shading(metric: Metric) -> tuple[tuple[Span, ...], tuple[tuple[str, float], ...]]:
    """The spans and reference lines a metric's strip draws, all read from the registry.

    The judged bands where it has them, else its display ``plot_bands``; a strandedness
    metric's call regions and Salmon's unstranded strip; ``plot_line`` as a dashed line.
    """
    references = ((text.fill(PAGE, "figures.strip.line", value=cut(metric.plot_line)),
                   metric.plot_line),) if metric.plot_line is not None else ()  # fmt: skip
    if metric.kind == "strandedness":
        calls = text.word(PAGE, "figures.strip.calls")
        spans = [Span("good", calls, low, high, ROLES["reference"])
                 for _, low, high in metric.call_regions()]  # fmt: skip
        if (tolerance := metric.plot_unstranded_tolerance) is not None:
            low, high = 0.5 - tolerance, 0.5 + tolerance
            label = text.fill(PAGE, "figures.strip.tolerance", low=cut(low), high=cut(high))
            spans.append(Span("good", label, low, high, ROLES["reference"]))
        return tuple(spans), references
    word = {b: text.word(PAGE, f"figures.strip.{b}") for b in ("good", "usable", "bad")}
    ranges = metric.spans() if metric.banded else metric.plot_spans()
    return tuple(Span(b, word[b], low, high) for b, low, high in ranges), references


def band_label(metric: Metric) -> tuple[str, ...]:
    """A band map's row in its lines: the metric's title, then its unit and band cut-offs.

    The rows had titles alone, so the map gave no unit and no cut-off (Rytis 2026-10-01). The
    table and hover join the lines with a space; the axes print them one under the other.
    """
    if metric.kind == "status":
        bands = ", ".join(text.fill(PAGE, "figures.band_map.verdict", verdict=v, band=b)
                          for v, b in FASTQC_BANDS.items())  # fmt: skip
    else:
        # Metric.band() reads every direction but "higher" as lower. One cut-off where the
        # source gives one: no value is usable, the rest is bad.
        side = "higher" if metric.direction == "higher" else "lower"
        key = side + ("" if metric.good == metric.usable else "_usable")
        bands = text.fill(PAGE, f"figures.band_map.cuts.{key}", good=cut(metric.good),
                          usable=cut(metric.usable))  # fmt: skip
    return tuple(line.format(title=metric.title, unit=metric.unit, bands=bands)
                 for line in text.words(PAGE, "figures.band_map.row"))  # fmt: skip


def cut(value: float | None) -> str:
    """A band cut-off as a reader reads it: thousands grouped, no exponent, no trailing zero."""
    return f"{value:,.10g}"


def band_map(values: pl.DataFrame, registry: Registry, key: str) -> HeatmapData | None:
    """Every unit's absolute band on every metric that has one; None where none has.

    A column per sample or run, in name order; a row per metric, labelled with its unit and
    its cut-offs.
    """
    metrics = [m for m in registry.metrics.values()
               if m.column in values.columns and (m.banded or m.kind == "status")]  # fmt: skip
    if not metrics:
        return None
    none = text.word(PAGE, "figures.band_map.none")
    lines = {m.id: band_label(m) for m in metrics}
    rows = [(unit, " ".join(lines[m.id]), m.band(row[m.column]) or none)
            for row in values.iter_rows(named=True) for m in metrics
            for unit in (row[key],)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=[key, "metric", "band"], orient="row").with_columns(
        value=pl.lit(None, pl.Float64))  # fmt: skip
    order = tuple(sorted(str(u) for u in values[key].unique()))
    unit = text.word(PAGE, f"figures.band_map.unit.{key}")
    return HeatmapData(
        frame=frame, row="metric", column=key, value="value", level="band",
        columns={key: ("none", text.word(PAGE, f"columns.{key}")),
                 "metric": ("none", text.word(PAGE, "figures.band_map.metric")),
                 "band": ("none", text.word(PAGE, "figures.band_map.band")),
                 "value": ("none", text.word(PAGE, "figures.band_map.value"))},
        row_label=text.word(PAGE, "figures.band_map.row_label"),
        column_label=text.word(PAGE, f"figures.band_map.hover.{key}"),
        column_title=text.fill(PAGE, "figures.band_map.x", unit=unit, n=len(order)),
        column_order=order, row_display={" ".join(v): "\n".join(v) for v in lines.values()},
        swatches=(*(Swatch(b, BAND_COLOURS[b]) for b in ("good", "usable", "bad")),
                  Swatch(none, None)),
        caption=text.fill(PAGE, "figures.band_map.caption", n=values.height, m=len(metrics)),
        width="double",
    )  # fmt: skip


# The value columns library size against genes detected plots, across and up.
SIZE_GENES = ("library_size", "genes_detected")


def flagged_by(flags: pl.DataFrame) -> pl.DataFrame:
    """Each flagged sample and every metric that flagged it, in the flags' order, joined by ", "."""
    return flags.group_by("sample", maintain_order=True).agg(
        flagged_by=pl.col("metric").unique(maintain_order=True).str.join(", "))  # fmt: skip


def density(values: pl.DataFrame, sheet: pl.DataFrame, flags: pl.DataFrame,
            colour: str | None, registry: Registry) -> DensityData | None:  # fmt: skip
    """Library size against genes detected, a density per group of 30 or more samples.

    The triangles are the samples flagged on a plotted metric: marking any flag pointed the
    lead at the dashed line while 102 of Step-1's 107 flagged samples were flagged on neither
    axis (2026-10-07). Every flag stays in ``flagged_by``, which the hover and the CSV show.
    """
    if not set(SIZE_GENES) <= set(values.columns):
        return None
    plotted = [m.id for m in registry.metrics.values() if m.column in SIZE_GENES]
    samples = sheet.unique("sample", keep="first", maintain_order=True).select(
        "sample", *([colour] if colour else [])
    )
    frame = values.select("sample", *SIZE_GENES).join(samples, on="sample", maintain_order="left")
    frame = frame.filter(pl.col("library_size") > 0)
    group = colour or "group"
    if colour is None:
        frame = frame.with_columns(group=pl.lit(text.word(PAGE, "figures.fates.all")))
    on_axes = flags.filter(pl.col("metric").is_in(plotted))["sample"].implode()
    frame = frame.join(flagged_by(flags), on="sample", how="left",
                       maintain_order="left").with_columns(
        pl.col(group).cast(pl.String), flagged=pl.col("sample").is_in(on_axes),
        flagged_by=pl.col("flagged_by").fill_null(""))  # fmt: skip
    levels = [str(v) for v in frame[group].unique(maintain_order=True)]
    data = DensityData(
        frame=frame, x="library_size", y="genes_detected", group=group,
        x_label=text.word(PAGE, "figures.density.x"), y_label=text.word(PAGE, "figures.density.y"),
        columns={"sample": ("none", text.word(PAGE, "columns.sample")),
                 "library_size": ("reads", text.word(PAGE, "figures.density.library_size")),
                 "genes_detected": ("genes", text.word(PAGE, "figures.density.genes_detected")),
                 group: ("none", text.fill(PAGE, "figures.fates.colour", column=group)),
                 "flagged": ("none", text.word(PAGE, "figures.density.flagged")),
                 "flagged_by": ("none", text.word(PAGE, "figures.density.flagged_by"))},
        colours=category_colours(levels), caption="", highlight="flagged",
        highlight_label=text.word(PAGE, "figures.density.flagged_label"),
        hover=("sample", "flagged_by"), log_x=True, width="double",
    )  # fmt: skip
    n = ", ".join(f"{g} {d.n}" for g, d in data.densities.items()) or text.word(
        PAGE, "figures.density.none")  # fmt: skip
    return replace(data, caption=text.fill(PAGE, "figures.density.caption", n=n,
                                           bandwidth=data.bandwidth(),
                                           minimum=DENSITY_MIN_POINTS, **HDR))  # fmt: skip


# The composition shares the scatter plots, across and up (P2Q.3).
COMPOSITION = ("mt_fraction", "rrna_fraction")


def composition(values: pl.DataFrame, sheet: pl.DataFrame, flags: pl.DataFrame,
                colour: str | None, registry: Registry) -> ScatterData | None:  # fmt: skip
    """Mitochondrial against nuclear rRNA share per sample, in %; None without gene types.

    A triangle is a sample flagged on either share; the hover names every metric that
    flagged it, as on the library-size density.
    """
    if not set(COMPOSITION) <= set(values.columns):
        return None
    plotted = [m.id for m in registry.metrics.values() if m.column in COMPOSITION]
    samples = sheet.unique("sample", keep="first", maintain_order=True).select(
        "sample", *([colour] if colour else [])
    )
    group = colour or "group"
    percent = values.select("sample", *(pl.col(c) * 100 for c in COMPOSITION))
    frame = percent.join(samples, on="sample", maintain_order="left")
    if colour is None:
        frame = frame.with_columns(group=pl.lit(text.word(PAGE, "figures.fates.all")))
    on_axes = flags.filter(pl.col("metric").is_in(plotted))["sample"].implode()
    shape = {True: text.word(PAGE, "figures.composition.flagged"),
             False: text.word(PAGE, "figures.composition.not_flagged")}  # fmt: skip
    frame = frame.join(flagged_by(flags), on="sample", how="left",
                       maintain_order="left").with_columns(
        pl.col(group).cast(pl.String), flagged_by=pl.col("flagged_by").fill_null(""),
        flag=pl.col("sample").is_in(on_axes).replace_strict(shape,
                                                            return_dtype=pl.String))  # fmt: skip
    levels = [str(v) for v in frame[group].unique(maintain_order=True)]
    words = "figures.composition"
    return ScatterData(
        frame=frame, x=COMPOSITION[1], y=COMPOSITION[0],
        columns={"sample": ("none", text.word(PAGE, "columns.sample")),
                 **{c: ("%", text.word(PAGE, f"{words}.{c}")) for c in COMPOSITION},
                 group: ("none", text.fill(PAGE, "figures.fates.colour", column=group)),
                 "flagged_by": ("none", text.word(PAGE, "figures.density.flagged_by")),
                 "flag": ("none", text.word(PAGE, f"{words}.flag"))},
        x_label=text.word(PAGE, f"{words}.x"), y_label=text.word(PAGE, f"{words}.y"),
        colour=group, colours=category_colours(levels), symbol="flag",
        caption=text.fill(PAGE, f"{words}.caption", n=frame.height),
        hover=("sample", "flagged_by"), width="double",
    )  # fmt: skip


def cross_sex(sexes: pl.DataFrame, n: int = 2) -> list[str]:
    """The samples whose expression carries the most of the other sex's marker, strongest first.

    The weaker of XIST and the chrY mean, in log2 CPM: high only where both are expressed, as in
    a transplanted sample carrying cells of donor and recipient.
    """
    weaker = pl.min_horizontal("xist_log2cpm", "chry_log2cpm")
    return sexes.sort(weaker, descending=True)["sample"].head(n).to_list()


def sex_scatter(sexes: pl.DataFrame) -> ScatterData:
    """XIST against the chrY genes, coloured by the recorded sex, disagreements named."""
    unknown = text.word(PAGE, "figures.sex.unknown")
    frame = sexes.with_columns(
        recorded=pl.when(pl.col("expected") == "").then(pl.lit(unknown)).otherwise("expected"),
        label=pl.when(pl.col("agrees") == False).then("sample").otherwise(pl.lit("")),  # noqa: E712
    )  # fmt: skip
    top = float(frame.select(pl.max_horizontal(pl.col("xist_log2cpm").max(),
                                               pl.col("chry_log2cpm").max())).item())  # fmt: skip
    levels = [str(v) for v in frame["recorded"].unique(maintain_order=True)]
    return ScatterData(
        frame=frame, columns={
            "sample": ("none", text.word(PAGE, "columns.sample")),
            **{c: (u, text.word(PAGE, f"tables.sex_check.columns.{c}")) for c, u in
               (("expected", "none"), ("xist_log2cpm", "log2 CPM"), ("chry_log2cpm", "log2 CPM"),
                ("expression_sex", "none"), ("agrees", "none"))},
            "recorded": ("none", text.word(PAGE, "figures.sex.recorded")),
            "label": ("none", text.word(PAGE, "figures.sex.label")),
        },
        x="xist_log2cpm", y="chry_log2cpm", x_label=text.word(PAGE, "figures.sex.x"),
        y_label=text.word(PAGE, "figures.sex.y"), colour="recorded",
        colours=category_colours(levels), label="label",
        caption=text.fill(PAGE, "figures.sex.caption", n=frame.height, xist=SEX_X,
                          genes=", ".join(SEX_Y), cross=", ".join(cross_sex(sexes))),  # fmt: skip
        lines=(Line(text.word(PAGE, "figures.sex.boundary"), (0.0, top), (0.0, top)),),
        hover=("sample", "expected", "expression_sex"),
    )  # fmt: skip


def linked_pca(scores: pl.DataFrame, share: tuple[float, float], sheet: pl.DataFrame,
               status: pl.DataFrame, dimensions: Sequence[str]) -> LinkedData:  # fmt: skip
    """The PCA beside the samples' labels and QC status, for brushing."""
    worst = status.group_by("sample").agg(
        qc_status=pl.col("qc_status").sort_by(pl.col("qc_status").replace_strict(STATUS_ORDER))
        .last())  # fmt: skip
    dims = [d for d in dimensions if d in sheet.columns and sheet[d].n_unique() <= MAX_CATEGORIES]
    samples = sheet.unique("sample", keep="first", maintain_order=True).select("sample", *dims)
    frame = scores.join(samples, on="sample", maintain_order="left").join(
        worst, on="sample", maintain_order="left").with_columns(
        [pl.col(d).cast(pl.String).fill_null("") for d in dims])  # fmt: skip
    return LinkedData(
        frame=frame,
        columns={"sample": ("none", text.word(PAGE, "columns.sample")),
                 "pc1": ("VST units", text.word(PAGE, "figures.linked.pc1")),
                 "pc2": ("VST units", text.word(PAGE, "figures.linked.pc2")),
                 "qc_status": ("none", text.word(PAGE, "figures.linked.qc_status")),
                 **{d: ("none", text.fill(PAGE, "figures.fates.colour", column=d)) for d in dims}},
        x="pc1", y="pc2", x_label=text.fill(PAGE, "figures.linked.x", share=share[0]),
        y_label=text.fill(PAGE, "figures.linked.y", share=share[1]),
        dimensions=(*dims, "qc_status"),
        dimension_labels=(*dims, text.word(PAGE, "figures.linked.qc")),
        colour_dimension="qc_status",
        colours={s: BAND_COLOURS[b] for s, b in STATUS_BANDS.items()},
        caption=text.fill(PAGE, "figures.linked.caption", n=frame.height),
        names=tuple(text.words(PAGE, "figures.linked.names")),
        summary=text.raw(PAGE, "figures.linked.summary"),
        empty=text.word(PAGE, "figures.linked.empty"),
        detail=dims[0] if dims else "qc_status", hover=("sample", "qc_status"),
    )  # fmt: skip


# A sample's worst run decides its status on the linked view.
STATUS_ORDER = {"pass": 0, "flag": 1, "exclude": 2}
STATUS_BANDS = {"pass": "good", "flag": "bad", "exclude": "usable"}
