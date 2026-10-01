"""The QC page's figures, built from the run's tables; every word from ``content/qc.toml``."""

from collections.abc import Sequence
from dataclasses import replace

import numpy as np
import polars as pl
from scipy.stats import median_abs_deviation

from ctx.de.run import SEX_X, SEX_Y
from ctx.qc import Fastp
from ctx.qc.registry import Metric, Registry
from ctx.report import text
from ctx.viz.contour import DENSITY_MIN_POINTS, HDR_MASSES, Mark, SurfaceData
from ctx.viz.curve import CurveData
from ctx.viz.density import DensityData
from ctx.viz.flow import FlowData, Node
from ctx.viz.heatmap import HeatmapData, Swatch
from ctx.viz.linked import LinkedData
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.strip import Range, Span, StripData
from ctx.viz.style import BAND_COLOURS, MAX_CATEGORIES, MAX_DISTINCT, ROLES, category_colours

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
    frame = fates.join(samples.select("sample", colour), on="sample", how="left").join(
        values.select("sample", *shares), on="sample", how="left")  # fmt: skip
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


def quality_figure(curves: pl.DataFrame, sheet: pl.DataFrame, metric: Metric,
                   colour: str | None) -> CurveData:  # fmt: skip
    """Mean quality along read 1 of every run, over the bands FastQC's own plot shades."""
    runs = sheet.select("run_accession", "sample", *([colour] if colour else []))
    frame = curves.join(runs, on="run_accession", how="left")
    if colour is None:
        colour = "group"
        frame = frame.with_columns(group=pl.lit(text.word(PAGE, "figures.fates.all")))
    frame = frame.with_columns(pl.col(colour).cast(pl.String))
    good, usable = metric.plot_bands or (None, None)
    spans = () if good is None else (
        Span("good", text.word(PAGE, "figures.quality.good"), good, None),
        Span("usable", text.word(PAGE, "figures.quality.usable"), usable, good),
        Span("bad", text.word(PAGE, "figures.quality.bad"), None, usable),
    )  # fmt: skip
    units = text.raw(PAGE, "figures.quality.units")
    columns = {c: (units.get(c, "none"), text.fill(PAGE, "figures.fates.colour", column=c)
                   if c == colour else text.word(PAGE, f"figures.quality.columns.{c}"))
               for c in frame.columns}  # fmt: skip
    levels = [str(v) for v in frame[colour].unique(maintain_order=True)]
    # Each run its own colour and dash while they fit (Rytis 2026-09-30); the group moves to
    # the hover. Beyond that, lines share their group's colour.
    distinct = frame["run_accession"].n_unique() <= MAX_DISTINCT
    hover = ("sample", colour) if distinct and colour != "group" else ("sample",)
    return CurveData(
        frame=frame, columns=columns, series="run_accession", x="position", y="quality",
        x_label=text.word(PAGE, "figures.quality.x"), y_label=text.word(PAGE, "figures.quality.y"),
        colour="run_accession" if distinct else colour, colours=category_colours(levels),
        distinct=distinct, hover=hover, spans=spans,
        caption=text.fill(PAGE, "figures.quality.caption", n=frame["run_accession"].n_unique(),
                          source=metric.plot_source or metric.source),
        width="double",
    )  # fmt: skip


def grid_figure(grid: pl.DataFrame, fastp: Fastp, run: str) -> SurfaceData | None:
    """Reads fastp keeps across -q and -l on one run; None when every cell keeps the same."""
    if grid["passed_pct"].n_unique() == 1:
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
    )  # fmt: skip
    # The band count is the one the surface draws, so the caption cannot name another.
    return replace(surface, caption=text.fill(
        PAGE, "figures.grid.caption", run=run, u=fastp.u, nq=grid["q"].n_unique(),
        nl=grid["l"].n_unique(), isoline=GRID_ISOLINE, bands=len(surface.edges()) - 1))  # fmt: skip


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
    frame = values.select(key, metric.column).drop_nulls().join(units, on=key, how="left")
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
    word = {b: text.word(PAGE, f"figures.strip.{b}") for b in ("good", "usable", "bad")}
    spans = tuple(Span(b, word[b], low, high) for b, low, high in metric.spans())
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
                          bands=text.fill(PAGE, "figures.strip.spans", source=metric.source)
                          if spans else text.word(PAGE, "figures.strip.no_spans")),
        flagged="flagged", flag_label=text.word(PAGE, "figures.strip.flagged_label"),
        hover=(key, "reason"), spans=spans, ranges=ranges,
        range_label=text.word(PAGE, "figures.strip.range"), log=log,
    )  # fmt: skip


def band_map(values: pl.DataFrame, registry: Registry, key: str) -> HeatmapData | None:
    """Every unit's absolute band on every metric that has one; None where none has."""
    metrics = [m for m in registry.metrics.values()
               if m.column in values.columns and (m.banded or m.kind == "status")]  # fmt: skip
    if not metrics:
        return None
    none = text.word(PAGE, "figures.band_map.none")
    rows = [(unit, m.title, m.band(row[m.column]) or none)
            for row in values.iter_rows(named=True) for m in metrics
            for unit in (row[key],)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=[key, "metric", "band"], orient="row").with_columns(
        value=pl.lit(None, pl.Float64))  # fmt: skip
    return HeatmapData(
        frame=frame, row="metric", column=key, value="value", level="band",
        columns={key: ("none", text.word(PAGE, f"columns.{key}")),
                 "metric": ("none", text.word(PAGE, "figures.band_map.metric")),
                 "band": ("none", text.word(PAGE, "figures.band_map.band")),
                 "value": ("none", text.word(PAGE, "figures.band_map.value"))},
        row_label=text.word(PAGE, "figures.band_map.metric"),
        column_label=text.word(PAGE, f"columns.{key}"),
        swatches=(*(Swatch(b, BAND_COLOURS[b]) for b in ("good", "usable", "bad")),
                  Swatch(none, None)),
        caption=text.fill(PAGE, "figures.band_map.caption", n=values.height, m=len(metrics)),
        width="double",
    )  # fmt: skip


def density(values: pl.DataFrame, sheet: pl.DataFrame, flags: pl.DataFrame,
            colour: str | None) -> DensityData | None:  # fmt: skip
    """Library size against genes detected, a density per group of 30 or more samples."""
    if not {"library_size", "genes_detected"} <= set(values.columns):
        return None
    samples = sheet.unique("sample", keep="first").select("sample", *([colour] if colour else []))
    frame = values.select("sample", "library_size", "genes_detected").join(samples, on="sample")
    frame = frame.filter(pl.col("library_size") > 0)
    group = colour or "group"
    if colour is None:
        frame = frame.with_columns(group=pl.lit(text.word(PAGE, "figures.fates.all")))
    flagged = pl.col("sample").is_in(flags["sample"].implode())
    frame = frame.with_columns(pl.col(group).cast(pl.String), flagged=flagged)
    levels = [str(v) for v in frame[group].unique(maintain_order=True)]
    data = DensityData(
        frame=frame, x="library_size", y="genes_detected", group=group,
        x_label=text.word(PAGE, "figures.density.x"), y_label=text.word(PAGE, "figures.density.y"),
        columns={"sample": ("none", text.word(PAGE, "columns.sample")),
                 "library_size": ("reads", text.word(PAGE, "figures.density.library_size")),
                 "genes_detected": ("genes", text.word(PAGE, "figures.density.genes_detected")),
                 group: ("none", text.fill(PAGE, "figures.fates.colour", column=group)),
                 "flagged": ("none", text.word(PAGE, "figures.density.flagged"))},
        colours=category_colours(levels), caption="", highlight="flagged",
        highlight_label=text.word(PAGE, "figures.strip.flagged_label"),
        hover=("sample", "flagged"), log_x=True, width="double",
    )  # fmt: skip
    n = ", ".join(f"{g} {d.n}" for g, d in data.densities.items()) or text.word(
        PAGE, "figures.density.none")  # fmt: skip
    return replace(data, caption=text.fill(PAGE, "figures.density.caption", n=n,
                                           bandwidth=data.bandwidth(),
                                           minimum=DENSITY_MIN_POINTS, **HDR))  # fmt: skip


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
                          genes=", ".join(SEX_Y)),  # fmt: skip
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
    samples = sheet.unique("sample", keep="first").select("sample", *dims)
    frame = scores.join(samples, on="sample").join(worst, on="sample").with_columns(
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
