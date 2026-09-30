"""The QC page's figures, built from the run's tables; every word from ``content/qc.toml``."""

from collections.abc import Sequence

import numpy as np
import polars as pl
from scipy.stats import median_abs_deviation

from ctx.qc import Fastp
from ctx.qc.registry import Metric, Registry
from ctx.report import text
from ctx.viz.contour import Mark, SurfaceData
from ctx.viz.curve import CurveData
from ctx.viz.flow import FlowData, Node
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.strip import Span
from ctx.viz.style import BAND_COLOURS, MAX_CATEGORIES, ROLES, category_colours

PAGE = "qc"
# The composition shares a bubble's area adds up, where the counts had gene types.
SHARES = ("mt_fraction", "mt_rrna_fraction", "rrna_fraction")
# Trim Galore's fates and Salmon's, in the order the Sankey stacks them.
TRIM_FATES = ("reads_too_short", "reads_too_long", "reads_too_many_n", "reads_discarded_untrimmed")
SALMON_FATES = ("num_mapped", "num_decoy_fragments", "not_assigned")


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
    median, low, high = _band(frame["useful_share"].to_numpy(), registry.mad_multiplier)
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
    caption = text.fill(PAGE, "figures.fates.caption", n=frame.height, median=median,
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
    return CurveData(
        frame=frame, columns=columns, series="run_accession", x="position", y="quality",
        x_label=text.word(PAGE, "figures.quality.x"), y_label=text.word(PAGE, "figures.quality.y"),
        colour=colour, colours=category_colours(levels), hover=("sample",), spans=spans,
        caption=text.fill(PAGE, "figures.quality.caption", n=frame["run_accession"].n_unique(),
                          source=metric.plot_source or metric.source),
        width="double",
    )  # fmt: skip


def grid_figure(grid: pl.DataFrame, fastp: Fastp, run: str) -> SurfaceData | None:
    """Reads fastp keeps across -q and -l on one run; None when every cell keeps the same."""
    if grid["passed_pct"].n_unique() == 1:
        return None
    units = text.raw(PAGE, "figures.grid.units")
    return SurfaceData(
        frame=grid, x="q", y="l", z="passed_pct",
        columns={c: (units.get(c, "none"), text.word(PAGE, f"figures.grid.columns.{c}"))
                 for c in grid.columns},  # fmt: skip
        x_label=text.word(PAGE, "figures.grid.x"), y_label=text.word(PAGE, "figures.grid.y"),
        z_label=text.word(PAGE, "figures.grid.z"),
        caption=text.fill(PAGE, "figures.grid.caption", run=run, u=fastp.u,
                          nq=grid["q"].n_unique(), nl=grid["l"].n_unique(), isoline=GRID_ISOLINE),
        isolines=(GRID_ISOLINE,), isoline_format=".0f",
        mark=Mark(float(fastp.q), float(fastp.l), text.fill(PAGE, "figures.grid.mark", q=fastp.q,
                                                              l=fastp.l)),
    )  # fmt: skip
