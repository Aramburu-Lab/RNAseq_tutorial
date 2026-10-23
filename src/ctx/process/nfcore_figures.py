"""M3's figure contracts on nf-core/rnaseq's alignment: every number nf-core's, the figures ctx's.

Bars, one per sample, stacked to 100 %: STAR's outcomes, RSeQC's strand shares, its read
distribution and junction kinds, and featureCounts' biotypes; one bar per sample of nf-core's
strand share beside its stranded threshold; curves, a line per sample: Qualimap's coverage
along transcripts, dupRadar's duplication against expression and Preseq's complexity. The
bars filter by the sheet's groups, and the strand bars also by nf-core's call (D-051); the
curves carry no control, as M2's FastQC curves carry none: each is a line per sample, and a
line filter would need a page control the curve contract lacks (D-116). Each describes and
tests nothing, so its caption says "No test"; every word is ``content/process_nfcore.toml``'s.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import polars as pl

from ctx.process.figures import OUTCOME_COLOURS
from ctx.process.nfcore import JUNCTION_KINDS, OUTCOMES
from ctx.qc.registry import Metric
from ctx.report import text
from ctx.viz.bar import BarData
from ctx.viz.curve import CurveData
from ctx.viz.scatter import Line
from ctx.viz.style import MAX_CATEGORIES, MAX_DISTINCT, category_colours

PAGE = "process_nfcore"
# RSeQC's three shares; the undetermined one in the reference grey, as no strand is its own.
STRANDS = ("rseqc_sense", "rseqc_antisense", "rseqc_undetermined")
# MultiQC's read_distribution categories by the part of a gene they lie in, the order the
# bars stack them: CDS, the two UTRs, introns, the 10 kb either side of a gene, the rest. A
# category none of these prefixes names is refused, never drawn as another.
REGIONS = {"CDS": "cds", "5'UTR": "utr5", "3'UTR": "utr3", "Intron": "introns",
           "TSS": "flanks", "TES": "flanks", "Other_intergenic": "intergenic"}  # fmt: skip
REGION_ORDER = ("cds", "utr5", "utr3", "introns", "flanks", "intergenic")
# Biotypes shown by name, the commonest over every sample; the rest summed as one, so a human
# annotation's dozens of biotypes keep within the page's colours.
BIOTYPES_SHOWN = MAX_CATEGORIES - 1
OTHER = "other"


def _columns(key: str, frame: pl.DataFrame,
             groups: Sequence[str] = ()) -> dict[str, tuple[str, str]]:  # fmt: skip
    """(unit, meaning) of every column, from ``figures.<key>``; each group column's own words."""
    units = text.raw(PAGE, f"figures.{key}.units")
    return {c: (units.get(c, "none"), text.fill(PAGE, "figures.group_column", column=c)
                if c in groups else text.word(PAGE, f"figures.{key}.columns.{c}"))
            for c in frame.columns}  # fmt: skip


def _labels(key: str, columns: Sequence[str]) -> dict[str, str]:
    """Words for columns on hover and on filters: ``figures.<key>.labels``, then the filters'."""
    words = {**text.raw(PAGE, "figures.filter_labels"), **text.raw(PAGE, f"figures.{key}.labels")}
    return {c: words.get(c, c) for c in columns}


@dataclass(frozen=True)
class Stack:
    """What one stacked bar figure stacks: its series column, levels, colours, hover, filters.

    ``levels`` name the series in legend order, each in ``figures.<key>.names`` unless
    ``names`` gives its words, as for biotypes, which the annotation names.
    """

    series: str
    levels: Sequence[str]
    colours: Mapping[str, str | None]
    hover: tuple[str, ...]
    filters: tuple[str, ...] = ()
    names: Mapping[str, str] | None = None


def _stacked(frame: pl.DataFrame, key: str, stack: Stack, groups: Sequence[str],
             caption: str) -> BarData:  # fmt: skip
    """One horizontal bar per sample, its series stacked to 100 %, filtered by the groups."""
    series, levels, colours = stack.series, stack.levels, stack.colours
    hover, filters = stack.hover, stack.filters
    return BarData(
        frame=frame, columns=_columns(key, frame, groups), category="sample", value="percent",
        category_label=text.word(PAGE, "figures.category"),
        value_label=text.word(PAGE, f"figures.{key}.value"), series=series,
        colours={lv: colours[lv] for lv in levels if lv in colours}, caption=caption,
        horizontal=True, stacked=True, hover=hover, width="double",
        names=dict(stack.names) if stack.names is not None
        else {lv: text.word(PAGE, f"figures.{key}.names.{lv}") for lv in levels},
        filters=(*filters, *groups), labels=_labels(key, [series, *hover, *filters, *groups]),
        summary_median=False,
    )  # fmt: skip


def mapping_bars(frame: pl.DataFrame, groups: Sequence[str], caption: str) -> BarData:
    """F1: STAR's outcome for every input read, STAR's own percent, one bar per sample.

    Args:
        frame: sample, outcome, reads, input_reads, percent, and the ``groups`` columns.
        groups: Sheet columns to filter by.
        caption: The filled caption.
    """
    stack = Stack("outcome", tuple(OUTCOMES), OUTCOME_COLOURS, ("reads", "input_reads"))
    return _stacked(frame, "mapping", stack, groups, caption)


def strand_bars(frame: pl.DataFrame, groups: Sequence[str], caption: str) -> BarData:
    """F2: RSeQC infer_experiment's sense, antisense and undetermined shares, per sample.

    Args:
        frame: sample, share, percent, call, provided, salmon_inferred, status, and the
            ``groups`` columns.
        groups: Sheet columns to filter by, beside nf-core's call.
        caption: The filled caption.
    """
    colours: dict[str, str | None] = {**category_colours(list(STRANDS[:2])), STRANDS[2]: None}
    stack = Stack("share", STRANDS, colours, ("call", "provided", "salmon_inferred", "status"),
                  ("call",))  # fmt: skip
    return _stacked(frame, "strand", stack, groups, caption)


def strand_check_bars(frame: pl.DataFrame, groups: Sequence[str], metric: Metric,
                      caption: str) -> BarData:  # fmt: skip
    """F2b: nf-core's rseqc_pct per sample, coloured by its call, beside its stranded threshold.

    The line is the registry's ``plot_line`` for ``rseqc_pct``, sourced to nf-core 3.27.0;
    the call is nf-core's own, read from its strand check table, never made here.
    """
    calls = frame["call"].unique(maintain_order=True).to_list()
    line = (
        (float(metric.plot_line), text.word(PAGE, "figures.strand_check.line"))
        if metric.plot_line is not None
        else None
    )
    return BarData(
        frame=frame, columns=_columns("strand_check", frame, groups), category="sample",
        value="rseqc_pct", category_label=text.word(PAGE, "figures.category"),
        value_label=text.word(PAGE, "figures.strand_check.value"), series="call",
        colours=category_colours(calls), caption=caption, horizontal=True, hline=line,
        hover=("status", "provided", "salmon_inferred"), width="double", filters=tuple(groups),
        labels=_labels("strand_check", ["call", "status", "provided", "salmon_inferred",
                                        *groups]),
    )  # fmt: skip


def distribution_frame(long: pl.DataFrame) -> pl.DataFrame:
    """RSeQC's tags per sample summed into ``REGIONS``, with each region's percent of the sample.

    Raises:
        ValueError: If a category matches no region: it would be drawn as another.
    """

    def region(category: str) -> str:
        """The region a MultiQC category lies in."""
        found = [r for prefix, r in REGIONS.items() if category.startswith(prefix)]
        if len(found) != 1:
            raise ValueError(f"RSeQC category {category!r} is in no region of {list(REGIONS)}")
        return found[0]

    mapped = long.with_columns(region=pl.col("category").map_elements(region, pl.String))
    summed = mapped.group_by("sample", "region", maintain_order=True).agg(pl.col("tags").sum())
    return summed.with_columns(percent=pl.col("tags") / pl.col("tags").sum().over("sample") * 100)


def distribution_bars(frame: pl.DataFrame, groups: Sequence[str], caption: str) -> BarData:
    """RSeQC read_distribution: each region's share of a sample's assigned tags."""
    levels = [r for r in REGION_ORDER if r in set(frame["region"])]
    stack = Stack("region", levels, category_colours(levels), ("tags",))
    return _stacked(frame, "distribution", stack, groups, caption)


def junction_bars(frame: pl.DataFrame, groups: Sequence[str], caption: str,
                  measure: str) -> BarData:  # fmt: skip
    """RSeQC junction_annotation: known, partly novel and novel, as RSeQC's percent, per sample.

    ``measure`` is ``events`` (spliced reads) or ``junctions`` (distinct junctions), one view
    each; ``frame`` holds both, and the bars draw one.
    """
    shown = frame.filter(pl.col("measure") == measure).drop("measure")
    stack = Stack("kind", JUNCTION_KINDS, category_colours(list(JUNCTION_KINDS)), ("count",))
    return _stacked(shown, "junctions", stack, groups, caption)


def biotype_frame(long: pl.DataFrame) -> pl.DataFrame:
    """The reads per biotype, the commonest named and the rest one ``other``, with each percent.

    The percent is of each sample's reads featureCounts counted on a biotype.
    """
    totals = long.group_by("biotype").agg(pl.col("reads").sum()).sort(
        ["reads", "biotype"], descending=[True, False])  # fmt: skip
    named = totals.filter(pl.col("reads") > 0)["biotype"].head(BIOTYPES_SHOWN).to_list()
    grouped = long.with_columns(pl.when(pl.col("biotype").is_in(named)).then("biotype")
                                .otherwise(pl.lit(OTHER)).alias("biotype"))  # fmt: skip
    summed = grouped.group_by("sample", "biotype", maintain_order=True).agg(pl.col("reads").sum())
    return summed.with_columns(percent=pl.col("reads") / pl.col("reads").sum().over("sample") * 100)


def biotype_bars(frame: pl.DataFrame, groups: Sequence[str], caption: str) -> BarData:
    """The biotypes featureCounts counted: each one's share of a sample's reads counted on one."""
    totals = frame.group_by("biotype").agg(pl.col("reads").sum()).sort(
        ["reads", "biotype"], descending=[True, False])  # fmt: skip
    named = [b for b in totals["biotype"] if b != OTHER]
    levels = named + [OTHER] * (OTHER in set(totals["biotype"]))
    names = {b: b for b in levels} | {OTHER: text.word(PAGE, "figures.biotypes.other")}
    colours: dict[str, str | None] = {**category_colours(named), OTHER: None}
    stack = Stack("biotype", levels, colours, ("reads",), names=names)
    return _stacked(frame, "biotypes", stack, groups, caption)


def curves(frame: pl.DataFrame, key: str, units: pl.DataFrame, colour: str | None,
           caption: str, thresholds: tuple[Line, ...] = ()) -> CurveData:  # fmt: skip
    """A line per sample of one of nf-core's curves; each its own colour and dash while they fit.

    Args:
        frame: sample, x, y, one row per point, in x order within a sample.
        key: The figure's words in content: ``coverage``, ``duplication`` or ``complexity``.
        units: sample and the colour column, one row per sample.
        colour: The sheet column colouring the lines when there are too many to tell apart.
        caption: The filled caption.
        thresholds: Dashed reference lines.
    """
    extra = [c for c in (colour,) if c and c in units.columns]
    joined = frame.join(units.select("sample", *extra).unique("sample", keep="first"),
                        on="sample", how="left", maintain_order="left")  # fmt: skip
    if not extra:
        colour = "group"
        joined = joined.with_columns(group=pl.lit(text.word(PAGE, "figures.all_samples")))
    assert colour is not None  # noqa: S101 - set just above when the sheet offers none
    joined = joined.with_columns(pl.col(colour).cast(pl.String))
    distinct = joined["sample"].n_unique() <= MAX_DISTINCT
    levels = [str(v) for v in joined[colour].unique(maintain_order=True)]
    return CurveData(
        frame=joined, columns=_columns(key, joined, (colour,)), series="sample", x="x", y="y",
        x_label=text.word(PAGE, f"figures.{key}.x"), y_label=text.word(PAGE, f"figures.{key}.y"),
        colour="sample" if distinct else colour, colours=category_colours(levels),
        distinct=distinct, hover=() if distinct else (colour,), thresholds=thresholds,
        width="double", caption=caption,
    )  # fmt: skip
