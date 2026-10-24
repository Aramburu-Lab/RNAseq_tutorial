"""What the M3 page on nf-core's alignment shows: each figure block and each small table.

Every block takes the run's frames and what the page states beside them (``Found``,
``Shown``), so ``ctx figure`` redraws a figure by calling the block that drew it on the frames
read back from the tables (``ctx.figures.draw_m3``). Every word is
``content/process_nfcore.toml``'s; every number is nf-core's, or arithmetic on nf-core's.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import polars as pl

from ctx.de.run import Writer
from ctx.process import nfcore
from ctx.process import nfcore_figures as fig
from ctx.process import tables as tab
from ctx.qc.blocks import exported
from ctx.qc.registry import Metric
from ctx.report import text
from ctx.report.builder import FigureBlock, FigureText, ResultTable
from ctx.viz.scatter import Line
from ctx.viz.style import MAX_DISTINCT

PAGE = "process_nfcore"
# The tables of what M3 measured, which a lite page keeps; the rest record its inputs (D-111).
RESULTS_TABLES = ("alignment", "strand_check", "concordance")
EMPTY = pl.DataFrame()


@dataclass(frozen=True)
class Found:
    """The run's frames, each one sample per row or one point per row (``nfcore_page.frames``).

    A redraw fills only what its figure reads; the rest stay empty.
    """

    alignment: pl.DataFrame = field(default_factory=pl.DataFrame)
    outcomes: Mapping[str, Mapping[str, int]] = field(default_factory=dict)
    strand: pl.DataFrame = field(default_factory=pl.DataFrame)
    distribution: pl.DataFrame = field(default_factory=pl.DataFrame)
    junctions: pl.DataFrame = field(default_factory=pl.DataFrame)
    coverage: pl.DataFrame = field(default_factory=pl.DataFrame)
    duplication: pl.DataFrame = field(default_factory=pl.DataFrame)
    complexity: pl.DataFrame = field(default_factory=pl.DataFrame)
    perfect: pl.DataFrame = field(default_factory=pl.DataFrame)
    biotypes: pl.DataFrame = field(default_factory=pl.DataFrame)
    comparison: pl.DataFrame = field(default_factory=pl.DataFrame)
    agreements: tab.Agreements = field(default_factory=lambda: tab.Agreements(EMPTY, EMPTY, EMPTY))
    # Per matrix, the genes only it holds and their share of its counts (nfcore.matrices).
    unshared: Mapping[str, tuple[int, float]] = field(default_factory=dict)
    steps: pl.DataFrame = field(default_factory=pl.DataFrame)


@dataclass(frozen=True)
class Shown:
    """What the page states beside the frames.

    Attributes:
        groups: sample, then each sheet column the bars filter by.
        colour: The sheet column that colours the curves past ``MAX_DISTINCT`` samples.
        strand_metric: The registry's ``rseqc_pct``, whose plot line is nf-core's threshold.
        prefix: Put before every caption of the counter figures, which ``process.blocks``
            draws: empty here.
    """

    groups: pl.DataFrame
    colour: str | None
    strand_metric: Metric
    prefix: str = ""

    @property
    def usable(self) -> list[str]:
        """The group columns the bars filter by."""
        return self.groups.columns[1:]


def _text(key: str, lead: str, **values: Any) -> FigureText:
    """Title, lead and reading guide of ``figures.<key>``, each line filled from ``values``."""
    reading = tuple(line.format_map(values) for line in text.words(PAGE, f"figures.{key}.reading"))
    return FigureText(title=text.word(PAGE, f"figures.{key}.title"), lead=lead, reading=reading,
                      tab=text.word(PAGE, f"figures.{key}.tab"))  # fmt: skip


def _median(values: pl.Series) -> float:
    """A column's median as a float."""
    found: Any = values.median()
    return float(found)


def _joined(frame: pl.DataFrame, shown: Shown) -> pl.DataFrame:
    """A per-sample frame with the sheet's group columns, for the bars' filters."""
    return frame.join(shown.groups, on="sample", how="left", maintain_order="left")


def mapping_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """F1: STAR's outcome for every input read of every sample."""
    rows = [(s, k, r[k], r["input"]) for s, r in found.outcomes.items() for k in nfcore.OUTCOMES]
    frame = pl.DataFrame(rows, schema={"sample": pl.String, "outcome": pl.String,
                                       "reads": pl.Int64, "input_reads": pl.Int64},
                         orient="row")  # fmt: skip
    pct = found.alignment.unpivot(index="sample", on=[f"{k}_pct" for k in nfcore.OUTCOMES],
                                  variable_name="outcome", value_name="percent")  # fmt: skip
    frame = frame.join(pct.with_columns(pl.col("outcome").str.strip_suffix("_pct")),
                       on=["sample", "outcome"], maintain_order="left")  # fmt: skip
    unique = found.alignment["unique_pct"]
    low = found.alignment.sort("unique_pct").row(0, named=True)
    values = {"n": found.alignment.height, "median": f"{_median(unique):.1f}",
              "low": f"{low['unique_pct']:.1f}", "low_sample": low["sample"],
              "reads": f"{_median(found.alignment['input_reads']):,.0f}"}  # fmt: skip
    caption = text.fill(PAGE, "figures.mapping.caption", **values)
    data = fig.mapping_bars(_joined(frame, shown), shown.usable, caption)
    return writer.figure(
        data, "aligned_mapping", _text("mapping", text.fill(PAGE, "figures.mapping.lead", **values))
    )


def strand_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """F2: RSeQC's sense, antisense and undetermined shares per sample, with nf-core's call."""
    strand = found.strand.rename({"rseqc_inferred": "call"})
    frame = strand.unpivot(
        index=["sample", "call", "provided", "salmon_inferred", "status"],
        on=list(fig.STRANDS),
        variable_name="share",
        value_name="percent",
    )
    counts = strand["call"].value_counts(sort=True)
    call, n_call = counts.row(0)
    values = {"n": strand.height, "call": call, "n_call": n_call,
              "failed": strand.filter(pl.col("status") != "pass").height}  # fmt: skip
    caption = text.fill(PAGE, "figures.strand.caption", **values)
    data = fig.strand_bars(_joined(frame, shown), shown.usable, caption)
    return writer.figure(
        data, "rseqc_strand", _text("strand", text.fill(PAGE, "figures.strand.lead", **values))
    )


def strand_check_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """F2b: nf-core's share for its call per sample, beside its stranded threshold."""
    frame = found.strand.select("sample", "rseqc_pct", "provided", "salmon_inferred", "status",
                                call="rseqc_inferred")  # fmt: skip
    metric = shown.strand_metric
    low: Any = frame["rseqc_pct"].min()
    values = {"n": frame.height, "low": f"{float(low):.1f}",
              "line": f"{metric.plot_line:g}", "source": metric.plot_source}  # fmt: skip
    caption = text.fill(PAGE, "figures.strand_check.caption", **values)
    data = fig.strand_check_bars(_joined(frame, shown), shown.usable, metric, caption)
    return writer.figure(
        data,
        "rseqc_share",
        _text("strand_check", text.fill(PAGE, "figures.strand_check.lead", **values), **values),
    )


def distribution_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """RSeQC read_distribution: where each sample's assigned tags lie in or near genes."""
    frame = fig.distribution_frame(found.distribution)
    cds = frame.filter(pl.col("region") == "cds")["percent"]
    values = {"n": frame["sample"].n_unique(), "cds": f"{_median(cds):.1f}" if cds.len() else
              text.word(PAGE, "undefined")}  # fmt: skip
    caption = text.fill(PAGE, "figures.distribution.caption", **values)
    data = fig.distribution_bars(_joined(frame, shown), shown.usable, caption)
    return writer.figure(
        data,
        "read_distribution",
        _text("distribution", text.fill(PAGE, "figures.distribution.lead", **values)),
    )


def junction_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """RSeQC junction_annotation: known, partly novel and novel, as events and as junctions.

    One view per measure (D-051): the spliced reads, or the distinct junctions they span.
    """
    frame = _joined(found.junctions, shown)
    options = []
    for measure in nfcore.JUNCTION_MEASURES:
        part = frame.filter(pl.col("measure") == measure)
        known = part.filter(pl.col("kind") == "known")["percent"]
        said = text.word(PAGE, f"figures.junctions.measures.{measure}")
        values = {"n": part["sample"].n_unique(), "known": f"{_median(known):.1f}",
                  "measure": said}  # fmt: skip
        caption = text.fill(PAGE, "figures.junctions.caption", **values)
        data = fig.junction_bars(frame, shown.usable, caption, measure)
        options.append(
            (
                said,
                data,
                _text("junctions", text.fill(PAGE, "figures.junctions.lead", **values)),
            )
        )
    block = writer.figure(options[0][1], "junction_annotation", options[0][2])
    return exported(block, options, text.word(PAGE, "figures.junctions.view"))


def _units(shown: Shown) -> pl.DataFrame:
    """The sample and the curves' colour column, one row per sample."""
    return shown.groups.select("sample", *([shown.colour] if shown.colour in shown.groups.columns
                                           else []))  # fmt: skip


def _styles(n: int) -> str:
    """The reading-guide line on how the curves are told apart."""
    return (text.word(PAGE, "figures.styles_distinct") if n <= MAX_DISTINCT
            else text.fill(PAGE, "figures.styles_shared", most=MAX_DISTINCT))  # fmt: skip


def coverage_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """Qualimap: each sample's coverage along its transcripts, 5' to 3'."""
    n = found.coverage["sample"].n_unique()
    caption = text.fill(PAGE, "figures.coverage.caption", n=n)
    data = fig.curves(found.coverage, "coverage", _units(shown), shown.colour, caption)
    words = _text("coverage", text.fill(PAGE, "figures.coverage.lead", n=n))
    return writer.figure(
        data, "gene_coverage", replace(words, reading=(*words.reading, _styles(n)))
    )


# dupRadar's fit holds 2,230 to 3,594 points a sample (Breast-tissue-BRCA, 2026-10-24): drawn
# whole, 103 samples made a 16.6 MB page, over the 15 MB cap. Drawn at about this many; the
# duplication_curve table keeps every point.
CURVE_POINTS = 200


def thinned(frame: pl.DataFrame, points: int = CURVE_POINTS) -> pl.DataFrame:
    """Every k-th point of each sample's curve, in its order, its first and last always kept.

    Args:
        frame: sample, x, y, each sample's points in x order.
        points: About how many points each sample keeps; never more than one over.
    """
    at, n = pl.int_range(pl.len()).over("sample"), pl.len().over("sample")
    step = ((n + points - 1) // points).clip(lower_bound=1)
    return frame.filter((at % step == 0) | (at == n - 1))


def duplication_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """dupRadar: each sample's duplication against expression, log10 reads per kilobase."""
    frame = thinned(found.duplication).with_columns(pl.col("x").log10())
    n = frame["sample"].n_unique()
    caption = text.fill(PAGE, "figures.duplication.caption", n=n, points=CURVE_POINTS)
    data = fig.curves(frame, "duplication", _units(shown), shown.colour, caption)
    words = _text("duplication", text.fill(PAGE, "figures.duplication.lead", n=n))
    return writer.figure(data, "duplication", replace(words, reading=(*words.reading, _styles(n))))


def complexity_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """Preseq: each sample's expected distinct molecules against molecules sequenced."""
    n = found.complexity["sample"].n_unique()
    perfect = found.perfect.sort("x")
    lines = (Line(text.word(PAGE, "figures.complexity.perfect"),
                  tuple(perfect["x"].to_list()), tuple(perfect["y"].to_list())),) \
        if perfect.height > 1 else ()  # fmt: skip
    caption = text.fill(PAGE, "figures.complexity.caption", n=n)
    data = fig.curves(found.complexity, "complexity", _units(shown), shown.colour, caption,
                      thresholds=lines)  # fmt: skip
    words = _text("complexity", text.fill(PAGE, "figures.complexity.lead", n=n))
    return writer.figure(
        data,
        "library_complexity",
        replace(words, reading=(*words.reading, _styles(n))),
    )


def biotype_figure(writer: Writer, found: Found, shown: Shown) -> FigureBlock:
    """The biotypes featureCounts counted: what each sample's reads on genes are."""
    frame = fig.biotype_frame(found.biotypes)
    top = frame.group_by("biotype").agg(pl.col("reads").sum()).sort(
        ["reads", "biotype"], descending=[True, False]).row(0, named=True)  # fmt: skip
    share = frame.filter(pl.col("biotype") == top["biotype"])["percent"]
    values = {"n": frame["sample"].n_unique(), "top": top["biotype"],
              "share": f"{_median(share):.1f}", "named": fig.BIOTYPES_SHOWN}  # fmt: skip
    caption = text.fill(PAGE, "figures.biotypes.caption", **values)
    data = fig.biotype_bars(_joined(frame, shown), shown.usable, caption)
    return writer.figure(
        data,
        "biotypes",
        _text("biotypes", text.fill(PAGE, "figures.biotypes.lead", **values), **values),
    )


def page_table(key: str, rows: Sequence[Sequence[str]], numeric: tuple[int, ...],
               **values: Any) -> ResultTable:  # fmt: skip
    """A small table on the page, its title and headings from ``tables.<key>.page``."""
    return ResultTable(text.fill(PAGE, f"tables.{key}.page.title", **values),
                       tuple(text.words(PAGE, f"tables.{key}.page.columns")),
                       tuple(tuple(r) for r in rows), numeric=numeric,
                       results=key in RESULTS_TABLES)  # fmt: skip


def _number(value: float | None, spec: str) -> str:
    """A number as the page prints it, or the content's word for one undefined."""
    return text.word(PAGE, "undefined") if value is None else format(value, spec)


def alignment_table(found: Found) -> ResultTable:
    """Per sample: Trim Galore's reads, STAR's input and unique share, splices, records."""
    rows = [(r["sample"], f"{r['trim_out']:,}", f"{r['input_reads']:,}",
             _number(r["unique_pct"], ".1f"), _number(r["splices_per_100_unique"], ".1f"),
             _number(r["splices_annotated_share"], ".3f"), f"{r['records']:,}",
             _number(r["primary_mapped_pct"], ".1f"))
            for r in found.alignment.iter_rows(named=True)]  # fmt: skip
    return page_table("alignment", rows, tuple(range(1, 8)))


def strand_table(found: Found) -> ResultTable:
    """Per sample: RSeQC's three shares, nf-core's share and call, Salmon's, and the check."""
    rows = [(r["sample"], _number(r["rseqc_sense"], ".1f"), _number(r["rseqc_antisense"], ".1f"),
             _number(r["rseqc_undetermined"], ".1f"), r["provided"], r["rseqc_inferred"],
             r["salmon_inferred"], r["status"])
            for r in found.strand.iter_rows(named=True)]  # fmt: skip
    return page_table("strand_check", rows, (1, 2, 3))


def concordance_table(found: Found) -> ResultTable:
    """The comparison over every sample's gene sums."""
    agree = found.agreements.concordance.filter(pl.col("unit") == "all")
    rows = [(text.word("process", f"comparisons.{r['comparison']}.name"), f"{r['n_genes']:,}",
             f"{r['n_identical']:,}", _number(r["spearman"], ".3f"),
             f"{_number(r['spearman_low'], '.3f')} to {_number(r['spearman_high'], '.3f')}",
             _number(r["bias"], ".3f"),
             f"{_number(r['loa_low'], '.3f')} to {_number(r['loa_high'], '.3f')}")
            for r in agree.iter_rows(named=True)]  # fmt: skip
    return page_table("concordance", rows, (1, 2, 3, 5))


def steps_table(found: Found) -> ResultTable:
    """What each of nf-core's processes behind the page cost, from its execution trace."""
    rows = [(r["process"], f"{r['tasks']:,}", f"{r['samples']:,}",
             _number(r["peak_rss"] / 2**20 if r["peak_rss"] is not None else None, ",.0f"),
             _number(r["realtime"] / 60 if r["realtime"] is not None else None, ",.1f"))
            for r in found.steps.iter_rows(named=True)]  # fmt: skip
    return page_table("steps", rows, (1, 2, 3, 4))


def reference_table(params: Mapping[str, object]) -> ResultTable:
    """The genome and annotation nf-core aligned to, as the run's params file names them."""
    rows = [(text.word(PAGE, f"tables.reference.rows.{k}"), str(params.get(k)))
            for k in ("fasta", "gtf", "gencode", "aligner", "pseudo_aligner")]  # fmt: skip
    return page_table("reference", rows, ())
