"""What the M3 page shows: each figure block, each small table and the notes beside them.

Every word is ``content/process.toml``'s; every number comes from ``page.Frames``. A caption
starts with ``Shown.prefix``, empty on the real page and "Made-up data." on the preview, so a
preview figure says so in every file it writes.
"""

from collections.abc import Sequence
from dataclasses import replace
from typing import TYPE_CHECKING, Any

import polars as pl

from ctx import concordance
from ctx.de.run import Writer
from ctx.process import Aligner, figures, views
from ctx.process import strand as strand_share
from ctx.process import tables as tab
from ctx.report import text
from ctx.report.builder import FigureBlock, FigureText, ResultTable
from ctx.report.lookups import lookup
from ctx.viz.aggregates import Aggregate, Lead

if TYPE_CHECKING:
    from ctx.process.page import Frames, Shown

PAGE = "process"
# How the folder table words a keyed folder: built by this run, or reused.
BUILT = {True: "built", False: "reused"}
# More contigs than this are counted, not listed, in the text and the reference table.
SHOWN_CONTIGS = 3


def _text(key: str, lead: str, **values: Any) -> FigureText:
    """Title, lead and reading guide of ``figures.<key>``, the title filled from ``values``."""
    return FigureText(title=text.fill(PAGE, f"figures.{key}.title", **values), lead=lead,
                      reading=tuple(line.format_map(values)
                                    for line in text.words(PAGE, f"figures.{key}.reading")),
                      tab=text.fill(PAGE, f"figures.{key}.tab", **values) if values
                      else "")  # fmt: skip


def bands(shown: "Shown") -> str:
    """Whether genome-wide mapping bands apply: not on a partial reference."""
    return text.word(PAGE, f"figures.mapping.bands.{'partial' if shown.partial else 'whole'}")


def mapping_lead() -> Lead:
    """F1's lead over the samples its filter shows: the median uniquely mapped and the lowest.

    Each from a bar's uniquely mapped segment, as a share of its input reads.
    """
    unique = {"outcome": "unique"}
    return Lead(text.raw(PAGE, "figures.mapping.lead"),
                (Aggregate("median", "median", "percent", unique),
                 Aggregate("low", "min", "percent", unique, row=("low_sample", "label"))),
                note=text.word(PAGE, "figures.mapping.note"))  # fmt: skip


def mapping_figure(writer: Writer, found: "Frames", shown: "Shown") -> FigureBlock:
    """F1: every sample as run, one stacked bar each, filtered by the sheet's groups."""
    frame = tab.outcomes_long(found.outcomes, shown.aligner).join(
        shown.groups, left_on="label", right_on="sample")  # fmt: skip
    caption = shown.prefix + text.fill(
        PAGE, "figures.mapping.caption", aligner=str(shown.aligner).upper(),
        n=len(found.outcomes), contigs=shown.genome_name, bands=bands(shown),
        reads=figures.reads(figures.stat(found.alignment["input_reads"], "median")),
    )  # fmt: skip
    data = figures.outcome_bars(frame, shown.aligner, "mapping", shown.usable, caption)
    return writer.figure(data, "mapping", replace(_text("mapping", ""), follows=mapping_lead()))


def settings_figure(writer: Writer, found: "Frames", shown: "Shown") -> FigureBlock:
    """F1b: the featured sample as run and under each setting, its few bars always labelled."""
    named = {text.word(PAGE, f"settings.{k}"): v for k, v in found.settings.items()}
    caption = shown.prefix + text.fill(PAGE, "figures.settings.caption",
                                       aligner=str(shown.aligner).upper(),
                                       sample=shown.featured, n=len(named) - 1)  # fmt: skip
    frame = tab.outcomes_long(named, shown.aligner)
    data = figures.outcome_bars(frame, shown.aligner, "settings", (), caption)
    unique = {k: 100 * r["unique"] / r["input"] if r["input"] else 0.0
              for k, r in found.settings.items()}  # fmt: skip
    # Under STAR the derived limit 1 and under HISAT2 the unspliced alignment are always there.
    low = min((k for k in unique if k != "main"), key=unique.__getitem__)
    lead = text.fill(PAGE, "figures.settings.lead", main=figures.number(unique["main"], ".1f"),
                     low=figures.number(unique[low], ".1f"),
                     setting=text.word(PAGE, f"settings.{low}"))  # fmt: skip
    return writer.figure(data, "featured_settings", _text("settings", lead))


def strand_lead() -> Lead:
    """F2's lead over the samples shown: their forward shares' range and the commonest call.

    Each sample is read once, from its -s 0 bar.
    """
    once = {"setting": "s0"}
    return Lead(text.raw(PAGE, "figures.strandedness.lead"),
                (Aggregate("low", "min", "forward_share", once),
                 Aggregate("high", "max", "forward_share", once),
                 Aggregate("call", "mode", "call", once, count="n_call"),
                 Aggregate("n", "count", where=once)),
                note=text.word(PAGE, "figures.strandedness.note"))  # fmt: skip


def strand_figure(writer: Writer, found: "Frames", shown: "Shown") -> FigureBlock:
    """F2: what featureCounts assigns at each -s, per sample, with the forward share and call.

    Each bar carries its sample's forward share whole: the lead takes its range to three
    significant digits, and a share held to three decimals printed 0.021 as 0.0210.
    """
    strand, inputs = found.strand, dict(found.alignment.select("sample", "input_reads").rows())
    rows = [{"sample": r["sample"], "setting": s, "reads": int(r[f"assigned_{s}"]),
             "percent": 100 * r[f"assigned_{s}"] / inputs[r["sample"]],
             "no_feature": int(r[f"no_feature_{s}"]), "ambiguous": int(r[f"ambiguous_{s}"]),
             "forward_share": r["forward_share"], "call": r["call"]}
            for r in strand.iter_rows(named=True) for s in ("s0", "s1", "s2")]  # fmt: skip
    frame = pl.DataFrame(rows).join(shown.groups, on="sample")
    metric = shown.strand_metric
    caption = shown.prefix + text.fill(
        PAGE, "figures.strandedness.caption", n=strand.height,
        share=f"{figures.stat(strand['forward_share'], 'median'):.3f}",
        stranded=metric.stranded, unstranded=metric.unstranded)  # fmt: skip
    data = figures.strand_bars(frame, shown.usable, caption)
    words = replace(_text("strandedness", ""), follows=strand_lead())
    if shown.aligner is Aligner.STAR:
        star = text.word(PAGE, "figures.strandedness.reading_star")
        words = replace(words, reading=(*words.reading, star))
    return writer.figure(data, "strand_assigned", words)


def share_figure(writer: Writer, found: "Frames", shown: "Shown") -> FigureBlock:
    """F2b: each sample's forward share and call, its thresholds moved live (WO8)."""
    metric = shown.strand_metric
    shares = found.strand.filter(pl.col("forward_share").is_not_null())["forward_share"]
    caption = shown.prefix + text.fill(
        PAGE, "figures.share.caption", n=shares.len(), stranded=metric.stranded,
        unstranded=metric.unstranded)  # fmt: skip
    data = strand_share.share_plot(found.strand, metric, caption)
    counts = strand_share.counted(data.frame["call"].to_list())
    lead = text.fill(PAGE, "figures.share.lead", n=data.frame.height, **counts)
    return writer.figure(data, "strand_share", _text("share", lead))


def agreement_figures(writer: Writer, found: "Frames", shown: "Shown"
                      ) -> tuple[list[FigureBlock], list[FigureBlock]]:  # fmt: skip
    """F3 and F4 for every comparison the aligner allows, in ``tables.COMPARISONS`` order."""
    table, sums = found.agreements.concordance, found.agreements.sums
    long = views.typed(found.comparison, sums)
    scatters, plots = [], []
    for name in table["comparison"].unique(maintain_order=True):
        rows = table.filter(pl.col("comparison") == name)
        every = rows.filter(pl.col("unit") == "all").row(0, named=True)
        words = figures.comparison_words(name)
        values = figures.scatter_values(every, rows.filter(pl.col("unit") != "all"),
                                        len(found.outcomes)) | words  # fmt: skip
        part = sums.filter(pl.col("comparison") == name)
        scatter = figures.counter_scatter(part, words, shown.prefix + text.fill(
            PAGE, "figures.counters.caption", **values))  # fmt: skip
        built = _text("counters", text.fill(PAGE, "figures.counters.lead", **values), **words)
        block = writer.figure(scatter, f"counters_{name}", built)
        scatters.append(_scatter_views(block, (scatter, built), long, name, words, shown))
        plot = figures.agreement_plot(part, every, words, shown.prefix + text.fill(
            PAGE, "figures.agreement.caption", **values))  # fmt: skip
        outside = text.word(PAGE, "undefined")
        if every["loa_low"] is not None:
            beyond = ~pl.col("diff_log2").is_between(every["loa_low"], every["loa_high"])
            outside = f"{plot.frame.filter(beyond).height:,}"
        built = _text("agreement", text.fill(PAGE, "figures.agreement.lead", **values,
                                             n_out=outside), **words)  # fmt: skip
        block = writer.figure(plot, f"agreement_{name}", built)
        plots.append(_plot_views(block, (plot, built), long, name, words, shown))
    return scatters, plots


def _scatter_views(block: FigureBlock, build: tuple[Any, FigureText], long: pl.DataFrame,
                   name: str, words: dict[str, str], shown: "Shown") -> FigureBlock:  # fmt: skip
    """F3 with a gene type and a sample as views (WO8): the build's every gene, every sample."""
    part = long.filter(pl.col("comparison") == name)
    found, axes, at = views.subsets(part)
    options = [(found[0].words(), *build)]
    for subset in found[1:]:
        data, filled = views.scatter(part, subset, words, shown.prefix)
        options.append((subset.words(), data, _text(
            "counters", text.fill(PAGE, "figures.counters.lead", **filled), **words)))  # fmt: skip
    return lookup(block, options, text.word(PAGE, "figures.counters.view"), axes, at)


def _plot_views(block: FigureBlock, build: tuple[Any, FigureText], long: pl.DataFrame,
                name: str, words: dict[str, str], shown: "Shown") -> FigureBlock:  # fmt: skip
    """F4 with its limits and a gene type as views (WO8): the build's 1.96 SD, every gene."""
    part = long.filter(pl.col("comparison") == name)
    options, at = [], []
    offered = views.types(part, [concordance.ALL])
    for i, empirical in enumerate((False, True)):
        for j, gene_type in enumerate(offered):
            if (i, j) == (0, 0):
                options.append((_plot_name(gene_type, empirical), *build))
            else:
                data, filled = views.plot(part, gene_type, empirical, words, shown.prefix)
                lead = text.fill(PAGE, f"figures.agreement.lead{'_empirical' * empirical}",
                                 **filled)  # fmt: skip
                options.append((_plot_name(gene_type, empirical), data,
                                _text("agreement", lead, **words)))  # fmt: skip
            at.append((i, j))
    axes = [(text.word(PAGE, "figures.agreement.axis_limits"),
             [text.word(PAGE, f"figures.agreement.limits_{k}") for k in ("sd", "empirical")]),
            (text.word(PAGE, "figures.counters.axis_type"),
             [views.type_words(t) for t in offered])]  # fmt: skip
    return lookup(block, options, text.word(PAGE, "figures.agreement.view"), axes, at)


def _plot_name(gene_type: str | None, empirical: bool) -> str:
    """How the live line names one F4 view."""
    return text.fill(PAGE, "figures.agreement.option", genes=views.type_words(gene_type),
                     limits=text.word(PAGE, "figures.agreement.limits_"
                                      + ("empirical" if empirical else "sd")))  # fmt: skip


# The tables of what M3 measured, which a lite page keeps; the rest record its inputs (D-111).
RESULTS_TABLES = ("alignment", "strandedness", "concordance")


def page_table(key: str, rows: Sequence[Sequence[str]], numeric: tuple[int, ...],
               **values: Any) -> ResultTable:  # fmt: skip
    """A small table on the page, its title and headings from ``tables.<key>.page``."""
    return ResultTable(text.fill(PAGE, f"tables.{key}.page.title", **values),
                       text.words(PAGE, f"tables.{key}.page.columns"),
                       tuple(tuple(r) for r in rows), numeric=numeric,
                       results=key in RESULTS_TABLES)  # fmt: skip


def _share(value: float | None) -> str:
    """A share to three decimals, or the content's word for none."""
    return figures.number(value, ".3f")


def _pct(value: float | None) -> str:
    """A percent to one decimal, or the content's word for none."""
    return figures.number(value, ".1f")


def _count(value: float | None) -> str:
    """A whole number with thousands separators, or the content's word for none."""
    return text.word(PAGE, "undefined") if value is None else f"{value:,.0f}"


def alignment_table(found: "Frames") -> ResultTable:
    """Per sample: uniquely mapped, STAR's splices, and flagstat's records beside reads."""
    rows = [(r["sample"], _pct(r["unique_pct"]), _pct(r["splices_per_100_unique"]),
             _share(r["splices_annotated_share"]), _count(r["splices_noncanonical"]),
             _count(r["records"]), _pct(r["mapped_records_pct"]),
             _pct(r["primary_mapped_pct"]))
            for r in found.alignment.iter_rows(named=True)]  # fmt: skip
    return page_table("alignment", rows, tuple(range(1, 8)))


def contigs(joined: str) -> str:
    """A reference's contigs as the page lists them: by name when few, else how many."""
    names = joined.split(",")
    return ", ".join(names) if len(names) <= SHOWN_CONTIGS else text.fill(
        PAGE, "genome_contigs", n=len(names))  # fmt: skip


def steps_table(steps: pl.DataFrame) -> ResultTable:
    """The cluster's cost per tool step, from GNU time's files (ED.6e-7)."""
    rows = [(text.word(PAGE, f"cluster.steps.{r['step']}"), f"{r['files']:,}",
             f"{r['wall_s_sum'] / 3600:,.2f}", f"{r['wall_s_max']:,.0f}",
             f"{r['max_rss_bytes'] / 2**30:,.2f}")
            for r in steps.iter_rows(named=True)]  # fmt: skip
    return page_table("cluster_steps", rows, (1, 2, 3, 4))


def result_tables(found: "Frames") -> dict[str, ResultTable]:
    """The page's tables, at most eight columns each; the files hold every column."""
    ref = [(r["reference"], contigs(r["contigs"]), f"{r['length_nt']:,}",
            f"{r['genes_subset']:,} of {r['genes_genome']:,}",
            figures.number(r["atg_share"], ".4f"))
           for r in found.reference.iter_rows(named=True)]  # fmt: skip
    strand = [(r["sample"], f"{r['assigned_s1']:,.0f}", f"{r['assigned_s2']:,.0f}",
               _share(r["forward_share"]), r["call"], r["used"], _share(r["compare_bias"]),
               r["compare_call"] or text.word(PAGE, "undefined"))
              for r in found.strand.iter_rows(named=True)]  # fmt: skip
    agree = found.agreements.concordance.filter(pl.col("unit") == "all")
    rows = [(text.word(PAGE, f"comparisons.{r['comparison']}.name"), f"{r['n_genes']:,}",
             f"{r['n_identical']:,}", _share(r["spearman"]),
             f"{_share(r['spearman_low'])} to {_share(r['spearman_high'])}", _share(r["bias"]),
             f"{_share(r['loa_low'])} to {_share(r['loa_high'])}")
            for r in agree.iter_rows(named=True)]  # fmt: skip
    out = {
        "reference": page_table("reference", ref, (2,)),
        "strandedness": page_table("strandedness", strand, (1, 2, 3, 6)),
        "concordance": page_table("concordance", rows, (1, 2, 3, 5)),
    }
    out["alignment"] = alignment_table(found)
    if found.folders is not None:
        out["folders"] = page_table("folders", [(
            text.word(PAGE, f"folders.{r['step']}"), r["folder"],
            text.word(PAGE, BUILT[r["built"]])) for r in found.folders.iter_rows(named=True)],
            ())  # fmt: skip
    if found.incomplete is not None:
        r = found.incomplete.row(0, named=True)
        out["incomplete"] = page_table("incomplete", [(
            f"{r['genes_kept']:,}", f"{r['genes_removed']:,}", f"{r['kept_reads_main']:,.0f}",
            f"{r['kept_reads_cut']:,.0f}", f"{r['gained']:,.0f}", f"{r['lost']:,.0f}",
            f"{r['removed_reads_main']:,.0f}", f"{r['unique_pct_main']:.1f} to "
            f"{r['unique_pct_cut']:.1f}")], tuple(range(7)))  # fmt: skip
    if found.context is not None:
        out["context"] = page_table("context", [(
            r["sample"], _pct(r["salmon_percent_mapped"]), _share(r["expected_share"]),
            _share(r["expected_read_share"]), _share(r["observed_assigned_share"]),
            figures.number(r["observed_over_expected"], ".2f"))
            for r in found.context.iter_rows(named=True)], (1, 2, 3, 4, 5))  # fmt: skip
    return out


def context_note(context: pl.DataFrame) -> str:
    """Salmon's median share beside the share assigned here, and how often the latter is higher.

    On a partial reference the share assigned here can exceed Salmon's: reads from contigs
    the genome lacks can align to similar genes on it and be counted there (D-039).
    """
    ratio = context["observed_over_expected"].drop_nulls()
    return text.fill(
        PAGE, "tables.context.page.note",
        share=figures.number(figures.stat(context["expected_read_share"], "median") * 100, ".1f"),
        n_above=int((ratio > 1).sum()), n=context.height,
        low=figures.number(figures.stat(ratio, "min") if ratio.len() else None, ".2f"),
        high=figures.number(figures.stat(ratio, "max") if ratio.len() else None, ".2f"),
    )  # fmt: skip
