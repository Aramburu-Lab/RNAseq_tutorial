"""The M3 page and its files: tables, figures, key numbers and the nine sections (D-045).

Every word is ``content/process.toml``'s; every number comes from ``ctx.process.tables``.
A page built for a declared subset (``[dataset] subset_of``) says first that it is the
optional laptop exercise (D-050). What the page shows of each frame is ``ctx.process.blocks``.
"""

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import polars as pl

from ctx import runtime, tools
from ctx.de.run import Writer
from ctx.process import (
    STRAND,
    Aligner,
    align,
    blocks,
    cluster,
    context,
    figure_tables,
    figures,
    reference,
)
from ctx.process import tables as tab
from ctx.qc.checkpoint import gene_types
from ctx.qc.registry import Metric
from ctx.report import text
from ctx.report.builder import (
    FigureBlock,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
)
from ctx.viz.style import MAX_CATEGORIES

if TYPE_CHECKING:
    from ctx.process.stages import Aligned, Built, Counted, Inputs, Options

PAGE = MODULE = "process"
# Sheet columns the bars may filter by, where the aligned samples hold two levels or more.
GROUPS = ("cell_type", "condition")
SHOWN_RUNS = 10
# The content table whose Input and Method a cluster's page reads in place of the page's own.
CLUSTER = "variants.cluster"
# The derived multimapper limit's row in the settings table (D-041).
LIMIT_ONE = "limit_1"


@dataclass(frozen=True)
class Context:
    """Everything one run hands the page."""

    writer: Writer
    inputs: "Inputs"
    options: "Options"
    indices: Mapping[str, "Built"]
    aligned: "Aligned"
    counted: "Counted"
    genome_genes: pl.DataFrame
    cards: tuple[ParameterCard, ...]
    command: str
    settings: runtime.Settings


@dataclass(frozen=True)
class Shown:
    """What the page shows beside the tables: the aligner, the featured sample, the groups.

    Attributes:
        aligner: Its outcomes are the mapping figures' series.
        featured: The sample aligned again under each setting.
        groups: sample, then each sheet column the bars filter by.
        genome_name: The reference as the text names it.
        strand_metric: The registry's ``gene_strand_bias``, whose bands the caption gives.
        genome_genes: Genes in the whole annotation file.
        partial: Whether the genome lacks a contig the annotation names, so genome-wide
            mapping bands do not apply.
        prefix: Put before every caption; the preview's says its data are made up.
        cluster: Whether the reads were aligned on the cluster (``--align-dir``).
    """

    aligner: Aligner
    featured: str
    groups: pl.DataFrame
    genome_name: str
    strand_metric: Metric
    genome_genes: int
    partial: bool = True
    prefix: str = ""
    cluster: bool = False

    @property
    def usable(self) -> list[str]:
        """The group columns the bars filter by."""
        return self.groups.columns[1:]


@dataclass(frozen=True)
class Frames:
    """The run's tables, before any is written.

    Attributes:
        runs: Each run's reads before and after fastp (``align.fastp_reads``), or None.
        folders: Each keyed tool folder the run used and whether it built it, or None.
    """

    alignment: pl.DataFrame
    outcomes: dict[str, dict[str, int]]
    settings: dict[str, dict[str, int]]
    strand: pl.DataFrame
    counting: pl.DataFrame
    gene_counts: pl.DataFrame
    comparison: pl.DataFrame
    agreements: tab.Agreements
    reference: pl.DataFrame
    incomplete: pl.DataFrame | None
    context: pl.DataFrame | None
    runs: pl.DataFrame | None = None
    folders: pl.DataFrame | None = None


def genome_name(ref: reference.Reference) -> str:
    """The genome as the text names it: its contigs when few, else how many there are."""
    contigs = ref.genome.contigs
    if len(contigs) <= blocks.SHOWN_CONTIGS:
        return ", ".join(contigs)
    return text.fill(PAGE, "genome_contigs", n=len(contigs))


def setting_reads(ctx: Context, rows: Mapping[str, Mapping[str, float]]
                  ) -> dict[str, dict[str, int]]:  # fmt: skip
    """The featured sample's outcomes under each setting, in the page's order.

    Limit 1 is derived from the main row (``align.derived_limit_one``), under STAR only.
    """
    main, featured = ctx.aligned.main, ctx.inputs.featured
    out = {"main": tab.outcome_reads(rows[featured], main.aligner, main.paired)}
    if main.aligner is Aligner.STAR:
        out[LIMIT_ONE] = tab.outcome_reads(align.derived_limit_one(rows[featured]),
                                           main.aligner, main.paired)  # fmt: skip
    for name, alignment in ctx.aligned.featured.items():
        row = tab.aligner_rows(ctx.aligned.reports[name], alignment)[featured]
        out[name] = tab.outcome_reads(row, alignment.aligner, alignment.paired)
    return out


def frames(ctx: Context) -> Frames:
    """Every table of the run, from the readers of ``align`` and ``count``."""
    main, report = ctx.aligned.main, ctx.aligned.reports["main"]
    rows = tab.aligner_rows(report, main)
    fastp = align.fastp_reads(report, list(align.runs(ctx.inputs.files)))
    alignment = tab.alignment_table(rows, main, ctx.inputs.files, fastp, tab.flagstat(report))
    matrices = {**tab.matrices(ctx.counted.plain), **tab.matrices(ctx.counted.variants)}
    used, ref = ctx.counted.used, ctx.indices["main"].reference
    genes = gene_types(ref.gtf)
    called = tab.at_called(matrices, used)
    inputs = dict(alignment.select("sample", "input_reads").rows())
    salmon = None
    if ctx.inputs.compare is not None:
        salmon = context.expected(ctx.inputs.compare, ctx.genome_genes, ref.genome.contigs,
                                  ctx.inputs.registry)  # fmt: skip
    star = tab.star_columns(main) if main.aligner is Aligner.STAR else None
    long = tab.comparisons(matrices, used, main if star is not None else None)
    counting = tab.counting(ctx.counted.counted, matrices, inputs)
    return Frames(
        alignment=alignment, outcomes={s: tab.outcome_reads(rows[s], main.aligner, main.paired)
                                       for s in main.samples},
        settings=setting_reads(ctx, rows),
        strand=tab.strandedness(ctx.counted.counted, ctx.counted.calls, used, star, salmon),
        counting=counting, gene_counts=tab.gene_counts(called, genes),
        comparison=long.join(genes.select("gene_id", "gene_name", "gene_type"), on="gene_id"),
        agreements=tab.agreements(long, genes), reference=reference_table(ctx),
        incomplete=incomplete(ctx, called), context=context_table(salmon, counting, used),
        runs=fastp, folders=folder_table(ctx),
    )  # fmt: skip


def context_table(salmon: pl.DataFrame | None, counting: pl.DataFrame,
                  used: Mapping[str, Any]) -> pl.DataFrame | None:  # fmt: skip
    """Salmon's expected share beside the share of input reads assigned here, and their ratio.

    The ratio is no test against a ceiling: on a partial reference the share assigned here
    can exceed Salmon's (D-039).
    """
    if salmon is None:
        return None
    observed = {s: counting.filter((pl.col("sample") == s) & (pl.col("counter") ==
                                   f"s{STRAND[u][1]}"))["assigned_pct"].item() / 100
                for s, u in used.items()}  # fmt: skip
    return salmon.with_columns(observed_assigned_share=pl.col("sample").replace_strict(
        observed, return_dtype=pl.Float64)).with_columns(
        observed_over_expected=pl.when(pl.col("expected_read_share") > 0).then(
            pl.col("observed_assigned_share") / pl.col("expected_read_share")))  # fmt: skip


def incomplete(ctx: Context, called: pl.DataFrame) -> pl.DataFrame | None:
    """The featured sample's genes on the whole and on the half reference, when cut was run."""
    cut = ctx.counted.cut
    if cut is None:
        return None
    sample = ctx.inputs.featured
    half = tab.matrices(cut)["cut"].select("gene_id", align.bam_name(sample))
    unique = []
    for name, alignment in (("main", ctx.aligned.main), ("cut", ctx.aligned.featured["cut"])):
        row = tab.aligner_rows(ctx.aligned.reports[name], alignment)[sample]
        reads = tab.outcome_reads(row, alignment.aligner, alignment.paired)
        unique.append(100 * reads["unique"] / reads["input"] if reads["input"] else 0.0)
    return tab.incomplete(called[sample], half, called, (unique[0], unique[1]))


def reference_table(ctx: Context) -> pl.DataFrame:
    """One row per reference: the genome, the annotation subset, the start codons, the index."""
    whole = len(ctx.genome_genes)
    rows = []
    for name, built in ctx.indices.items():
        ref = built.reference
        types = gene_types(ref.gtf)
        star = ctx.options.aligner is Aligner.STAR
        rows.append({
            "reference": text.word(PAGE, f"references.{name}"),
            "file": ctx.inputs.genome_name, "sha256": ref.source.sha256,
            "bytes": ref.source.bytes, "contigs": ",".join(ref.genome.contigs),
            "length_nt": ref.genome.length, "n_bases": ref.genome.n_bases,
            "sa_index_bases": reference.sa_index_bases(ref.genome.length) if star else None,
            "genes_subset": ref.annotation.genes,
            "protein_coding_subset": types.filter(pl.col("gene_type") == "protein_coding").height,
            "genes_genome": whole, "annotation_release": ref.annotation.release,
            "annotation_sha256": ref.annotation_source_sha256,
            "start_codons": ref.start_codons.records, "atg_share": ref.start_codons.share,
            "index_key": built.index.name, "built": built.built, "build_s": built.seconds,
        })  # fmt: skip
    # Typed whatever the aligner: SAindexNbases is null under HISAT2.
    return pl.DataFrame(rows, infer_schema_length=None).with_columns(
        pl.col("sa_index_bases").cast(pl.Int64), pl.col("atg_share").cast(pl.Float64)
    )


def folder_table(ctx: Context) -> pl.DataFrame:
    """Each content-keyed folder the run read, in the order it ran, and whether it built it.

    A rerun with one counting option changed reuses the trim and the alignments (D-037).
    """
    aligned, counted = ctx.aligned, ctx.counted
    rows = [(f"index_{k}", b.index.name, b.built) for k, b in ctx.indices.items()]
    rows.append(("trim", aligned.trimmed.folder.name, aligned.trimmed.built))
    rows += [(f"align_{k}", a.folder.name, a.built)
             for k, a in {"main": aligned.main, **aligned.featured}.items()]  # fmt: skip
    counts = {"plain": counted.plain, "variants": counted.variants, "cut": counted.cut}
    rows += [(f"count_{k}", c.folder.name, c.built) for k, c in counts.items() if c is not None]
    return pl.DataFrame(rows, schema={"step": pl.String, "folder": pl.String,
                                      "built": pl.Boolean}, orient="row")  # fmt: skip


def write_tables(writer: Writer, found: Frames, skipped: pl.DataFrame) -> None:
    """Every table as Parquet, CSV and README, its columns described in the content."""
    named = {"alignment": found.alignment, "settings": settings_table(found.settings),
             "strandedness": found.strand, "counting": found.counting,
             "gene_comparison": found.comparison, "concordance": found.agreements.concordance,
             "discordant": found.agreements.discordant, "reference": found.reference,
             "skipped": skipped}  # fmt: skip
    optional = {"incomplete": found.incomplete, "context": found.context,
                "folders": found.folders}  # fmt: skip
    named |= {k: v for k, v in optional.items() if v is not None}
    for what, frame in named.items():
        writer.table(frame, what, described(f"tables.{what}", frame.columns))
    samples = found.gene_counts.columns[3:]
    columns = described("tables.gene_counts", found.gene_counts.columns[:3])
    columns |= {s: ("reads", text.fill(PAGE, "tables.gene_counts.sample", sample=s))
                for s in samples}  # fmt: skip
    writer.table(found.gene_counts, "gene_counts", columns)


def described(key: str, columns: Sequence[str]) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of each column from ``<key>.units`` and ``<key>.columns``."""
    units = text.raw(PAGE, f"{key}.units")
    return {c: (units.get(c, "none"), text.word(PAGE, f"{key}.columns.{c}")) for c in columns}


def settings_table(found: Mapping[str, Mapping[str, int]]) -> pl.DataFrame:
    """The featured sample per setting: input reads and the shared outcomes, reads and %."""
    rows = []
    for name, reads in found.items():
        both = tab.shared(reads)
        rows.append({"setting": name, "derived": name == LIMIT_ONE, "input_reads": both["input"],
                     **{k: both[k] for k in ("unique", "multi", "unmapped")},
                     **{f"{k}_pct": 100 * both[k] / both["input"] if both["input"] else None
                        for k in ("unique", "multi", "unmapped")}})  # fmt: skip
    return pl.DataFrame(rows)


def groups_of(sheet: pl.DataFrame, samples: Sequence[str]) -> pl.DataFrame:
    """The aligned samples and the sheet columns worth filtering their bars by.

    A column is worth it with two levels at least and no more than a figure has colours; a
    missing value is shown as a word.
    """
    kept = sheet.unique("sample", keep="first").filter(pl.col("sample").is_in(samples))
    usable = [g for g in GROUPS if g in kept.columns
              and 1 < kept[g].n_unique() <= MAX_CATEGORIES]  # fmt: skip
    missing = text.word(PAGE, "figures.no_value")
    return kept.select("sample", *(pl.col(g).cast(pl.String).fill_null(missing) for g in usable))


def notes(ctx: Context, found: Frames) -> list[Result]:
    """What the page says before any result: a subset's page, the cap, the runs left unread."""
    said: list[Result] = []
    if (parent := ctx.inputs.spec.parent()) is not None:
        said.append(text.fill(PAGE, "notes.subset", parent=parent.name))
    if (cap := ctx.options.reads_per_run) is not None and found.runs is not None:
        said.append(capped(cap, found.runs, ctx.aligned.main.paired))
    if ctx.inputs.skipped.height:
        said.append(text.fill(PAGE, "notes.skipped", n=ctx.inputs.skipped.height,
                              runs=", ".join(ctx.inputs.skipped["run_accession"]
                                             .head(SHOWN_RUNS))))  # fmt: skip
    return said


def capped(cap: int, runs: pl.DataFrame, paired: bool) -> str:
    """The cap and the reads each run kept under it: fastp's reads before filtering.

    fastp counts both mates of a pair, the cap one pair, so a paired run's count is halved.
    """
    kept = runs["raw_reads"] // (2 if paired else 1)
    return text.fill(PAGE, "notes.capped", cap=f"{cap:,}", seed=align.SEED, n=kept.len(),
                     low=f"{figures.stat(kept, 'min'):,.0f}",
                     high=f"{figures.stat(kept, 'max'):,.0f}",
                     median=figures.reads(figures.stat(kept, "median")))  # fmt: skip


def results(writer: Writer, found: Frames, shown: Shown, first: Sequence[Result]
            ) -> tuple[Result, ...]:  # fmt: skip
    """The Results section in order: ``first``, the reference, then each group of figures."""
    tables = blocks.result_tables(found)
    align_items: list[FigureBlock | ResultTable | str] = [
        blocks.mapping_figure(writer, found, shown), tables["alignment"],
        text.word(PAGE, "tables.alignment.page.note"),
        blocks.settings_figure(writer, found, shown),
    ]  # fmt: skip
    align_items += [tables["incomplete"]] if "incomplete" in tables else []
    if found.context is not None:
        align_items += [tables["context"], blocks.context_note(found.context)]
    strand_items: list[FigureBlock | ResultTable | str] = [
        blocks.strand_figure(writer, found, shown),
        blocks.share_figure(writer, found, shown),
        tables["strandedness"],
    ]
    scatters, plots = blocks.agreement_figures(writer, found, shown)
    called = found.gene_counts.select(pl.sum_horizontal(found.gene_counts.columns[3:]) > 0).sum()
    lead = text.fill(PAGE, "groups.counters_lead", with_reads=f"{called.item():,}",
                     subset=f"{found.gene_counts.height:,}",
                     genome=f"{shown.genome_genes:,}")  # fmt: skip
    return (
        *first, tables["reference"], *([tables["folders"]] if "folders" in tables else []),
        Group(text.word(PAGE, "groups.align"), tuple(align_items)),
        Group(text.word(PAGE, "groups.strand"), tuple(strand_items)),
        Group(text.word(PAGE, "groups.counters"), tuple(scatters), lead=lead,
              tabs=len(scatters) > 1),
        Group(text.word(PAGE, "groups.agreement"), tuple(plots), tabs=len(plots) > 1),
        tables["concordance"],
    )  # fmt: skip


def key_numbers(found: Frames) -> tuple[KeyNumber, ...]:
    """Samples aligned, median uniquely mapped, the strand call, genes the counters disagree on."""
    alignment, strand = found.alignment, found.strand
    unique = alignment.sort("unique_pct").row(0, named=True)
    median = figures.stat(alignment["unique_pct"], "median")
    call, n_call = Counter(strand["call"]).most_common(1)[0]
    multimap = found.agreements.concordance.filter(
        (pl.col("comparison") == "multimap") & (pl.col("unit") == "all")
    ).row(0, named=True)
    return (
        KeyNumber(text.word(PAGE, "keys.samples"), f"{alignment.height:,}", text.fill(
            PAGE, "keys.samples_detail",
            reads=figures.reads(figures.stat(alignment["input_reads"], "median")))),
        KeyNumber(text.word(PAGE, "keys.unique"), f"{median:.1f} %",
                  text.fill(PAGE, "keys.unique_detail", low=f"{unique['unique_pct']:.1f}",
                            sample=unique["sample"])),
        KeyNumber(text.word(PAGE, "keys.strand"), call,
                  text.fill(PAGE, "keys.strand_detail", n_call=n_call, n=strand.height)),
        KeyNumber(text.word(PAGE, "keys.differ"),
                  f"{multimap['n_genes'] - multimap['n_identical']:,}",
                  text.fill(PAGE, "keys.differ_detail", n=f"{multimap['n_genes']:,}")),
    )  # fmt: skip


def partial(found: Frames) -> bool:
    """Whether the main reference lacks a gene of the annotation, so a contig is missing."""
    main = found.reference.row(0, named=True)
    return bool(main["genes_subset"] < main["genes_genome"])


def prose(shown: Shown, values: Mapping[str, Any]) -> dict[str, str]:
    """Input, Method and Interpretation filled; the reference note says whether bands apply."""
    names = {"genome_name": shown.genome_name, "aligner": str(shown.aligner).upper()}
    note = text.fill(PAGE, f"reference_note.{'partial' if shown.partial else 'whole'}", **names)
    where = f"{CLUSTER}." if shown.cluster else ""
    return {
        "input": text.fill(PAGE, f"{where}input", **values, **names),
        "method": text.fill(PAGE, f"{where}method", **names),
        "interpretation": text.fill(PAGE, "interpretation", **values, **names,
                                    reference_note=note),
    }  # fmt: skip


def sections(ctx: Context, found: Frames, shown: Shown) -> dict[str, Any]:
    """The fixed sections, those that state this run's numbers filled from it."""
    ref = ctx.indices["main"].reference
    low, high = (figures.stat(found.alignment["unique_pct"], h) for h in ("min", "max"))
    filled = prose(shown, {
        "samples": found.alignment.height, "runs": ctx.inputs.files["run_accession"].n_unique(),
        "genome": ctx.inputs.genome_name, "length": f"{ref.genome.length:,}",
        "release": ref.annotation.release or text.word(PAGE, "no_release"),
        "low": f"{low:.1f}", "high": f"{high:.1f}",
    })  # fmt: skip
    return text.skeleton(PAGE, CLUSTER if shown.cluster else None) | filled


def build(ctx: Context) -> ModulePage:
    """Write every table and figure, and return the page."""
    found = frames(ctx)
    write_tables(ctx.writer, found, ctx.inputs.skipped)
    ref = ctx.indices["main"].reference
    groups = groups_of(ctx.inputs.sheet, ctx.aligned.main.samples)
    shown = Shown(
        ctx.options.aligner,
        ctx.inputs.featured,
        groups,
        genome_name(ref),
        ctx.inputs.registry.metrics["gene_strand_bias"],
        len(ctx.genome_genes),
        partial=partial(found),
        cluster=ctx.inputs.source is not None,
    )
    # What the figures draw that no table above holds, so ctx figure redraws them (F.1).
    figure_tables.write(ctx.writer, found, shown, ctx.inputs.registry)
    body = results(ctx.writer, found, shown, notes(ctx, found))
    said: tuple[str, ...] = ()
    if ctx.inputs.source is not None:
        steps = cluster.steps(cluster_dirs(ctx))
        ctx.writer.table(steps, "cluster_steps", described("tables.cluster_steps", steps.columns))
        body = (*body, blocks.steps_table(steps))
        said = cluster_notes(ctx, steps)
    usage, budget = runtime.measure(ctx.settings.started), runtime.budget(MODULE)
    aligner = "STAR" if ctx.options.aligner is Aligner.STAR else "hisat2"
    pins = {
        name: tools.TOOLS[name].version
        for name in ("fastp", "seqkit", aligner, "samtools", "featureCounts", "multiqc")
    }
    inputs = [ctx.inputs.sheet_path.name, ctx.inputs.genome_name, ctx.inputs.gtf.name]
    inputs += [p.name for p in ctx.inputs.compare.sources] if ctx.inputs.compare else []
    return ModulePage(
        title=text.fill(PAGE, "title", dataset=ctx.inputs.spec.name),
        **sections(ctx, found, shown), parameters=ctx.cards, results=body,
        key_numbers=key_numbers(found),
        reproduce=Reproduction(
            command=ctx.command, wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self, peak_rss_children=usage.peak_rss_children,
            profile=f"{ctx.settings.profile}, {ctx.settings.threads} threads",
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions={**runtime.versions(), **pins}, inputs=tuple(inputs), notes=said,
        ),
    )  # fmt: skip


def cluster_dirs(ctx: Context) -> list[Path]:
    """Each folder of cluster time files the page reads: every sample's, then the featured's."""
    source = ctx.inputs.source
    if source is None:
        return []
    return [
        *(source.sample_dir(s) for s in ctx.aligned.main.samples),
        *(source.featured_dir(ctx.inputs.featured, s) for s in ctx.aligned.featured),
    ]


def cluster_notes(ctx: Context, steps: pl.DataFrame) -> tuple[str, ...]:
    """The Reproduce section's word on the cluster: versions.txt's tools, the cost, the blocks."""
    source = ctx.inputs.source
    if source is None:
        return ()
    versions = ", ".join(f"{t} {v}" for t, v in source.versions.items())
    wall = float(steps["wall_s_sum"].sum()) if steps.height else 0.0
    rss = figures.stat(steps["max_rss_bytes"], "max") if steps.height else 0.0
    return (text.fill(PAGE, "cluster.reproduce.tools", versions=versions, index=source.index),
            text.fill(PAGE, "cluster.reproduce.cost", files=int(steps["files"].sum()),
                      hours=f"{wall / 3600:.2f}", gib=f"{rss / 2**30:.2f}"),
            text.word(PAGE, "cluster.reproduce.blocks"))  # fmt: skip
