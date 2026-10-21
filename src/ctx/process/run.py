"""M3 (`ctx process`): trim, align and count every sample, and test the settings that move counts.

Order of operations (the M3 plan, section 7): every input problem is found and listed before
anything is aligned, the reference prepared on the way since its contigs and start codons are
among them; its index is built or reused; each run is capped and trimmed; every sample is
aligned at the main setting and the featured sample again under each setting; featureCounts
counts at -s 0, 1 and 2, the strand is called, and the counting variants follow at each
sample's called strand; the checks I1 to I6 must hold before any table, figure or page is
written (D-038). Tools write into content-keyed folders under ``tools/`` (D-037).
"""

import json
import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import date
from functools import partial
from pathlib import Path

import polars as pl

from ctx import runtime, tools
from ctx.dataset import Dataset, DatasetError, load
from ctx.de.run import Writer, newest
from ctx.gather import gdc
from ctx.gather.counts import refuse_sheet as refuse_matrix
from ctx.io.naming import MODULES, stem
from ctx.process import (
    DEFAULT_SETTINGS,
    MULTIMAP_RANGE,
    SETTINGS,
    STRAND,
    Aligner,
    ProcessError,
    Strandedness,
    align,
    context,
    count,
    reference,
)
from ctx.process import page as pages
from ctx.process import tables as tab
from ctx.qc import Fastp, QcError, reads
from ctx.qc import multiqc_io as mqc_io
from ctx.qc import registry as reg
from ctx.qc.checkpoint import gene_types
from ctx.report import text
from ctx.report.builder import ParameterCard, write

MODULE = PAGE = "process"
GTF_RESOURCE = "annotation_gtf"
START_CODONS = "reference_start_codons"
STRAND_METRIC = "gene_strand_bias"
UNSPLICED = "unspliced"
# The registry metrics this command reads, and those it reads with --compare-dataset.
REGISTRY_METRICS = (START_CODONS, STRAND_METRIC)
COMPARE_METRICS = (context.PERCENT, context.STRAND)


@dataclass(frozen=True)
class Options:
    """What `ctx process` was asked to do, beside the dataset.

    Attributes:
        dataset_file: The dataset file; its ``annotation_gtf`` resource is the annotation.
        genome: The genome FASTA, plain or gzipped.
        aligner: STAR or HISAT2 (D-036).
        strandedness: The strand to count at; auto calls it per sample from the reads.
        multimap: STAR's --outFilterMultimapNmax; None keeps STAR's 10.
        settings: The featured sample's extra alignments, as tokens of ``SETTINGS``.
        featured: The featured sample; None takes the sheet's first.
        reads_per_run: Reads kept per run before fastp; None keeps every read.
        compare: A dataset whose Salmon counts and metrics cover the same samples.
        fetchngs_dir: The outdir the sheet's FASTQ paths are relative to; None is
            ``<root>/data/raw/<dataset>``.
        registry: The QC registry; None resolves it.
        fastp: fastp's filter, M2's own, so M2 and M3 trim alike.
    """

    dataset_file: Path
    genome: Path
    aligner: Aligner = Aligner.STAR
    strandedness: Strandedness = Strandedness.AUTO
    multimap: int | None = None
    settings: tuple[str, ...] = (DEFAULT_SETTINGS,)
    featured: str | None = None
    reads_per_run: int | None = None
    compare: str | None = None
    fetchngs_dir: Path | None = None
    registry: Path | None = None
    fastp: Fastp = field(default_factory=Fastp)


@dataclass(frozen=True)
class Sheet:
    """The gather sheet, where it is, its readable runs' files and the runs left unread."""

    sheet: pl.DataFrame
    path: Path
    files: pl.DataFrame
    skipped: pl.DataFrame


@dataclass(frozen=True)
class Inputs:
    """Every input, checked: the dataset, its sheet and reads, the annotation and the reference.

    Attributes:
        genes: The whole annotation's gene records (``checkpoint.gene_types``), read once.
    """

    spec: Dataset
    sheet: pl.DataFrame
    sheet_path: Path
    files: pl.DataFrame
    skipped: pl.DataFrame
    gtf: Path
    genes: pl.DataFrame
    registry: reg.Registry
    featured: str
    compare: context.Tables | None
    reference: reference.Reference


@dataclass(frozen=True)
class Built:
    """A reference and its index, and how long the index took (0 when reused)."""

    reference: reference.Reference
    index: Path
    built: bool
    seconds: float


@dataclass(frozen=True)
class Aligned:
    """The trimmed runs, the main alignment, the featured sample's settings and their reports."""

    trimmed: align.Trimmed
    main: align.Alignment
    featured: dict[str, align.Alignment]
    reports: dict[str, mqc_io.MultiqcRun]


@dataclass(frozen=True)
class Counted:
    """Every count folder, its MultiQC runs, the strand calls and the strand each sample used."""

    plain: count.Counts
    variants: count.Counts
    cut: count.Counts | None
    runs: dict[str, mqc_io.MultiqcRun]
    counted: pl.DataFrame
    calls: pl.DataFrame
    used: dict[str, Strandedness]


def option_problems(options: Options) -> list[str]:
    """Problems in the options alone, all of them."""
    problems = []
    low, high = MULTIMAP_RANGE
    if options.multimap is not None and options.aligner is Aligner.HISAT2:
        problems.append("--multimap-max is STAR's --outFilterMultimapNmax; HISAT2 takes none")
    elif options.multimap is not None and not low <= options.multimap <= high:
        problems.append(f"--multimap-max {options.multimap} is outside {low} to {high}")
    if unknown := [s for s in options.settings if s not in SETTINGS]:
        problems.append(f"--settings: unknown {unknown}; known are {list(SETTINGS)}")
    if options.reads_per_run is not None and options.reads_per_run < 1:
        problems.append(f"--reads-per-run {options.reads_per_run} keeps no read")
    if not options.genome.is_file():
        problems.append(f"--genome {options.genome} is not a file")
    return problems


def annotation(root: Path, spec: Dataset) -> tuple[Path | None, list[str]]:
    """The dataset's annotation on disk, its SHA-256 the one the dataset file declares."""
    found = [r for r in spec.resources if r.name == GTF_RESOURCE]
    if not found:
        return None, [f"{spec.path.name} declares no {GTF_RESOURCE} resource"]
    path = root / found[0].dest
    if not path.is_file():
        return None, [f"the annotation is not at {path}; run ctx setup with this dataset file"]
    if (digest := reference.sha256(path)) != found[0].sha256:
        return None, [f"{path} has SHA-256 {digest}, not the {found[0].sha256} "
                      f"{spec.path.name} declares"]  # fmt: skip
    return path, []


def runs_of(path: Path) -> tuple[pl.DataFrame | None, list[str]]:
    """Gather's sheet at ``path`` as the runs ctx process aligns, or why not.

    A GDC sheet or a count matrix's sheet holds no reads and is refused; any other passes
    ``gdc.chosen``. The M3 preview reads its sheet here too, so it takes the runs this does.
    """
    sheet = pl.read_parquet(path)
    if refused := gdc.refuse_sheet(sheet, path, MODULE) or refuse_matrix(sheet, path, MODULE):
        return None, refused
    return gdc.chosen(sheet), []


def _sheet(root: Path, dataset: str, options: Options) -> tuple[Sheet | None, list[str]]:
    """The gather sheet and its readable runs, or why not; a mix of layouts is a problem.

    A mixed sheet is still returned with its problem, so the checks that read its samples
    (``--featured-sample``, ``--compare-dataset``) run beside it.
    """
    folder = root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    try:
        path = newest(folder, dataset, "samplesheet")
        sheet, refused = runs_of(path)
        if sheet is None:
            return None, refused
        files, skipped = reads.fastq_files(sheet, options.fetchngs_dir or root / "data" / "raw"
                                           / dataset)  # fmt: skip
    except (DatasetError, QcError) as error:
        return None, [str(error)]
    found = Sheet(sheet, path, files, skipped)
    mates = files.group_by("run_accession").len()["len"].unique().to_list()
    if len(mates) > 1:
        return found, ["the sheet mixes single-end and paired-end runs; one alignment takes one"]
    return found, []


def _spec(options: Options, dataset: str) -> tuple[Dataset | None, list[str]]:
    """The dataset file, read and checked, or why not; a file with another slug is a problem.

    A subset's parent file is read here too: the page names the parent by its label, and a
    missing one would otherwise stop the run after every alignment.
    """
    if not options.dataset_file.is_file():
        return None, [f"--dataset-file {options.dataset_file} is not a file"]
    try:
        spec = load(options.dataset_file)
        spec.parent()
    except DatasetError as error:
        return None, [str(error)]
    if spec.slug != dataset:
        return spec, [f"{options.dataset_file}: slug {spec.slug!r}, not {dataset!r}"]
    return spec, []


def _registry(options: Options) -> tuple[reg.Registry | None, list[str]]:
    """The QC registry, holding every metric this command reads, or why not."""
    try:
        registry = reg.load(reg.resolve(options.registry))
    except reg.RegistryError as error:
        return None, [str(error)]
    needed = [*REGISTRY_METRICS, *(COMPARE_METRICS if options.compare else ())]
    if missing := [m for m in needed if m not in registry.metrics]:
        return None, [f"the QC registry {registry.path} has no {missing}; ctx process reads them"]
    return registry, []


def _compare(root: Path, options: Options, samples: list[str], registry: reg.Registry | None,
             genes: pl.DataFrame | None) -> tuple[context.Tables | None, list[str]]:  # fmt: skip
    """The other dataset's tables and their problems; checked without a sheet too, for none."""
    if not options.compare or registry is None:
        return None, []
    tables, found = context.read(root, options.compare, samples, registry)
    if tables is not None and genes is not None:
        found += context.unknown_genes(tables, genes)
    return tables, found


def _probe(options: Options) -> list[str]:
    """STAR's refusal, naming D-036 and both routes, when STAR is asked for and cannot map."""
    if options.aligner is not Aligner.STAR:
        return []
    try:
        align.star_probe()
    except tools.ToolError as error:
        return [str(error)]
    return []


def resolve(settings: runtime.Settings, dataset: str, options: Options) -> Inputs:
    """Every input checked, every problem listed together before anything is aligned.

    The reference is prepared under ``tools/index/`` on the way, since its contigs and start
    codons are among the checks; nothing else is written.

    Raises:
        ProcessError: Listing every problem found.
    """
    problems = option_problems(options)
    spec, found = _spec(options, dataset)
    problems += found
    registry, found = _registry(options)
    problems += found
    gtf, found = annotation(settings.root, spec) if spec else (None, [])
    problems += found
    genes = gene_types(gtf) if gtf is not None else None
    sheet, found = _sheet(settings.root, dataset, options)
    problems += found
    samples = sheet.files["sample"].unique(maintain_order=True).to_list() if sheet else []
    featured = options.featured or (samples[0] if samples else "")
    if sheet and featured not in samples:
        problems.append(f"--featured-sample {featured} is not a sample with readable reads")
    compare, found = _compare(settings.root, options, samples, registry, genes)
    problems += found
    problems += _probe(options)
    main = None
    if gtf is not None and options.genome.is_file():
        main, found = _reference(options.genome, gtf, settings.root, registry)
        problems += found
    if (problems or spec is None or registry is None or sheet is None or gtf is None
            or genes is None or main is None):  # fmt: skip
        raise ProcessError("ctx process cannot start:\n  " + "\n  ".join(problems))
    return Inputs(spec, sheet.sheet, sheet.path, sheet.files, sheet.skipped, gtf, genes,
                  registry, featured, compare, main)  # fmt: skip


def _out(root: Path) -> Path:
    """The module's output folder."""
    return root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"


def _reference(genome: Path, gtf: Path, root: Path, registry: reg.Registry | None
               ) -> tuple[reference.Reference | None, list[str]]:  # fmt: skip
    """The whole reference, or why it cannot be used: no shared contig, or foreign start codons.

    Without a usable registry the contigs are still checked; its problem is listed already.
    """
    try:
        ref = reference.prepare(genome, gtf, _out(root) / "tools" / "index")
    except ProcessError as error:
        return None, [str(error)]
    if registry is None:
        return ref, []
    share, metric = ref.start_codons.share, registry.metrics[START_CODONS]
    if share is None or metric.band(share) == "bad":
        return None, [f"{ref.start_codons.atg} of {ref.start_codons.records} start codons read "
                      f"ATG ({share}), below the registry's usable band {metric.usable}: the "
                      "genome and the annotation are not one assembly (D-039)"]  # fmt: skip
    return ref, []


def build_index(ref: reference.Reference, aligner: Aligner, parent: Path, threads: int) -> Built:
    """The aligner's index of a reference, timed."""
    started = time.perf_counter()
    builder = reference.star_index if aligner is Aligner.STAR else reference.hisat2_index
    index, built = builder(ref, parent, threads)
    return Built(ref, index, built, time.perf_counter() - started if built else 0.0)


def align_all(inputs: Inputs, options: Options, indices: Mapping[str, Built],
              threads: int, config: Path) -> Aligned:  # fmt: skip
    """Cap and trim every run, align every sample, then the featured sample per setting."""
    root = config.parent
    tools_ds = root / "tools" / inputs.spec.slug
    trimmed = align.trim(inputs.files, tools_ds, fastp=options.fastp,
                         cap=options.reads_per_run, threads=threads)  # fmt: skip
    runs = align.runs(inputs.files)
    main_reads = align.by_sample(inputs.files, {r: trimmed.reads(r) for r in runs})
    main = align.align(main_reads, tools_ds, align.AlignOptions(
        options.aligner, indices["main"].index, threads, multimap=options.multimap))  # fmt: skip
    featured = {inputs.featured: main_reads[inputs.featured]}
    jobs: dict[str, tuple[dict[str, list[tuple[Path, ...]]], align.AlignOptions]] = {}
    base = align.AlignOptions(options.aligner, indices["main"].index, threads, options.multimap)
    if "cut" in options.settings:
        jobs["cut"] = (featured, align.AlignOptions(options.aligner, indices["cut"].index,
                                                    threads, options.multimap))  # fmt: skip
    if "untrimmed" in options.settings:
        raw = align.by_sample(inputs.files, {r: trimmed.raw(r) for r in runs})
        jobs["untrimmed"] = ({inputs.featured: raw[inputs.featured]}, base)
    if options.aligner is Aligner.HISAT2:
        jobs[UNSPLICED] = (featured, align.AlignOptions(Aligner.HISAT2, indices["main"].index,
                                                        threads, spliced=False))  # fmt: skip
    done = {name: align.align(given, tools_ds, setting) for name, (given, setting) in jobs.items()}
    return Aligned(trimmed, main, done, summaries(inputs, config, main, done, trimmed, threads))


def summaries(inputs: Inputs, config: Path, main: align.Alignment,  # noqa: PLR0913 - one each
              featured: Mapping[str, align.Alignment], trimmed: align.Trimmed, jobs: int
              ) -> dict[str, mqc_io.MultiqcRun]:  # fmt: skip
    """The main MultiQC run, with fastp's reports, and one per featured setting.

    ``jobs`` of them run at a time (``mqc_io.at_once``).
    """
    folder = config.parent / "multiqc" / inputs.spec.slug
    work = {"main": partial(align.summarise, main, folder / "main", config, trimmed)}
    work |= {name: partial(align.summarise, a, folder / name, config)
             for name, a in featured.items()}  # fmt: skip
    return mqc_io.at_once(work, jobs)


def multiqc_config(writer: Writer, registry: reg.Registry) -> Path:
    """The MultiQC config the registry implies, raw counts unscaled, beside the tables."""
    config = writer.out / f"{writer.stem('multiqc_config')}.yml"
    writer.out.mkdir(parents=True, exist_ok=True)
    settings = reg.multiqc_config(registry, text.word(PAGE, "multiqc.comment"))
    config.write_text(json.dumps(settings, indent=1) + "\n")
    writer.files.append(config)
    return config


def count_all(inputs: Inputs, options: Options, indices: Mapping[str, Built],
              aligned: Aligned, threads: int, config: Path) -> Counted:  # fmt: skip
    """Count at -s 0, 1 and 2, call the strand, then count the variants at each sample's."""
    tools_ds = config.parent / "tools" / inputs.spec.slug
    folder = config.parent / "multiqc" / inputs.spec.slug
    main, gtf = aligned.main, indices["main"].reference.gtf
    plain = count.count(main, gtf, tools_ds, count.plain_settings(main), threads)
    runs = count.summaries(plain, main, folder, config, threads)
    counted = count.assigned(runs)
    calls = count.strand_calls(counted, inputs.registry.metrics[STRAND_METRIC])
    used = count.called_strands(calls, options.strandedness)
    variants = count.count(main, gtf, tools_ds, count.variant_settings(main, used), threads)
    variant_runs = count.summaries(variants, main, folder, config, threads)
    cut = None
    if "cut" in aligned.featured:
        half, sample = aligned.featured["cut"], inputs.featured
        options_cut = count.FcOptions((STRAND[used[sample]][1],), half.paired)
        cut = count.count(half, indices["cut"].reference.gtf, tools_ds, {"cut": options_cut},
                          threads)  # fmt: skip
        variant_runs |= count.summaries(cut, half, folder, config, threads)
    counted = pl.concat([counted, count.assigned(variant_runs)])
    return Counted(plain, variants, cut, {**runs, **variant_runs}, counted, calls, used)


def checks(inputs: Inputs, indices: Mapping[str, Built], aligned: Aligned,
           counted: Counted) -> None:  # fmt: skip
    """I1 to I6 on the main alignment, the featured settings and every count folder (D-038).

    Raises:
        ProcessError: Listing every failed check.
    """
    main, report = aligned.main, aligned.reports["main"]
    fastp = align.fastp_reads(report, list(align.runs(inputs.files)))
    files = inputs.files
    mine = files.filter(pl.col("sample") == inputs.featured)
    found = []
    for name, alignment in {"main": main, **aligned.featured}.items():
        given = files if name == "main" else mine
        column = "raw_reads" if name == "untrimmed" else "passed_reads"
        said = align.aligner_reads(aligned.reports[name], alignment)
        expected = count.expected_input(fastp, given, column, paired=alignment.paired)
        found += count.input_problems(said, expected)
        found += count.primary_problems(aligned.reports[name], said, alignment)
    if main.aligner is Aligner.STAR:
        found += count.star_parser_problems(main, report, align.aligner_reads(report, main))
    folders = [(counted.plain, main, indices["main"]), (counted.variants, main, indices["main"])]
    if counted.cut is not None:
        folders.append((counted.cut, aligned.featured["cut"], indices["cut"]))
    for counts, alignment, built in folders:
        found += count.assigned_problems(counts, counted.counted, alignment)
        found += count.order_problems(counts, alignment)
        found += count.gene_set_problems(built.reference.gtf, tab.matrices(counts))
    count.check(found)


def run(settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str,
        dataset: str, options: Options) -> list[Path]:  # fmt: skip
    """Trim, align and count the dataset's reads, check them, and write the M3 page and files.

    Returns:
        The page, then every table and figure file.

    Raises:
        ProcessError: Listing every input problem, when nothing has been aligned and only the
            reference under ``tools/index/`` may have been prepared; or every failed check
            after the alignments, when the tool folders, the MultiQC reports and the MultiQC
            config are written and no table, figure or page is.
        ToolError: If a tool fails.
        QcError: If MultiQC fails.
    """
    inputs = resolve(settings, dataset, options)
    parent, threads = _out(settings.root) / "tools" / "index", settings.threads
    indices = {"main": build_index(inputs.reference, options.aligner, parent, threads)}
    if "cut" in options.settings:
        half = reference.prepare(options.genome, inputs.gtf, parent, half=True)
        indices["cut"] = build_index(half, options.aligner, parent, threads)
    day = date.today()
    writer = Writer(_out(settings.root), dataset, day, releases(inputs, indices, options))
    tools_readme(writer.out)
    config = multiqc_config(writer, inputs.registry)
    aligned = align_all(inputs, options, indices, threads, config)
    counted = count_all(inputs, options, indices, aligned, threads, config)
    checks(inputs, indices, aligned, counted)
    page = pages.build(pages.Context(
        writer=writer, inputs=inputs, options=options, indices=indices, aligned=aligned,
        counted=counted, genome_genes=inputs.genes, cards=cards, command=command,
        settings=settings))  # fmt: skip
    out = writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"
    return [write(page, out), *writer.files]


def releases(inputs: Inputs, indices: Mapping[str, Built], options: Options) -> dict[str, str]:
    """What every README records: genome, annotation, index, tool pins and inputs."""
    ref = inputs.reference
    tool_names = ("fastp", "seqkit", "samtools", "featureCounts", "multiqc",
                  "STAR" if options.aligner is Aligner.STAR else "hisat2")  # fmt: skip
    release = ref.annotation.release or text.word(PAGE, "no_release")
    return {
        "genome": f"{options.genome.name}, SHA-256 {ref.source.sha256}",
        "annotation": f"{inputs.gtf.name}, {release}, SHA-256 {ref.annotation_source_sha256}",
        **{f"index:{k}": b.index.name for k, b in indices.items()},
        **{f"tool:{t}": tools.TOOLS[t].version for t in tool_names},
        "input:samplesheet": inputs.sheet_path.name,
        **({f"input:{inputs.compare.dataset}": ", ".join(p.name for p in inputs.compare.sources)}
           if inputs.compare else {}),
    }  # fmt: skip


def tools_readme(out: Path) -> None:
    """``tools/README.md``: what the content-keyed folders hold and how they are keyed (D-037)."""
    folder = out / "tools"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "README.md").write_text(text.word(PAGE, "tools_readme") + "\n")
