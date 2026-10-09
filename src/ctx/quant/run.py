"""M4 (`ctx quant`): what Salmon did with every sample's reads, and how they became gene counts.

Two paths. The cohort path reads an nf-core/rnaseq outdir and copy Q1. The laptop path
(``--genome``) runs Salmon here on M3's trimmed reads (``ctx.quant.salmon``,
``ctx.quant.laptop``), into content-keyed folders under ``outputs/04_quant/tools/``.

Order of operations: every input problem is found and listed before anything runs; on the
laptop path Salmon passes its probe, then the reference, the index and each sample's
quantification are built or reused; MultiQC reads Salmon's logs; the port of tximport
summarises every ``quant.sf`` there is; the checks of ``ctx.quant.checks`` must all hold
before any table, figure or page is written. Writes only under ``outputs/04_quant/``.
"""

from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

import polars as pl

from ctx import runtime, tools
from ctx.dataset import Contrast, Dataset, DatasetError, load
from ctx.de.run import Writer, newest
from ctx.gather import rnaseq
from ctx.io.naming import MODULES, stem
from ctx.process.reference import Reference, prepare
from ctx.process.run import annotation
from ctx.qc.checkpoint import gene_types
from ctx.quant import (
    CountsFrom,
    QuantError,
    SalmonOptions,
    checks,
    laptop,
    outdir,
    salmon,
    tables,
    tximport,
)
from ctx.quant import page as pages
from ctx.quant.salmon import Quantified
from ctx.report.builder import ParameterCard, write

MODULE = "quant"


@dataclass(frozen=True)
class Options:
    """What `ctx quant` was asked to do, beside the dataset.

    Attributes:
        dataset_file: The dataset file; its ``annotation_gtf`` gives gene types, its
            ``[de] featured`` contrast the fold-change figure.
        rnaseq_dir: The nf-core/rnaseq outdir ``ctx gather`` read; None on the laptop path.
        salmon_dir: Copy Q1, the salmon folder with every ``quant.sf``; None is
            ``<root>/data/raw/<dataset>_salmon`` where it exists, else no transcript figure.
        counts_from: Which gene counts ctx writes from its port of tximport.
        genome: The laptop path's genome FASTA, as ``ctx process`` was given it.
        salmon: The laptop path's Salmon setting.
    """

    dataset_file: Path
    rnaseq_dir: Path | None
    salmon_dir: Path | None = None
    counts_from: CountsFrom = CountsFrom.LENGTH_SCALED_TPM
    genome: Path | None = None
    salmon: SalmonOptions = field(default_factory=SalmonOptions)


@dataclass(frozen=True)
class Inputs:
    """Every input, checked.

    Attributes:
        samples: The sheet's samples, in sheet order, one row each.
        genes: The annotation's gene records (``checkpoint.gene_types``).
        counts: The length-scaled matrix as ``ctx gather`` reads it, integer; None on the
            laptop path, whose counts come from its own Salmon run.
        salmon_dir: Copy Q1, or None when it is not there or on the laptop path.
        contrast: The featured contrast, or None.
        reads: The laptop path's trimmed reads per sample; empty on the cohort path.
        trim: M3's trim folder those reads are in, or None.
    """

    spec: Dataset
    sheet_path: Path
    samples: pl.DataFrame
    rnaseq_dir: Path | None
    gtf: Path
    genes: pl.DataFrame
    counts: rnaseq.Counts | None
    salmon_dir: Path | None
    contrast: Contrast | None
    reads: laptop.Reads = field(default_factory=dict)
    trim: Path | None = None


@dataclass(frozen=True)
class Measured:
    """What the run measured before any check.

    Attributes:
        salmon: MultiQC's Salmon table; ``fates`` and ``sums`` are made from it.
        port: The port's gene tables, or None without any ``quant.sf``.
        counts: The length-scaled integer matrix the figures read.
        salmon_dir: The folder holding every ``quant.sf`` and the merged tables, or None.
        quantified: The laptop path's keyed folders, or None.
        ref: The laptop path's reference, or None.
    """

    salmon: pl.DataFrame
    fates: pl.DataFrame
    sums: pl.DataFrame
    port: tximport.GeneTables | None
    counts: rnaseq.Counts
    salmon_dir: Path | None
    quantified: Quantified | None = None
    ref: Reference | None = None


def _out(root: Path) -> Path:
    """The module's output folder."""
    return root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"


def _spec(options: Options, dataset: str) -> tuple[Dataset | None, list[str]]:
    """The dataset file, read and checked, or why not."""
    if not options.dataset_file.is_file():
        return None, [f"--dataset-file {options.dataset_file} is not a file"]
    try:
        spec = load(options.dataset_file)
    except DatasetError as error:
        return None, [str(error)]
    if spec.slug != dataset:
        return spec, [f"{options.dataset_file}: slug {spec.slug!r}, not {dataset!r}"]
    return spec, []


def _sheet(root: Path, dataset: str) -> tuple[Path | None, pl.DataFrame | None, list[str]]:
    """The newest gather sheet, one row per run, or why not."""
    try:
        path = newest(root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather", dataset,
                      "samplesheet")  # fmt: skip
    except DatasetError as error:
        return None, None, [str(error)]
    return path, pl.read_parquet(path), []


def _counts(folder: Path, samples: set[str]) -> tuple[rnaseq.Counts | None, list[str]]:
    """The rnaseq outdir's length-scaled matrix, checked as ``ctx gather`` checks it."""
    if not folder.is_dir():
        return None, [f"--rnaseq-dir {folder} is not a folder"]
    try:
        return rnaseq.read(folder, samples), []
    except rnaseq.RnaseqError as error:
        return None, [str(error)]


def salmon_dir(root: Path, dataset: str, given: Path | None) -> Path | None:
    """Copy Q1: the folder given, else the default where it exists, else None."""
    if given is not None:
        return given
    default = root / "data" / "raw" / f"{dataset}_salmon"
    return default if default.is_dir() else None


def q1_problems(folder: Path | None, samples: list[str]) -> list[str]:
    """Everything copy Q1 must hold that is missing: each sample's quant.sf, the tables."""
    if folder is None:
        return []
    if not folder.is_dir():
        return [f"--salmon-dir {folder} is not a folder"]
    problems = outdir.quant_problems(folder, samples)
    for what in outdir.GENE_TABLES:
        problems += outdir.table_problems(folder / outdir.MERGED.format(what), outdir.GENE_KEYS,
                                          samples)  # fmt: skip
    for what in outdir.TRANSCRIPT_TABLES:
        problems += outdir.table_problems(folder / outdir.MERGED.format(what),
                                          outdir.TRANSCRIPT_KEYS, samples)  # fmt: skip
    try:
        outdir.tx2gene(folder)
    except QuantError as error:
        problems.append(str(error))
    return problems


def featured_contrast(spec: Dataset) -> tuple[Contrast | None, list[str]]:
    """The dataset's featured DE contrast, or None when it names none."""
    if spec.de is None:
        return None, []
    found = [c for c in spec.contrasts if c.name == spec.de.featured]
    if not found:
        return None, [f"{spec.path.name}: [de] featured {spec.de.featured!r} is no contrast"]
    return found[0], []


def _path(options: Options) -> list[str]:
    """One path, and only its own options: an outdir to read, or a genome to run Salmon on."""
    if (options.rnaseq_dir is None) == (options.genome is None):
        return ["give --rnaseq-dir (an nf-core/rnaseq outdir) or --genome (Salmon here on M3's "
                "trimmed reads), one of the two"]  # fmt: skip
    if options.genome is not None and options.salmon_dir is not None:
        return ["--salmon-dir reads copy Q1 of an nf-core outdir; the laptop path writes its own"]
    if options.genome is not None and not options.genome.is_file():
        return [f"--genome {options.genome} is not a file"]
    return []


def _required(*names: str) -> list[str]:
    """Each pinned tool that is missing or reports another version."""
    problems = []
    for name in names:
        try:
            tools.require(name)
        except tools.ToolError as error:
            problems.append(str(error))
    return problems


def resolve(settings: runtime.Settings, dataset: str, options: Options) -> Inputs:
    """Every input checked, every problem listed together before anything runs.

    Raises:
        QuantError: Listing every problem found.
    """
    root = settings.root
    spec, problems = _spec(options, dataset)
    sheet_path, runs, found = _sheet(root, dataset)
    problems += found + _path(options)
    gtf, found = annotation(root, spec) if spec else (None, [])
    problems += found
    sheet = runs.unique("sample", keep="first", maintain_order=True) if runs is not None else None
    samples = sheet["sample"].to_list() if sheet is not None else []
    counts, q1, trim = None, None, None
    reads: laptop.Reads = {}
    if options.rnaseq_dir is not None and sheet is not None:
        counts, found = _counts(options.rnaseq_dir, set(samples))
        q1 = salmon_dir(root, dataset, options.salmon_dir)
        problems += found + q1_problems(q1, samples)
    elif options.genome is not None and runs is not None:
        reads, trim, found = laptop.trimmed_reads(root, dataset, runs)
        problems += found + _required("salmon")
    contrast, found = featured_contrast(spec) if spec else (None, [])
    problems += found + _required("multiqc")
    if problems or spec is None or sheet is None or sheet_path is None or gtf is None:
        raise QuantError("ctx quant cannot start:\n  " + "\n  ".join(problems))
    return Inputs(spec, sheet_path, sheet, options.rnaseq_dir, gtf, gene_types(gtf), counts, q1,
                  contrast, reads, trim)  # fmt: skip


def run_salmon(inputs: Inputs, options: Options, out: Path,
               threads: int) -> tuple[Reference, Quantified]:  # fmt: skip
    """The laptop path: the probe, the reference, Salmon's index, then each sample.

    Raises:
        ToolError: If Salmon fails its probe or a run.
    """
    if options.genome is None:
        raise QuantError("the laptop path needs --genome")
    salmon.salmon_probe()
    ref = prepare(options.genome, inputs.gtf, out / "tools" / "index")
    built = salmon.index(ref, out / "tools" / "index", threads, decoys=options.salmon.decoys)
    return ref, salmon.quantify(inputs.reads, built, out / "tools" / inputs.spec.slug,
                                options.salmon, threads)  # fmt: skip


def _outdir(inputs: Inputs) -> Path:
    """The cohort path's nf-core/rnaseq outdir."""
    if inputs.rnaseq_dir is None:
        raise QuantError("the cohort path needs --rnaseq-dir")
    return inputs.rnaseq_dir


def measure(inputs: Inputs, options: Options, out: Path, threads: int) -> Measured:
    """Salmon where asked, MultiQC on its logs, the port on every quant.sf, fates and sums."""
    samples = inputs.samples["sample"].to_list()
    ref, quantified = (run_salmon(inputs, options, out, threads) if inputs.reads
                       else (None, None))  # fmt: skip
    logs = quantified.folder if quantified else _outdir(inputs) / "salmon"
    salmon_table = outdir.salmon_table(outdir.salmon_run(logs, out / "multiqc" / inputs.spec.slug,
                                                         samples))  # fmt: skip
    fates = tables.fates(salmon_table)
    folder = quantified.folder if quantified else inputs.salmon_dir
    port = (tximport.summarize(outdir.quants(folder, samples), outdir.tx2gene(folder))
            if folder is not None else None)  # fmt: skip
    if quantified is not None and port is not None:
        names = dict(outdir.tx2gene(quantified.folder).select("gene_id", "gene_name").iter_rows())
        counts = laptop.counts(port, names, releases_laptop(inputs, options, ref, quantified))
        scaled = port.counts_from(CountsFrom.LENGTH_SCALED_TPM).sum(axis=0)
        sums = dict(zip(port.samples, map(float, scaled), strict=True))
    else:
        if inputs.counts is None:
            raise QuantError("the cohort path read no length-scaled matrix")
        counts = inputs.counts
        sums = outdir.column_sums(logs, "gene_counts_length_scaled", samples)
    return Measured(salmon_table, fates, tables.sums(fates, sums), port, counts, folder,
                    quantified, ref)  # fmt: skip


def verify(inputs: Inputs, measured: Measured) -> None:
    """Every check of ``ctx.quant.checks``, before anything is written.

    The port is held to nf-core's merged tables on the cohort path only: the laptop path's
    tables are the port's own, and ``tests/test_quant_tximport.py`` holds the port to R.

    Raises:
        QuantError: Listing every failed check.
    """
    targets = tables.library(measured.salmon).select("sample", "targets")
    problems = checks.fate_problems(measured.fates) + checks.sum_problems(measured.sums, targets)
    problems += checks.ceiling_problems(measured.counts.matrix["gene_id"].to_list(),
                                        inputs.genes["gene_id"].to_list())  # fmt: skip
    if measured.port is not None and inputs.salmon_dir is not None:
        problems += checks.tximport_problems(measured.port, inputs.salmon_dir)
    checks.check(problems)


def releases_laptop(inputs: Inputs, options: Options, ref: Reference | None,
                    quantified: Quantified) -> dict[str, str]:  # fmt: skip
    """The laptop path's provenance: Salmon, the genome, the reads' and Salmon's folders."""
    genome = (f"{options.genome.name}, SHA-256 {ref.source.sha256}" if options.genome and ref
              else "none")  # fmt: skip
    return {"salmon": tools.TOOLS["salmon"].version, "genome": genome,
            "reads": inputs.trim.name if inputs.trim else "none",
            "salmon_index": quantified.index.name,
            "salmon_quant": quantified.folder.name}  # fmt: skip


def releases(inputs: Inputs, measured: Measured) -> dict[str, str]:
    """What every README records: the counts' provenance, the annotation and the inputs."""
    return {
        **measured.counts.provenance,
        "annotation": inputs.gtf.name, "tool:multiqc": tools.TOOLS["multiqc"].version,
        "input:samplesheet": inputs.sheet_path.name,
        "input:salmon_dir": str(measured.salmon_dir) if measured.salmon_dir else "none",
    }  # fmt: skip


def run(settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str,
        dataset: str, options: Options) -> list[Path]:  # fmt: skip
    """Read or run Salmon, check its outputs, and write the M4 page, tables and figures.

    Returns:
        The page, then every table and figure file.

    Raises:
        QuantError: Listing every input problem, when nothing has run; or every failed
            check, when only MultiQC's and Salmon's folders are written.
        ToolError: If Salmon fails its probe or a run.
    """
    inputs = resolve(settings, dataset, options)
    out, day = _out(settings.root), date.today()
    measured = measure(inputs, options, out, settings.threads)
    verify(inputs, measured)
    writer = Writer(out, dataset, day, releases(inputs, measured))
    page = pages.build(pages.Context(writer, inputs, options, measured, cards, command,
                                     settings))  # fmt: skip
    path = out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"
    return [write(page, path), *writer.files]
