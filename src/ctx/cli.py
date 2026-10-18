"""Command-line entry point: one subcommand per teaching module.

Nothing heavy is imported here. Each command resolves its runtime profile first (``profiled``),
and only then imports what it needs, so thread limits hold for every library.
"""

import functools
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Annotated, Any

import typer

from ctx import __version__, runtime
from ctx.gather import own_fastq
from ctx.genecard.hpa import RELEASE as HPA_RELEASE
from ctx.process import DEFAULT_SETTINGS, Aligner, Strandedness
from ctx.qc import Fastp
from ctx.quant import CountsFrom, SalmonOptions
from ctx.report.cards import HELP, parameter_cards

app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="Cancer Transcriptomics Teaching Suite.",
)


def _print_version(value: bool) -> None:
    """Print the installed version and stop, when --version was given.

    Args:
        value: True when the flag was passed.

    Raises:
        typer.Exit: After printing, so no subcommand runs.
    """
    if value:
        typer.echo(__version__)
        raise typer.Exit


@app.callback()
def main(
    ctx: typer.Context,
    profile: Annotated[str | None, typer.Option(help=HELP["--profile"])] = None,
    root: Annotated[Path | None, typer.Option(help=HELP["--root"])] = None,
    _version: Annotated[
        bool,
        typer.Option(
            "--version", callback=_print_version, is_eager=True, help="Print the version and exit."
        ),
    ] = False,
) -> None:
    """Cancer Transcriptomics Teaching Suite."""
    # Only stored here: click runs this before a subcommand's --help, so resolving the
    # profile now would make `ctx setup --help` fail on a login node.
    ctx.obj = {"profile": profile, "root": root}


def profiled(command: Callable[..., Any]) -> Callable[..., Any]:
    """Resolve the profile and root before the command body, and so before any heavy import.

    Args:
        command: A subcommand taking ``ctx: typer.Context``.

    Returns:
        The wrapped command; it finds the settings in ``ctx.obj["settings"]``.
    """

    @functools.wraps(command)
    def wrapper(ctx: typer.Context, *args: Any, **kwargs: Any) -> Any:
        """Store the resolved settings in the context, then run the command."""
        options = ctx.find_root().obj or {}
        ctx.obj = {
            **options,
            "settings": runtime.start(options.get("profile"), options.get("root")),
        }
        return command(ctx, *args, **kwargs)

    return wrapper


@app.command()
@profiled
def setup(
    ctx: typer.Context,
    dataset_file: Annotated[Path | None, typer.Option(help=HELP["--dataset-file"])] = None,
) -> None:
    """Report this machine: versions, hardware, profile, fonts and colour contrast."""
    from ctx import setup as module  # noqa: PLC0415 - heavy; imported after the runtime profile
    from ctx.dataset import DatasetError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset_file)
    except (module.FetchError, DatasetError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def gather(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    fetchngs_dir: Annotated[Path | None, typer.Option(help=HELP["--fetchngs-dir"])] = None,
    dataset_file: Annotated[Path | None, typer.Option(help=HELP["--dataset-file"])] = None,
    rnaseq_dir: Annotated[Path | None, typer.Option(help=HELP["--rnaseq-dir"])] = None,
    gdc_dir: Annotated[list[Path] | None, typer.Option(help=HELP["--gdc-dir"])] = None,
    mirna_dir: Annotated[Path | None, typer.Option(help=HELP["--mirna-dir"])] = None,
    isoform_dir: Annotated[Path | None, typer.Option(help=HELP["--isoform-dir"])] = None,
    counts: Annotated[Path | None, typer.Option(help=HELP["--counts"])] = None,
) -> None:
    """Check a fetchngs outdir, GDC folders or a count matrix; write the sample sheet (M1c)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.gather import counts_run, gdc_run, harmonise  # noqa: PLC0415 - after the profile
    from ctx.gather.counts import CountsError  # noqa: PLC0415 - same reason
    from ctx.gather.curation import CurationError  # noqa: PLC0415 - same reason
    from ctx.gather.fetchngs import FetchngsError  # noqa: PLC0415 - same module, same reason
    from ctx.gather.gdc import GdcError  # noqa: PLC0415 - same reason
    from ctx.gather.rnaseq import RnaseqError  # noqa: PLC0415 - same reason
    from ctx.gather.trace import TraceError  # noqa: PLC0415 - same reason
    from ctx.qc.registry import RegistryError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    try:
        if counts is not None:
            if fetchngs_dir or rnaseq_dir or gdc_dir or mirna_dir or isoform_dir:
                raise DatasetError("--counts is a route of its own: no other input folder")
            if dataset_file is None:
                raise DatasetError("--counts takes --dataset-file: the sample table and [counts]")
            written = counts_run.run(ctx.obj["settings"], cards, command, dataset,
                                     matrix=counts, dataset_file=dataset_file)  # fmt: skip
        elif gdc_dir:
            if fetchngs_dir or rnaseq_dir or dataset_file is None:
                raise DatasetError("--gdc-dir takes --dataset-file, and neither --fetchngs-dir "
                                   "nor --rnaseq-dir")  # fmt: skip
            written = gdc_run.run(ctx.obj["settings"], cards, command, dataset, gdc_dirs=gdc_dir,
                                  dataset_file=dataset_file, mirna_dir=mirna_dir,
                                  isoform_dir=isoform_dir)  # fmt: skip
        elif mirna_dir or isoform_dir:
            raise DatasetError("--mirna-dir and --isoform-dir take --gdc-dir")
        else:
            written = harmonise.run(
                ctx.obj["settings"], cards, command, dataset, fetchngs_dir=fetchngs_dir,
                dataset_file=dataset_file, rnaseq_dir=rnaseq_dir,
            )  # fmt: skip
    except (FetchngsError, DatasetError, CurationError, RnaseqError, RegistryError,
            GdcError, TraceError, CountsError) as error:  # fmt: skip
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
def reads(
    run_folder: Annotated[Path, typer.Option(help=HELP["--run-folder"])],
    fastq_dir: Annotated[Path, typer.Option(help=HELP["--fastq-dir"])],
    table: Annotated[Path, typer.Option(help=HELP["--table"])],
    sample_column: Annotated[
        str, typer.Option(help=HELP["--sample-column"])
    ] = own_fastq.DEFAULT_COLUMN,
) -> None:
    """Lay out your own FASTQ files as a fetchngs download, for nf-core/rnaseq and gather."""
    # Standard library only and no page, so no profile to resolve first; on Kebnekaise the
    # same function runs as a file (ctx.gather.own_fastq's Usage).
    try:
        written = own_fastq.write(run_folder, fastq_dir, table, sample_column)
    except own_fastq.ReadsError as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)
    typer.echo(own_fastq.summary(written[0]))


def _columns(value: str) -> tuple[str, ...]:
    """A comma-separated option as its column names."""
    return tuple(c.strip() for c in value.split(",") if c.strip())


@app.command()
@profiled
def qc(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    registry: Annotated[Path | None, typer.Option(help=HELP["--registry"])] = None,
    strata: Annotated[str, typer.Option(help=HELP["--strata"])] = "cell_type",
    peers: Annotated[str, typer.Option(help=HELP["--peers"])] = "cell_type,condition",
    mad: Annotated[float | None, typer.Option(help=HELP["--mad"])] = None,
    rnaseq_dir: Annotated[Path | None, typer.Option(help=HELP["--rnaseq-dir"])] = None,
    reads: Annotated[bool, typer.Option("--reads/--no-reads", help=HELP["--reads"])] = False,
    fetchngs_dir: Annotated[Path | None, typer.Option(help=HELP["--fetchngs-dir"])] = None,
    fastp_q: Annotated[int, typer.Option(help=HELP["--fastp-q"])] = Fastp.q,
    fastp_u: Annotated[int, typer.Option(help=HELP["--fastp-u"])] = Fastp.u,
    fastp_l: Annotated[int, typer.Option(help=HELP["--fastp-l"])] = Fastp.l,
    cut_tail: Annotated[
        bool, typer.Option("--cut-tail/--no-cut-tail", help=HELP["--cut-tail"])
    ] = Fastp.cut_tail,
    cut_window: Annotated[int, typer.Option(help=HELP["--cut-window"])] = Fastp.window,
    cut_mean: Annotated[int, typer.Option(help=HELP["--cut-mean"])] = Fastp.mean,
    grid_run: Annotated[str | None, typer.Option(help=HELP["--grid-run"])] = None,
) -> None:
    """Judge every sample by the QC registry and test the flags for confounding (M2)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.qc import QcError  # noqa: PLC0415 - same reason
    from ctx.qc import run as module  # noqa: PLC0415 - heavy; imported after the runtime profile
    from ctx.qc.registry import RegistryError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    fastp = Fastp(fastp_q, fastp_u, fastp_l, cut_tail, cut_window, cut_mean)
    options = module.Options(dataset_file, registry, _columns(strata), _columns(peers), mad,
                             rnaseq_dir, reads, fetchngs_dir, fastp, grid_run)  # fmt: skip
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset, options)
    except (DatasetError, QcError, RegistryError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def process(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    genome: Annotated[Path, typer.Option(help=HELP["--genome"])],
    aligner: Annotated[Aligner, typer.Option(help=HELP["--aligner"])] = Aligner.STAR,
    strandedness: Annotated[
        Strandedness, typer.Option(help=HELP["--strandedness"])
    ] = Strandedness.AUTO,
    multimap_max: Annotated[int | None, typer.Option(help=HELP["--multimap-max"])] = None,
    settings: Annotated[str, typer.Option(help=HELP["--settings"])] = DEFAULT_SETTINGS,
    featured_sample: Annotated[str | None, typer.Option(help=HELP["--featured-sample"])] = None,
    reads_per_run: Annotated[int | None, typer.Option(help=HELP["--reads-per-run"])] = None,
    compare_dataset: Annotated[str | None, typer.Option(help=HELP["--compare-dataset"])] = None,
    fetchngs_dir: Annotated[Path | None, typer.Option(help=HELP["--fetchngs-dir"])] = None,
    registry: Annotated[Path | None, typer.Option(help=HELP["--registry"])] = None,
    fastp_q: Annotated[int, typer.Option(help=HELP["--fastp-q"])] = Fastp.q,
    fastp_u: Annotated[int, typer.Option(help=HELP["--fastp-u"])] = Fastp.u,
    fastp_l: Annotated[int, typer.Option(help=HELP["--fastp-l"])] = Fastp.l,
    cut_tail: Annotated[
        bool, typer.Option("--cut-tail/--no-cut-tail", help=HELP["--cut-tail"])
    ] = Fastp.cut_tail,
    cut_window: Annotated[int, typer.Option(help=HELP["--cut-window"])] = Fastp.window,
    cut_mean: Annotated[int, typer.Option(help=HELP["--cut-mean"])] = Fastp.mean,
) -> None:
    """Trim, align and count every sample's reads, and test the settings that move counts (M3)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.process import ProcessError  # noqa: PLC0415 - same reason
    from ctx.process import run as module  # noqa: PLC0415 - heavy; after the runtime profile
    from ctx.qc import QcError  # noqa: PLC0415 - same reason
    from ctx.qc.registry import RegistryError  # noqa: PLC0415 - same reason
    from ctx.tools import ToolError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    fastp = Fastp(fastp_q, fastp_u, fastp_l, cut_tail, cut_window, cut_mean)
    options = module.Options(
        dataset_file, genome, aligner, strandedness, multimap_max, _columns(settings),
        featured_sample, reads_per_run, compare_dataset, fetchngs_dir, registry, fastp,
    )  # fmt: skip
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset, options)
    except (DatasetError, QcError, ProcessError, RegistryError, ToolError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def quant(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    rnaseq_dir: Annotated[Path | None, typer.Option(help=HELP["--rnaseq-dir"])] = None,
    salmon_dir: Annotated[Path | None, typer.Option(help=HELP["--salmon-dir"])] = None,
    counts_from_abundance: Annotated[
        CountsFrom, typer.Option(help=HELP["--counts-from-abundance"])
    ] = CountsFrom.LENGTH_SCALED_TPM,
    genome: Annotated[Path | None, typer.Option(help=HELP["--genome"])] = None,
    libtype: Annotated[str, typer.Option(help=HELP["--libtype"])] = SalmonOptions.libtype,
    gc_bias: Annotated[
        bool, typer.Option("--gc-bias/--no-gc-bias", help=HELP["--gc-bias"])
    ] = False,
    seq_bias: Annotated[
        bool, typer.Option("--seq-bias/--no-seq-bias", help=HELP["--seq-bias"])
    ] = False,
    decoys: Annotated[bool, typer.Option("--decoys/--no-decoys", help=HELP["--decoys"])] = True,
) -> None:
    """Read or run Salmon on every sample, and see how transcripts become gene counts (M4)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.quant import QuantError  # noqa: PLC0415 - same reason
    from ctx.quant import run as module  # noqa: PLC0415 - heavy; after the runtime profile
    from ctx.tools import ToolError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    options = module.Options(dataset_file, rnaseq_dir, salmon_dir, counts_from_abundance, genome,
                             SalmonOptions(libtype, gc_bias, seq_bias, decoys))  # fmt: skip
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset, options)
    except (DatasetError, QuantError, ToolError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def de(
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    min_count: Annotated[int, typer.Option(help=HELP["--min-count"])] = 10,
    lfc_shrink: Annotated[
        bool, typer.Option("--lfc-shrink/--no-lfc-shrink", help=HELP["--lfc-shrink"])
    ] = True,
) -> None:
    """Fit the dataset's contrasts with DESeq2 and write the M5 tables, figures and page."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.de import run as module  # noqa: PLC0415 - heavy; imported after the runtime profile
    from ctx.de.deseq import ContrastError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    try:
        written = module.run(
            ctx.obj["settings"], cards, command, dataset, dataset_file=dataset_file,
            min_count=min_count, lfc_shrink=lfc_shrink,
        )  # fmt: skip
    except (DatasetError, ContrastError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def function(
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    genes_from: Annotated[Path | None, typer.Option(help=HELP["--genes-from"])] = None,
) -> None:
    """Enrichment of every contrast (ORA and preranked GSEA); the M6 tables, figures and page."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.function import run as module  # noqa: PLC0415 - heavy; after the runtime profile
    from ctx.report.genes import GeneListError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    try:
        written = module.run(
            ctx.obj["settings"], cards, command, dataset, dataset_file=dataset_file,
            genes_from=genes_from,
        )  # fmt: skip
    except (DatasetError, GeneListError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def ncrna(
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    biotypes: Annotated[str, typer.Option(help=HELP["--biotypes"])] = "protein_coding,lncRNA",
    min_count: Annotated[int, typer.Option(help=HELP["--min-count"])] = 10,
) -> None:
    """Non-coding genes by biotype, and each lncRNA beside its nearest gene; the M8 page."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.ncrna import NcrnaError  # noqa: PLC0415 - same reason
    from ctx.ncrna import run as module  # noqa: PLC0415 - heavy; after the runtime profile

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    options = module.Options(dataset_file, _columns(biotypes), min_count)
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset, options)
    except (DatasetError, NcrnaError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def clinical(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    endpoint: Annotated[str | None, typer.Option(help=HELP["--endpoint"])] = None,
    split: Annotated[str, typer.Option(help=HELP["--split"])] = "median",
    horizon: Annotated[float, typer.Option(help=HELP["--horizon"])] = 60.0,
    covariates: Annotated[str, typer.Option(help=HELP["--covariates"])] = "age,stage",
    cox_strata: Annotated[str, typer.Option(help=HELP["--cox-strata"])] = "",
    genes_from: Annotated[Path | None, typer.Option(help=HELP["--genes-from"])] = None,
    min_log2cpm: Annotated[float, typer.Option(help=HELP["--min-log2cpm"])] = 1.0,
) -> None:
    """Relate expression to TCGA-CDR endpoints and clinical variables; the M10 page."""
    from ctx.clinical import ClinicalError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.clinical import run as module  # noqa: PLC0415 - heavy; after the runtime profile
    from ctx.dataset import DatasetError  # noqa: PLC0415 - same reason
    from ctx.report.genes import GeneListError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    options = module.Options(dataset_file, endpoint, split, horizon, _columns(covariates),
                             _columns(cox_strata), genes_from, min_log2cpm)  # fmt: skip
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset, options)
    except (DatasetError, ClinicalError, GeneListError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def proteo(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    method: Annotated[str, typer.Option(help=HELP["--method"])] = "spearman",
    min_quantified: Annotated[float, typer.Option(help=HELP["--min-quantified"])] = 0.7,
    lfc_min: Annotated[float, typer.Option(help=HELP["--lfc-min"])] = 1.0,
    gene_class: Annotated[str, typer.Option(help=HELP["--gene-class"])] = "all",
) -> None:
    """Compare RNA with protein per gene, or as fold changes across studies (M7)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.proteo import ProteoError  # noqa: PLC0415 - same reason
    from ctx.proteo import run as module  # noqa: PLC0415 - heavy; after the runtime profile

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    options = module.Options(dataset_file, method, min_quantified, lfc_min, gene_class)
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset, options)
    except (DatasetError, ProteoError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def genecard(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    genes: Annotated[str, typer.Option(help=HELP["--genes"])] = "",
    genes_from: Annotated[Path | None, typer.Option(help=HELP["--genes-from"])] = None,
    contrast: Annotated[str | None, typer.Option(help=HELP["--contrast"])] = None,
    hpa_release: Annotated[str, typer.Option(help=HELP["--hpa-release"])] = HPA_RELEASE,
    offline: Annotated[bool, typer.Option("--offline/--online", help=HELP["--offline"])] = False,
    refresh: Annotated[
        bool, typer.Option("--refresh/--no-refresh", help=HELP["--refresh"])
    ] = False,
) -> None:
    """Gene cards: this dataset's results beside GTEx, HPA and Ensembl data (M9)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.genecard import GenecardError  # noqa: PLC0415 - same reason
    from ctx.genecard import run as module  # noqa: PLC0415 - heavy; after the runtime profile
    from ctx.report.genes import GeneListError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    options = module.Options(dataset_file, _columns(genes), genes_from, contrast, hpa_release,
                             offline, refresh)  # fmt: skip
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset, options)
    except (DatasetError, GenecardError, GeneListError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def report(
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path | None, typer.Option(help=HELP["--dataset-file"])] = None,
) -> None:
    """Write the hub linking every module's newest page, and decide the sanity rules."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.report import hub  # noqa: PLC0415 - imported after the runtime profile

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    try:
        written = hub.run(ctx.obj["settings"], cards, dataset, dataset_file)
    except DatasetError as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def run(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=HELP["--dataset"])],
    dataset_file: Annotated[Path, typer.Option(help=HELP["--dataset-file"])],
    fetchngs_dir: Annotated[Path | None, typer.Option(help=HELP["--fetchngs-dir"])] = None,
    rnaseq_dir: Annotated[Path | None, typer.Option(help=HELP["--rnaseq-dir"])] = None,
    gdc_dir: Annotated[list[Path] | None, typer.Option(help=HELP["--gdc-dir"])] = None,
    salmon_dir: Annotated[Path | None, typer.Option(help=HELP["--salmon-dir"])] = None,
    mirna_dir: Annotated[Path | None, typer.Option(help=HELP["--mirna-dir"])] = None,
    isoform_dir: Annotated[Path | None, typer.Option(help=HELP["--isoform-dir"])] = None,
    counts: Annotated[Path | None, typer.Option(help=HELP["--counts"])] = None,
    check: Annotated[bool, typer.Option("--check/--no-check", help=HELP["--check"])] = False,
) -> None:
    """Run every module the dataset file declares, in course order, each as its own command."""
    from ctx import pipeline  # noqa: PLC0415 - imported after the runtime profile

    settings = ctx.obj["settings"]
    gdc_dirs = tuple(gdc_dir or ())
    options = pipeline.Options(dataset, dataset_file, fetchngs_dir, rnaseq_dir, gdc_dirs,
                               salmon_dir, mirna_dir, isoform_dir, counts)  # fmt: skip
    problems, planned = pipeline.check(settings.root, options)
    for problem in problems:
        typer.echo(f"- {problem.text}", err=True)
        if problem.lines:
            typer.echo(f"  add these lines to the dataset file:\n{problem.lines}", err=True)
    stops = [p for p in problems if p.blocks]
    typer.echo(f"{len(problems)} problems, {len(stops)} of them stopping the run; ctx run runs: "
               f"{', '.join(planned) or 'nothing'}")  # fmt: skip
    if check or stops:
        raise typer.Exit(1 if (problems if check else stops) else 0)
    profile = ctx.find_root().obj.get("profile")
    raise typer.Exit(pipeline.execute(settings.root, profile, pipeline.planned(options)))


@app.command()
@profiled
def guide(
    ctx: typer.Context,
    dataset_file: Annotated[Path | None, typer.Option(help=HELP["--dataset-file"])] = None,
) -> None:
    """Write the custom-data guide: download, count and run a study, from its dataset file."""
    from ctx.dataset import DatasetError, load  # noqa: PLC0415 - imported after the runtime profile
    from ctx.report import guide as module  # noqa: PLC0415 - same reason

    group = ctx.find_root().command
    # The page teaches `ctx run`, so its cards are run's options and guide's own.
    shared = parameter_cards(group, group.commands["run"], ctx.command)  # type: ignore[attr-defined]
    cards = tuple({card.name: card for card in shared}.values())
    try:
        spec = load(dataset_file or module.DEFAULT)
    except (DatasetError, OSError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    typer.echo(module.run(ctx.obj["settings"], cards, runtime.retypable(sys.argv[1:]), spec))
