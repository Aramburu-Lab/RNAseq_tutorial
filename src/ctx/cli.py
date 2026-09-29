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
from ctx.report.builder import ParameterCard

app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="Cancer Transcriptomics Teaching Suite.",
)

# (meaning, what changes if you move it), per CLI option; the default is read from click. The
# meaning doubles as the --help text, so it describes the option on every command.
CARDS: dict[str, tuple[str, str]] = {
    "--profile": (
        "Where this runs: laptop (at most 4 threads) or hpc (inside a Slurm job). Unset: "
        "CTX_PROFILE, else laptop.",
        "Changes threads and memory, never results. hpc reads both from Slurm and needs --root.",
    ),
    "--root": (
        "The project directory every output is written under. Unset: CTX_ROOT, else the "
        "working directory on a laptop.",
        "Moves every output. Nothing is ever written inside the installed package.",
    ),
    "--dataset": (
        "Short name of the dataset, lowercase words joined by underscores; every output name "
        "starts with it.",
        "Renames every output and selects which dataset's earlier outputs this command reads.",
    ),
    "--fetchngs-dir": (
        "The nf-core/fetchngs outdir to read. Unset: <root>/data/raw/<dataset>.",
        "Points gather at another download. gather only reads it, never writes into it.",
    ),
    "--dataset-file": (
        "The dataset's TOML file, conf/datasets/<slug>.toml: its resources, sample curation, "
        "contrasts, gene sets and sanity rules.",
        "Changes what this command reads from the dataset; each page says what it uses.",
    ),
    "--min-count": (
        "Reads a gene needs at the median library size, in at least the smallest group's "
        "samples, to be tested (edgeR filterByExpr's min.count).",
        "Higher tests fewer, better-measured genes, so BH corrects for fewer tests, but weakly "
        "expressed biology drops out.",
    ),
    "--lfc-shrink": (
        "Also estimate each gene's log2 fold change with apeGLM shrinkage, which pulls it "
        "toward 0 as far as its counts leave it uncertain, so a weakly expressed gene cannot "
        "top a ranking on noise. --no-lfc-shrink turns it off.",
        "Off: the shrunk columns stay empty, the MA plot shows the raw estimate alone, and M5 "
        "skips the shrinkage step of each contrast. The test, p-values, calls, volcano and "
        "enrichment use the raw estimate either way, so none of them change.",
    ),
    "--rnaseq-dir": (
        "An nf-core/rnaseq outdir (Salmon only) whose gene matrix to read. Unset: no counts.",
        "Set: gather checks the release and the samples, rounds the length-scaled counts to "
        "integers and writes them. It only reads the outdir.",
    ),
}
# Where one option does different things on different commands, the effect is per command.
CARD_EFFECTS: dict[tuple[str, str], str] = {
    ("setup", "--dataset-file"): "Set: fetches each resource and keeps it only if its size "
    "and SHA-256 match. Unset: fetches nothing.",
    ("gather", "--dataset-file"): "Set: fills donor, cell type, condition, replicate and "
    "batch from the SDRF, and stops on any mismatch. Unset: they stay empty.",
    ("report", "--dataset-file"): "Set: the hub also decides and shows the file's sanity rules.",
    ("gather", "--dataset"): "Renames every output. Unless --fetchngs-dir is set, gather "
    "reads <root>/data/raw/<dataset>.",
    ("de", "--dataset-file"): "Selects the contrasts, their designs and expected n, the "
    "featured contrast, the lessons and the PCA subsets this page fits and draws.",
    ("function", "--dataset-file"): "Selects the gene sets, the featured set and the "
    "running score this page tests.",
}
_UNCARDED = {"--help", "--version"}


def parameter_cards(*commands: Any) -> tuple[ParameterCard, ...]:
    """One card per CLI option of the given commands, with the default read from click.

    Args:
        commands: click commands (typer bundles its own click, so typed loosely): the group
            and the subcommand whose options a page documents.

    Returns:
        The cards, in the order click lists the options.

    Raises:
        ValueError: If an option has no entry in ``CARDS``; a parameter without a card
            is a teaching page that is silently incomplete.
    """
    options = [p for c in commands for p in c.params if p.opts and p.opts[0] not in _UNCARDED]
    if missing := [p.opts[0] for p in options if p.opts[0] not in CARDS]:
        raise ValueError(f"no parameter card for {missing}")
    command = commands[-1].name if commands else ""
    cards = []
    for option in options:
        meaning, effect = CARDS[option.opts[0]]
        effect = CARD_EFFECTS.get((command, option.opts[0]), effect)
        if option.required:
            default = "required"
        elif option.secondary_opts:  # an on/off switch, such as --lfc-shrink/--no-lfc-shrink
            default = "on" if option.default else "off"
        else:
            default = "unset" if option.default is None else str(option.default)
        cards.append(ParameterCard(option.opts[0], meaning, default, effect))
    return tuple(cards)


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
    profile: Annotated[str | None, typer.Option(help=CARDS["--profile"][0])] = None,
    root: Annotated[Path | None, typer.Option(help=CARDS["--root"][0])] = None,
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
    dataset_file: Annotated[Path | None, typer.Option(help=CARDS["--dataset-file"][0])] = None,
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
def gather(
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=CARDS["--dataset"][0])],
    fetchngs_dir: Annotated[Path | None, typer.Option(help=CARDS["--fetchngs-dir"][0])] = None,
    dataset_file: Annotated[Path | None, typer.Option(help=CARDS["--dataset-file"][0])] = None,
    rnaseq_dir: Annotated[Path | None, typer.Option(help=CARDS["--rnaseq-dir"][0])] = None,
) -> None:
    """Check one nf-core/fetchngs outdir and write the harmonised sample sheet (M1c)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.gather import harmonise  # noqa: PLC0415 - heavy; imported after the runtime profile
    from ctx.gather.curation import CurationError  # noqa: PLC0415 - same reason
    from ctx.gather.fetchngs import FetchngsError  # noqa: PLC0415 - same module, same reason
    from ctx.gather.rnaseq import RnaseqError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    try:
        written = harmonise.run(
            ctx.obj["settings"], cards, command, dataset, fetchngs_dir=fetchngs_dir,
            dataset_file=dataset_file, rnaseq_dir=rnaseq_dir,
        )  # fmt: skip
    except (FetchngsError, DatasetError, CurationError, RnaseqError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def de(
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=CARDS["--dataset"][0])],
    dataset_file: Annotated[Path, typer.Option(help=CARDS["--dataset-file"][0])],
    min_count: Annotated[int, typer.Option(help=CARDS["--min-count"][0])] = 10,
    lfc_shrink: Annotated[
        bool, typer.Option("--lfc-shrink/--no-lfc-shrink", help=CARDS["--lfc-shrink"][0])
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
    dataset: Annotated[str, typer.Option(help=CARDS["--dataset"][0])],
    dataset_file: Annotated[Path, typer.Option(help=CARDS["--dataset-file"][0])],
) -> None:
    """Enrichment of every contrast (ORA and preranked GSEA); the M6 tables, figures and page."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.function import run as module  # noqa: PLC0415 - heavy; after the runtime profile

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    try:
        written = module.run(
            ctx.obj["settings"], cards, command, dataset, dataset_file=dataset_file
        )
    except DatasetError as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def report(
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=CARDS["--dataset"][0])],
    dataset_file: Annotated[Path | None, typer.Option(help=CARDS["--dataset-file"][0])] = None,
) -> None:
    """Write the hub linking every module's newest page, and decide the sanity rules."""
    from ctx.report import hub  # noqa: PLC0415 - imported after the runtime profile

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    for path in hub.run(ctx.obj["settings"], cards, dataset, dataset_file):
        typer.echo(path)
