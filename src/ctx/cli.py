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
from ctx.report import text
from ctx.report.builder import ParameterCard

app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="Cancer Transcriptomics Teaching Suite.",
)

# (meaning, what changes if you move it), per CLI option, from content/cards.toml; the default
# is read from click. The meaning doubles as the --help text, so it describes the option on
# every command.
CARDS: dict[str, tuple[str, str]] = {
    option: (card["meaning"], card["effect"])
    for option, card in text.raw("cards", "options").items()
}
# Where one option does different things on different commands, the effect is per command.
CARD_EFFECTS: dict[tuple[str, str], str] = {
    (command, option): effect
    for command, effects in text.raw("cards", "effects").items()
    for option, effect in effects.items()
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


def _columns(value: str) -> tuple[str, ...]:
    """A comma-separated option as its column names."""
    return tuple(c.strip() for c in value.split(",") if c.strip())


@app.command()
@profiled
def qc(  # noqa: PLR0913, PLR0917 - typer reads one parameter per CLI option
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=CARDS["--dataset"][0])],
    dataset_file: Annotated[Path, typer.Option(help=CARDS["--dataset-file"][0])],
    registry: Annotated[Path | None, typer.Option(help=CARDS["--registry"][0])] = None,
    strata: Annotated[str, typer.Option(help=CARDS["--strata"][0])] = "cell_type",
    peers: Annotated[str, typer.Option(help=CARDS["--peers"][0])] = "cell_type,condition",
    mad: Annotated[float | None, typer.Option(help=CARDS["--mad"][0])] = None,
) -> None:
    """Judge every sample by the QC registry and test the flags for confounding (M2)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.qc import QcError  # noqa: PLC0415 - same reason
    from ctx.qc import run as module  # noqa: PLC0415 - heavy; imported after the runtime profile
    from ctx.qc.registry import RegistryError  # noqa: PLC0415 - same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = runtime.retypable(sys.argv[1:])
    options = module.Options(dataset_file, registry, _columns(strata), _columns(peers), mad)
    try:
        written = module.run(ctx.obj["settings"], cards, command, dataset, options)
    except (DatasetError, QcError, RegistryError) as error:
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
