"""Command-line entry point: one subcommand per teaching module.

Nothing heavy is imported here. Each command resolves its runtime profile first (``profiled``),
and only then imports what it needs, so thread limits hold for every library.
"""

import functools
import shlex
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

# (meaning, what changes if you move it), per CLI option; the default is read from click.
CARDS: dict[str, tuple[str, str]] = {
    "--profile": (
        "Where this runs: laptop, or hpc inside a Slurm job. Unset means CTX_PROFILE, else laptop.",
        "hpc takes threads and memory from the Slurm allocation and needs an explicit --root; "
        "laptop uses at most 4 threads. The results do not change, only the resources.",
    ),
    "--root": (
        "The project directory outputs are written under. Unset means CTX_ROOT, else the "
        "working directory (laptop only).",
        "Moves every output. Nothing is ever written inside the installed package.",
    ),
    "--dataset": (
        "Short name for this download: lowercase words joined by underscores. It names the "
        "outputs, and the default fetchngs outdir <root>/data/raw/<dataset>.",
        "Changes the output names and where gather looks when --fetchngs-dir is unset.",
    ),
    "--fetchngs-dir": (
        "The nf-core/fetchngs outdir to read. Unset means <root>/data/raw/<dataset>.",
        "Points gather at another download. gather only reads it and never writes into it.",
    ),
    "--dataset-file": (
        "The cohort's dataset file (conf/datasets/<slug>.toml): its SDRF join, value maps, "
        "batch rule and expected counts. Unset leaves the curation columns empty.",
        "Set, gather joins the SDRF (fetched by ctx setup) to every run, fills donor, cell "
        "type, condition, replicate and batch, and stops with a problem table if anything "
        "fails to join, map or match the declared counts.",
    ),
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
    cards = []
    for option in options:
        meaning, effect = CARDS[option.opts[0]]
        default = (
            "required"
            if option.required
            else ("unset" if option.default is None else str(option.default))
        )
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
def setup(ctx: typer.Context) -> None:
    """Report this machine: versions, hardware, profile, fonts and colour contrast."""
    from ctx import setup as module  # noqa: PLC0415 - heavy; imported after the runtime profile

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    written = module.run(ctx.obj["settings"], cards, "ctx " + shlex.join(sys.argv[1:]))
    for path in written:
        typer.echo(path)


@app.command()
@profiled
def gather(
    ctx: typer.Context,
    dataset: Annotated[str, typer.Option(help=CARDS["--dataset"][0])],
    fetchngs_dir: Annotated[Path | None, typer.Option(help=CARDS["--fetchngs-dir"][0])] = None,
    dataset_file: Annotated[Path | None, typer.Option(help=CARDS["--dataset-file"][0])] = None,
) -> None:
    """Check one nf-core/fetchngs outdir and write the harmonised sample sheet (M1c)."""
    from ctx.dataset import DatasetError  # noqa: PLC0415 - imported after the runtime profile
    from ctx.gather import harmonise  # noqa: PLC0415 - heavy; imported after the runtime profile
    from ctx.gather.curation import CurationError  # noqa: PLC0415 - same reason
    from ctx.gather.fetchngs import FetchngsError  # noqa: PLC0415 - same module, same reason

    cards = parameter_cards(ctx.find_root().command, ctx.command)
    command = "ctx " + shlex.join(sys.argv[1:])
    try:
        written = harmonise.run(
            ctx.obj["settings"], cards, command, dataset, fetchngs_dir, dataset_file
        )
    except (FetchngsError, DatasetError, CurationError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)
