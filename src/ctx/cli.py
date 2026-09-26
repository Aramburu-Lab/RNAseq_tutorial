"""Command-line entry point: one subcommand per teaching module."""

from typing import Annotated

import typer

from ctx import __version__

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
    version: Annotated[
        bool,
        typer.Option(
            "--version", callback=_print_version, is_eager=True, help="Print the version and exit."
        ),
    ] = False,
) -> None:
    """Cancer Transcriptomics Teaching Suite."""
