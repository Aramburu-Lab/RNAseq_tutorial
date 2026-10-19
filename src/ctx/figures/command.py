"""``ctx figure <id> --params <file>``: redraw one figure from its inputs alone.

It writes into the reader's own folder, never into a ctx root, where the redrawn files would
take a module's output names (CONVENTIONS.md Outputs). The files are the ones the module
writes for that figure, PDF and SVG per theme, its page and its data, named by the same rule
and dated the day it runs; the page links no glossary and carries no Reproduce panel, since
none sits beside it (``ctx.report.glossary.page_path``).
"""

import importlib
from datetime import date
from pathlib import Path

import typer

from ctx.figures.kind import NoFigureError
from ctx.figures.params import ParamsError, Request, load_params
from ctx.figures.record import FOLDERS, sha256
from ctx.io.naming import stem


class RootError(ValueError):
    """An output folder inside a ctx root."""


def inside_root(folder: Path) -> Path | None:
    """The ctx root holding ``folder``, or None: a module folder, or a folder with ``outputs/``.

    A root is recognised by its module folders, as ``glossary.page_path`` recognises a page's.
    """
    for place in (folder.resolve(), *folder.resolve().parents):
        if place.name in FOLDERS.values():
            return place.parent.parent
        outputs = place / "outputs"
        if outputs.is_dir() and any((outputs / f).is_dir() for f in FOLDERS.values()):
            return place
    return None


def draw(request: Request, out: Path, day: date) -> list[Path]:
    """Draw the requested figure into ``out`` with the module's own figure function.

    Returns:
        The files written, as ``ctx.viz.save.write_figure`` lists them.

    Raises:
        RootError: If ``out`` is inside a ctx root; nothing is written.
    """
    from ctx.viz.save import write_figure  # noqa: PLC0415 - heavy; after the runtime profile

    if (root := inside_root(out)) is not None:
        raise RootError(f"{out}: inside the ctx root {root}; ctx figure writes outside any root")
    module, function = request.kind.draw.split(":")
    data, words = getattr(importlib.import_module(module), function)(
        request.inputs, request.params, request.spec
    )
    what = request.kind.what(request.params)
    figure_stem = stem(request.kind.module, f"{request.spec.slug}_{what}", 1, day=day)
    releases = {f"sha256:{p.name}": sha256(p) for p in request.inputs.values()}
    return write_figure(data, out, figure_stem, title=words.title, data_releases=releases,
                        text=words)  # fmt: skip


def main(figure: str, params: Path, out: Path) -> None:
    """The command: check the params file, draw, and print each file written.

    Raises:
        typer.Exit: With code 1 after printing the problems, when nothing could be drawn.
    """
    try:
        written = draw(load_params(params, figure), out, date.today())
    except (ParamsError, RootError, NoFigureError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(path)
