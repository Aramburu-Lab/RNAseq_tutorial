"""A module's writer that writes nothing: it keeps each figure it is handed, by description.

A drawer that calls a module's own block function, the one that builds a figure and its words
and hands both to ``Writer.figure``, passes this instead: the figure and its words are the
module's, built by its one definition, and nothing lands in a module folder. A table the
block writes for its page's views is dropped; the redraw writes the figure's files alone.
"""

from collections.abc import Mapping
from datetime import date
from pathlib import Path
from typing import Any

import polars as pl

from ctx.de.run import Writer
from ctx.figures.kind import NoFigureError
from ctx.report.builder import FigureBlock, FigureText
from ctx.viz.save import Figure, figure_block

# Any day: a stem built here names no file, only the block the module's code reads back.
DAY = date(2000, 1, 1)


class Recorder(Writer):
    """``Writer`` for one module folder, keeping figures instead of writing them."""

    def __init__(self, folder: str, dataset: str) -> None:
        """A writer named for ``folder`` (``02_qc``), so its stems carry the module's token."""
        super().__init__(Path(folder), dataset, DAY, {})
        self.drawn: dict[str, tuple[Figure, FigureText]] = {}

    def table(self, frame: pl.DataFrame, what: str, columns: Mapping[str, tuple[str, str]],
              csv_decimals: int | None = None) -> None:  # fmt: skip
        """Drop a table: a redraw writes only the figure's own files."""

    def figure(self, data: Figure, what: str, text: FigureText) -> FigureBlock:
        """Keep the figure and its words, and return its block as the module's code expects."""
        data = self.marked(data)
        self.drawn[what] = (data, text)
        return figure_block(data, self.stem(what), text)

    def one(self, what: str) -> tuple[Any, FigureText]:
        """The figure the module drew as ``<dataset>_<what>``.

        Raises:
            NoFigureError: When the module's code drew no such figure from these inputs.
        """
        if what not in self.drawn:
            raise NoFigureError(f"{what}: the module draws no such figure from these inputs")
        return self.drawn[what]
