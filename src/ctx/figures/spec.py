"""A figure's stable ID and its ``FigureSpec``: what it is, what drew it, what it read.

The ID is ``<cohort>_m<N>_<figure>``, N the module's number in the course (M1 Gather to M10
Clinical), so ``breast_m5_volcano`` names the volcano of every breast dataset; which dataset
and contrast it draws are its parameters. Output names still follow the naming lint.
"""

import re
from dataclasses import dataclass
from typing import Any

# The course's module numbers (TODO.md Module specification), not the output folders' NN.
NUMBERS = {
    "gather": 1, "qc": 2, "process": 3, "quant": 4, "de": 5, "function": 6, "proteo": 7,
    "ncrna": 8, "genecard": 9, "clinical": 10,
}  # fmt: skip
ID = re.compile(r"(?P<cohort>[a-z][a-z0-9]*)_m(?P<number>[1-9][0-9]?)_(?P<figure>[a-z0-9_]+)")


def figure_id(cohort: str, module: str, figure: str) -> str:
    """The stable ID of one figure kind of one cohort.

    Raises:
        ValueError: If the parts do not form an ID ``parse`` reads back.
    """
    built = f"{cohort}_m{NUMBERS[module]}_{figure}"
    if parse(built) != (cohort, module, figure):
        raise ValueError(f"not a figure ID: {built!r}")
    return built


def parse(figure: str) -> tuple[str, str, str] | None:
    """(cohort, module, figure kind) of an ID, or None for a string that is not one."""
    match = ID.fullmatch(figure)
    if match is None:
        return None
    module = next((m for m, n in NUMBERS.items() if n == int(match["number"])), None)
    return (match["cohort"], module, match["figure"]) if module else None


@dataclass(frozen=True)
class Input:
    """One file a figure reads, as the reader copies it into a clean folder.

    Attributes:
        name: The key in the params file's ``[inputs]``.
        file: The file's name in the reader's folder.
        source: Where the panel copies it from: under ``$ROOT`` or ``$CTX_REPO``, or a path.
        base: ``root``, ``repo`` or ``path``: what ``source`` is relative to.
        sha256: The SHA-256 of its bytes when the figure was drawn.
    """

    name: str
    file: str
    source: str
    base: str
    sha256: str


@dataclass(frozen=True)
class FigureSpec:
    """Everything the Reproduce panel, the params file and the README say about one figure.

    Attributes:
        figure_id: ``<cohort>_m<N>_<figure>``.
        dataset: The dataset slug in its output names.
        module: The naming token of the module that drew it.
        params: The figure's parameter dataclass, filled.
        inputs: Every file it reads, in the params file's order.
        plotting: The function that draws it, ``module.path:function``.
        plotting_file: That function's file in the repository.
        stem: The figure's file stem; ``ctx figure`` writes the same name, dated the day it runs.
        data_sha256: The SHA-256 of its data CSV, which a redraw must equal.
        root: The ctx root the inputs sit in, as a page shows a path.
        commit: Where the drawing ctx came from (``runtime.source``).
        version: The drawing ctx's version.
        upstream: The command that wrote the inputs and drew the figure.
    """

    figure_id: str
    dataset: str
    module: str
    params: Any
    inputs: tuple[Input, ...]
    plotting: str
    plotting_file: str
    stem: str
    data_sha256: str
    root: str
    commit: str
    version: str
    upstream: str
