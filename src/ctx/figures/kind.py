"""What a figure kind is: its parameters, its inputs and its drawer (``ctx.figures.kinds``).

Each module's kinds sit in ``kinds_m<N>``; this file holds what they share. An input is a
module's table by ``(module, what after the dataset slug)``, or one of the sources below.
"""

import re
from collections.abc import Callable, Collection
from dataclasses import dataclass, field, fields
from typing import Any

from ctx.dataset import Dataset

# The dataset file the module read, and the QC registry gather's depth standard comes from.
DATASET_FILE = ("dataset_file", "")
QC_REGISTRY = ("qc_registry", "")
# The registry ``ctx qc`` read, which its figure options table names (``ctx.qc.figure_tables``).
QC_REGISTRY_RUN = ("qc_registry_run", "")
# A file ``ctx setup`` fetched, by its resource name in the dataset file: ("resource", name).
RESOURCE = "resource"

Inputs = dict[str, tuple[str, str]]


class NoFigureError(ValueError):
    """Parameters under which the module draws no figure, saying why."""


def resource(name: str) -> tuple[str, str]:
    """The input that is the dataset file's resource ``name``, under the root it was fetched to."""
    return (RESOURCE, name)


def _no_check(_: Any, __: Dataset) -> list[str]:
    """A kind whose parameters the dataset file bounds no further."""
    return []


@dataclass(frozen=True)
class Kind:
    """One sort of figure: its module, name, parameters, inputs and drawer.

    Attributes:
        module: The naming token of the module that writes it.
        figure: The kind's name in the ID.
        params: Its parameter dataclass.
        plotting: The module's own function that builds the figure, ``module:function``.
        draw: The drawer ``ctx figure`` calls, ``module:function``.
        made: Every figure of this kind the dataset file makes, as filled parameters.
        what: A figure's output description after the dataset slug.
        inputs: Its inputs by name, each ``(module, what)`` of a table or a special source.
        check: Problems with parameter values the dataset file does not declare.
        unless: Tables, as inputs name them, whose presence for the dataset means the module
            wrote another figure under this kind's stem: such a figure has no spec.
        parsed: For a kind whose figures the dataset file cannot list (a gene the reader
            named on the command line): the parameters of a description, or None for one
            that is not this kind's. ``what`` of the result must give the description back.
    """

    module: str
    figure: str
    params: type
    plotting: str
    draw: str
    made: Callable[[Dataset], list[Any]]
    what: Callable[[Any], str]
    inputs: Callable[[Any, Dataset], Inputs]
    check: Callable[[Any, Dataset], list[str]] = field(default=_no_check)
    unless: tuple[tuple[str, str], ...] = ()
    parsed: Callable[[str, Dataset], Any] | None = None


def slug(name: str) -> str:
    """A name as an output-name token: lower case, each run of other characters one ``_``.

    M5's lessons (``ctx.de.run._lessons``) and M6's gene sets (``ctx.function.run._slug``)
    carry their names in output names by this rule.
    """
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def unknown(name: str, value: object, known: Collection[Any], what: str, dataset: str) -> list[str]:
    """One problem when ``value`` is not among ``known``, else none."""
    return [] if value in known else [f"{name}: {value!r} is not {what} of {dataset}"]


def field_names(params: type) -> tuple[str, ...]:
    """A parameter dataclass's fields, in the order the params file writes them."""
    return tuple(f.name for f in fields(params))


@dataclass(frozen=True)
class Plain:
    """A figure the dataset file makes once: its dataset alone."""

    dataset: str


@dataclass(frozen=True)
class OneContrast:
    """A figure of one contrast."""

    dataset: str
    contrast: str


def featured(spec: Dataset) -> list[OneContrast]:
    """The featured contrast alone: M5 and M6 draw many of their figures for it only."""
    return [OneContrast(spec.slug, spec.de.featured)] if spec.de else []


def once(spec: Dataset) -> list[Plain]:
    """The one figure of a kind every dataset file makes."""
    return [Plain(spec.slug)]
