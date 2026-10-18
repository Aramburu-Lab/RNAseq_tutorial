"""Which figures have a spec: each kind's parameter dataclass, its inputs and its drawer.

A kind names every figure of its sort the dataset file can make (``made``), so a module's
figure stem is matched against what the file declares, never parsed by guesswork: contrast
names and dataset slugs both hold underscores. Drawers are named by path and imported when
used, since the modules that draw import ``ctx.viz.save``, which imports this file.
"""

import re
from collections.abc import Callable
from dataclasses import dataclass, fields
from typing import Any

from ctx.dataset import Dataset
from ctx.de.stats import ALPHA, LFC_MIN

# An input is a module's table, by (module, what after the dataset slug); or one of these.
DATASET_FILE = ("dataset_file", "")
QC_REGISTRY = ("qc_registry", "")


class NoFigureError(ValueError):
    """Parameters under which the module draws no figure, saying why."""


@dataclass(frozen=True)
class Sheet:
    """M1's figures: drawn from the curated sample sheet."""

    dataset: str


@dataclass(frozen=True)
class OneContrast:
    """A figure of one contrast."""

    dataset: str
    contrast: str


@dataclass(frozen=True)
class TwoContrasts:
    """One contrast's fold changes against another's."""

    dataset: str
    x: str
    y: str


@dataclass(frozen=True)
class Lesson:
    """One lesson's call counts at one setting of the call rule (``ctx.de.figures.Setting``)."""

    dataset: str
    lesson: str
    basis: str = "padj"
    alpha: float = ALPHA
    lfc_min: float = LFC_MIN


@dataclass(frozen=True)
class OneSet:
    """One gene set in every contrast."""

    dataset: str
    gene_set: str


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
    """

    module: str
    figure: str
    params: type
    plotting: str
    draw: str
    made: Callable[[Dataset], list[Any]]
    what: Callable[[Any], str]
    inputs: Callable[[Any, Dataset], dict[str, tuple[str, str]]]


def lesson_slug(title: str) -> str:
    """A lesson's title as M5's output names carry it (``ctx.de.run._lessons``)."""
    return re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_")


def set_slug(name: str) -> str:
    """A gene set's name as M6's output names carry it (``ctx.function.run._slug``)."""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def _contrasts(spec: Dataset) -> list[OneContrast]:
    """Every declared contrast as a figure's parameters."""
    return [OneContrast(spec.slug, c.name) for c in spec.contrasts]


def _featured(spec: Dataset) -> list[OneContrast]:
    """The featured contrast alone: M5 and M6 draw these figures for it only."""
    return [OneContrast(spec.slug, spec.de.featured)] if spec.de else []


def _sheet(_: Any, __: Dataset) -> dict[str, tuple[str, str]]:
    """M1 reads its own sheet and the dataset file."""
    return {"samplesheet": ("gather", "samplesheet"), "dataset_file": DATASET_FILE}


def _one(p: OneContrast, _: Dataset) -> dict[str, tuple[str, str]]:
    """One contrast's results, the size factors that name its samples, the dataset file."""
    return {"results": ("de", f"{p.contrast}_results"), "size_factors": ("de", "size_factors"),
            "dataset_file": DATASET_FILE}  # fmt: skip


def _lesson_names(p: Lesson, spec: Dataset) -> tuple[str, ...]:
    """The contrasts of one lesson, as the dataset file lists them; none for an unknown one."""
    lessons = spec.de.lessons if spec.de else {}
    return next((names for title, names in lessons.items() if title == p.lesson), ())


def _lesson_inputs(p: Lesson, spec: Dataset) -> dict[str, tuple[str, str]]:
    """Every contrast of the lesson's results, the size factors and the dataset file."""
    found = {f"results_{n}": ("de", f"{n}_results") for n in _lesson_names(p, spec)}
    return {**found, "size_factors": ("de", "size_factors"), "dataset_file": DATASET_FILE}


def _pair_inputs(p: TwoContrasts, _: Dataset) -> dict[str, tuple[str, str]]:
    """Both contrasts' results, the size factors and the dataset file."""
    return {"results_x": ("de", f"{p.x}_results"), "results_y": ("de", f"{p.y}_results"),
            "size_factors": ("de", "size_factors"), "dataset_file": DATASET_FILE}  # fmt: skip


def _gsea_inputs(p: OneContrast, _: Dataset) -> dict[str, tuple[str, str]]:
    """M6's GSEA table of the contrast and the dataset file."""
    return {"gsea": ("function", f"{p.contrast}_gsea"), "dataset_file": DATASET_FILE}


def _ora_inputs(p: OneContrast, _: Dataset) -> dict[str, tuple[str, str]]:
    """M6's ORA table, M5's results that give its background, and the dataset file."""
    return {"ora": ("function", f"{p.contrast}_ora"), "results": ("de", f"{p.contrast}_results"),
            "dataset_file": DATASET_FILE}  # fmt: skip


def _set_inputs(_: OneSet, spec: Dataset) -> dict[str, tuple[str, str]]:
    """Every contrast's GSEA table, in the dataset file's order, and the dataset file."""
    found = {f"gsea_{c.name}": ("function", f"{c.name}_gsea") for c in spec.contrasts}
    return {**found, "dataset_file": DATASET_FILE}


def _lessons(spec: Dataset) -> list[Lesson]:
    """Every lesson of ``[de.lessons]`` at the build's setting."""
    return [Lesson(spec.slug, title) for title in (spec.de.lessons if spec.de else {})]


def _pairs(spec: Dataset) -> list[TwoContrasts]:
    """Every pair of ``[de] scatters``."""
    return [TwoContrasts(spec.slug, x, y) for x, y in (spec.de.scatters if spec.de else ())]


def _sets(spec: Dataset) -> list[OneSet]:
    """The featured gene set of ``[function]``."""
    return [OneSet(spec.slug, spec.function.featured_set)] if spec.function else []


KINDS = (
    Kind("gather", "composition", Sheet, "ctx.gather.harmonise:_composition",
         "ctx.figures.draw_m1:composition", lambda s: [Sheet(s.slug)],
         lambda _: "composition", _sheet),
    Kind("gather", "library_sizes", Sheet, "ctx.gather.harmonise:_library_sizes",
         "ctx.figures.draw_m1:library_sizes", lambda s: [Sheet(s.slug)],
         lambda _: "library_sizes",
         lambda p, s: {**_sheet(p, s), "qc_registry": QC_REGISTRY}),
    Kind("de", "volcano", OneContrast, "ctx.de.run:volcano", "ctx.figures.draw_m5:volcano",
         _featured, lambda p: f"volcano_{p.contrast}", _one),
    Kind("de", "ma", OneContrast, "ctx.de.run:_ma", "ctx.figures.draw_m5:ma", _featured,
         lambda p: f"ma_{p.contrast}", _one),
    Kind("de", "phist", OneContrast, "ctx.de.run:p_histogram", "ctx.figures.draw_m5:phist",
         _featured, lambda p: f"phist_{p.contrast}", _one),
    Kind("de", "lfc", TwoContrasts, "ctx.de.figures:lfc_scatter", "ctx.figures.draw_m5:lfc",
         _pairs, lambda p: f"lfc_{p.x}_vs_{p.y}", _pair_inputs),
    Kind("de", "calls", Lesson, "ctx.de.figures:calls_bars", "ctx.figures.draw_m5:calls",
         _lessons, lambda p: f"calls_{lesson_slug(p.lesson)}", _lesson_inputs),
    Kind("function", "gsea", OneContrast, "ctx.function.run:_gsea_bars",
         "ctx.figures.draw_m6:gsea", _featured, lambda p: f"gsea_{p.contrast}", _gsea_inputs),
    Kind("function", "ora", OneContrast, "ctx.function.run:_ora_bars", "ctx.figures.draw_m6:ora",
         _featured, lambda p: f"ora_{p.contrast}", _ora_inputs),
    Kind("function", "nes", OneSet, "ctx.function.run:_featured_bars",
         "ctx.figures.draw_m6:nes", _sets, lambda p: f"nes_{set_slug(p.gene_set)}", _set_inputs),
)  # fmt: skip


def kind(module: str, figure: str) -> Kind | None:
    """The kind of one module's figure name, or None when it has no spec yet."""
    return next((k for k in KINDS if (k.module, k.figure) == (module, figure)), None)


def matched(module: str, what: str, spec: Dataset) -> tuple[Kind, Any] | None:
    """The kind and parameters of a figure a module writes as ``<dataset>_<what>``, or None."""
    for k in KINDS:
        if k.module == module:
            for params in k.made(spec):
                if k.what(params) == what:
                    return k, params
    return None


def field_names(params: type) -> tuple[str, ...]:
    """A parameter dataclass's fields, in the order the params file writes them."""
    return tuple(f.name for f in fields(params))
