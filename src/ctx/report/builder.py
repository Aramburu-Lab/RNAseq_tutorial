"""Render self-contained HTML pages: plotly.js inlined once, no external request, light/dark.

Two page kinds share one base template, so there is one theme toggle: a figure page (the
figure, its caption and its data) and a module report with the eight sections of the report
skeleton (CONVENTIONS.md, HTML reports).
"""

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any

from jinja2 import Environment, PackageLoader, select_autoescape
from markupsafe import Markup

from ctx.viz.style import THEMES, css_variables, plotly_templates

MAX_BYTES = 15 * 10**6
# The report skeleton, in the order CONVENTIONS.md fixes; ids are the section anchors.
SECTIONS = (
    ("question", "The big question"),
    ("concept", "The concept in plain English"),
    ("parameters", "Parameters"),
    ("results", "Results"),
    ("practice", "What it means in practice"),
    ("think", "Think about it"),
    ("mistakes", "Common mistakes"),
    ("reproduce", "Reproduce this page"),
)


@dataclass(frozen=True)
class ParameterCard:
    """One CLI parameter: what it means, its default, and what moving it changes."""

    name: str
    meaning: str
    default: str
    effect: str


@dataclass(frozen=True)
class SelfCheck:
    """A think-about-it prompt with its answer hidden in ``<details>``."""

    prompt: str
    answer: str


@dataclass(frozen=True)
class FigureBlock:
    """A plotly figure with its caption, its data rows and the name of its CSV."""

    figure: Any
    caption: str
    rows: Sequence[Mapping[str, Any]]
    data_file: str


@dataclass(frozen=True)
class Reproduction:
    """How to get this page again, and what it cost."""

    command: str
    wall_seconds: float
    peak_rss_self: int
    peak_rss_children: int
    profile: str
    budget: str
    versions: Mapping[str, str]


@dataclass(frozen=True)
class FigurePage:
    """A page holding one figure."""

    title: str
    figure: FigureBlock


@dataclass(frozen=True)
class ModulePage:
    """A module report: every section of the skeleton is required and non-empty."""

    title: str
    question: str
    concept: str
    parameters: tuple[ParameterCard, ...]
    results: tuple[FigureBlock | str, ...]
    practice: str
    think: tuple[SelfCheck, ...]
    mistakes: tuple[str, ...]
    reproduce: Reproduction

    def __post_init__(self) -> None:
        """Refuse a page with an empty section; the skeleton is the teaching contract."""
        if empty := [f.name for f in fields(self) if not getattr(self, f.name)]:
            raise ValueError(f"empty report sections: {empty}")


_ENV = Environment(
    loader=PackageLoader("ctx.report", "templates"),
    autoescape=select_autoescape(default=True, default_for_string=True),
    trim_blocks=True,
    lstrip_blocks=True,
)


def _json_for_script(value: Any) -> Markup:
    """JSON safe inside a ``<script>`` element: every ``<`` is escaped as a unicode escape."""
    return Markup(json.dumps(value, default=str).replace("<", "\\u003c"))


def _figure_html(block: FigureBlock, index: int) -> Markup:
    """The plotly div and its drawing script, without plotly.js itself."""
    html = block.figure.to_html(
        include_plotlyjs=False, full_html=False, div_id=f"ctx-figure-{index}",
        config={"displaylogo": False, "responsive": True},
    )  # fmt: skip
    return Markup(html)


def _prepare(block: FigureBlock, index: int) -> dict[str, Any]:
    """Everything the figure macro prints for one figure."""
    return {
        "html": _figure_html(block, index),
        "caption": block.caption,
        "rows": _json_for_script(list(block.rows)),
        "data_file": block.data_file,
        "index": index,
    }


def render(page: FigurePage | ModulePage) -> str:
    """Render a page to one self-contained HTML string.

    Args:
        page: A figure page or a module report.

    Returns:
        The HTML.
    """
    import plotly.offline  # noqa: PLC0415 - heavy; imported after the runtime profile

    context: dict[str, Any] = {
        "page": page,
        "sections": SECTIONS,
        # Our own constants, and CSS inside <style> must not be HTML-escaped.
        "css": {name: Markup(css_variables(theme)) for name, theme in THEMES.items()},
        "templates": _json_for_script(plotly_templates()),
        "plotly_js": Markup(plotly.offline.get_plotlyjs()),
    }
    if isinstance(page, FigurePage):
        context["figure"] = _prepare(page.figure, 0)
        return _ENV.get_template("figure.html.j2").render(**context)
    context["results"] = [
        ("figure", _prepare(item, i)) if isinstance(item, FigureBlock) else ("text", item)
        for i, item in enumerate(page.results)
    ]
    return _ENV.get_template("module.html.j2").render(**context)


def write(page: FigurePage | ModulePage, path: Path) -> Path:
    """Render and write a page, refusing one over 15 MB before anything is written.

    Args:
        page: The page.
        path: Destination ``.html`` file.

    Returns:
        The written path.

    Raises:
        ValueError: If the page exceeds ``MAX_BYTES``.
    """
    html = render(page).encode()
    if len(html) > MAX_BYTES:
        raise ValueError(f"{path.name}: {len(html)} bytes, over the {MAX_BYTES}-byte cap")
    path.write_bytes(html)
    return path
