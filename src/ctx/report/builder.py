"""Render self-contained HTML pages: plotly.js inlined once, no external request, light/dark.

Three page kinds share one base template, so there is one theme toggle and one look: a figure
page (one figure with its reading guide, caption and data), a module report with the eight
sections of the report skeleton (CONVENTIONS.md, HTML reports), and the course hub.
"""

import json
import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from jinja2 import Environment, PackageLoader, select_autoescape
from markupsafe import Markup, escape

from ctx.io.naming import figure_files
from ctx.report import text
from ctx.viz.style import (
    NARROW_PX,
    RAIL_FROM_PX,
    THEMES,
    css_layout,
    css_variables,
    plotly_templates,
)

MAX_BYTES = 15 * 10**6
# The report skeleton, in the order CONVENTIONS.md fixes; ids are the section anchors, and
# the headings are in content/ui.toml.
SECTIONS = tuple(
    (section, text.word("ui", f"sections.{section}"))
    for section in (
        "purpose", "input", "method", "parameters", "results", "interpretation", "questions",
        "pitfalls", "reproduce",
    )
)  # fmt: skip
# The ModulePage field behind each section. The frozen M5 and M6 pages still pass the old
# names; P2S.4 renames the fields once the Step-1 record is taken (D-045).
FIELDS = {"purpose": "question", "method": "concept", "interpretation": "practice",
          "questions": "think", "pitfalls": "mistakes"}  # fmt: skip
MAX_TABLE_COLUMNS = 8  # CONVENTIONS.md, HTML reports
MAX_KEY_NUMBERS = 4
# The packages a reader asks about first; the full list sits behind a disclosure.
KEY_PACKAGES = ("python", "ctx", "polars", "numpy", "scipy", "pydeseq2", "plotly", "matplotlib")
PLOTLY_CONFIG = {
    "displaylogo": False,
    "responsive": True,
    # No toolbar: over a phone-width plot it covers the legend, and drag-to-zoom and
    # double-click-to-reset work without it. A PNG export would carry no statistics anyway:
    # the caption sits beside the plot, not in it.
    "displayModeBar": False,
}
HEX = re.compile(r"#[0-9A-Fa-f]{6}")
TITLE = re.compile(r"(?:M(?P<n>\d+)[a-z]? )?(?P<headline>[^:]+?)(?:: (?P<dataset>.+))?")

# The top bar names the course on every page.
SUITE = text.word("ui", "suite")


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
class FigureText:
    """The words around a figure on a page, each filled from the run's own numbers.

    Attributes:
        title: What the figure shows, naming the contrast or subset.
        lead: One line: the finding, or what to look for.
        reading: "How to read this": what a mark is and the rule that colours it.
        tab: Tab label inside a tabbed group; the title when empty.
        back: File name of the module page, linked from the standalone figure page.
    """

    title: str = ""
    lead: str = ""
    reading: tuple[str, ...] = ()
    tab: str = ""
    back: str = ""


@dataclass(frozen=True)
class FigureBlock:
    """A plotly figure with its caption, its data rows and the name of its CSV.

    Attributes:
        figure: The plotly Figure.
        caption: The statistics caption, printed once, under the plot.
        rows: The rows the page draws, embedded for the reader and the page script.
        data_file: The CSV beside the page.
        text: Title, lead and reading guide.
        stem: The figure's file stem, which names its PDF, SVG and README; empty for none.
        controls: Settings for the page script, with a ``kind``: "calls" for live calling,
            "brush" for linked brushing, "filter" for bars filtered by column; or None.
    """

    figure: Any
    caption: str
    rows: Sequence[Mapping[str, Any]]
    data_file: str
    text: FigureText = FigureText()
    stem: str = ""
    controls: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class ResultTable:
    """A small table of results; at most ``MAX_TABLE_COLUMNS`` columns, cells as text."""

    title: str
    columns: tuple[str, ...]
    rows: tuple[tuple[str, ...], ...]
    numeric: tuple[int, ...] = ()
    swatch: int | None = None

    def __post_init__(self) -> None:
        """Refuse a table the page could only print wrongly."""
        if len(self.columns) > MAX_TABLE_COLUMNS:
            raise ValueError(f"{self.title}: {len(self.columns)} columns; at most 8")
        if ragged := [r for r in self.rows if len(r) != len(self.columns)]:
            raise ValueError(f"{self.title}: rows of the wrong length: {ragged[:2]}")
        if self.swatch is not None and not all(HEX.fullmatch(r[self.swatch]) for r in self.rows):
            raise ValueError(f"{self.title}: a swatch cell is not #RRGGBB")


@dataclass(frozen=True)
class Group:
    """Results that answer one question, under one heading; tabs for alternative views."""

    title: str
    items: tuple["FigureBlock | ResultTable | str", ...]
    lead: str = ""
    tabs: bool = False

    def __post_init__(self) -> None:
        """Refuse an empty group, and tabs over anything but two or more figures."""
        if not self.items:
            raise ValueError(f"group {self.title!r} is empty")
        figures = all(isinstance(i, FigureBlock) for i in self.items)
        if self.tabs and not (figures and len(self.items) > 1):
            raise ValueError(f"group {self.title!r}: tabs need two or more figures only")


@dataclass(frozen=True)
class KeyNumber:
    """One number that summarises the run, computed in it."""

    label: str
    value: str
    detail: str = ""


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
    inputs: tuple[str, ...] = ()


@dataclass(frozen=True)
class FigurePage:
    """A page holding one figure."""

    title: str
    figure: FigureBlock


@dataclass(frozen=True)
class HubLink:
    """One module's newest page, linked relative to the hub."""

    module: str
    title: str
    href: str
    written: str


@dataclass(frozen=True)
class HubRule:
    """One sanity rule: its verdict, and each condition with the value it measured."""

    rule: str
    text: str
    verdict: str
    conditions: tuple[tuple[str, str, str], ...]


@dataclass(frozen=True)
class HubPage:
    """The course hub: one link per module page of one dataset."""

    title: str
    intro: str
    parameters: tuple[ParameterCard, ...]
    links: tuple[HubLink, ...]
    missing: tuple[str, ...]
    rules: tuple[HubRule, ...] = ()
    key_numbers: tuple[KeyNumber, ...] = ()


@dataclass(frozen=True)
class CourseLink:
    """One page on the course home page, and whether its data are real or made up."""

    label: str
    href: str
    real: bool


@dataclass(frozen=True)
class CourseRow:
    """One module of the course, with a link to each dataset's page of it."""

    number: str
    name: str
    links: tuple[CourseLink, ...]


@dataclass(frozen=True)
class CoursePage:
    """The course home page of a kit: every module in course order, and each dataset's hub."""

    title: str
    rows: tuple[CourseRow, ...]
    hubs: tuple[CourseLink, ...]
    key_numbers: tuple[KeyNumber, ...] = ()


Result = Group | FigureBlock | ResultTable | str


@dataclass(frozen=True)
class ModulePage:
    """A module report: every section of the skeleton is required and non-empty."""

    title: str
    question: str
    concept: str
    parameters: tuple[ParameterCard, ...]
    results: tuple[Result, ...]
    practice: str
    think: tuple[SelfCheck, ...]
    mistakes: tuple[str, ...]
    reproduce: Reproduction
    key_numbers: tuple[KeyNumber, ...] = ()
    terms: tuple[tuple[str, str], ...] = ()
    input: str = ""

    def __post_init__(self) -> None:
        """Refuse a page with an empty section; the skeleton is the teaching contract.

        Key numbers and terms are optional: a page may have nothing to summarise. Input is
        optional until the frozen M5 and M6 pages gain one after the Step-1 record (D-045).
        """
        sections = (FIELDS.get(s, s) for s, _ in SECTIONS if s != "input")
        required = ("title", *sections)
        if empty := [name for name in required if not getattr(self, name)]:
            raise ValueError(f"empty report sections: {empty}")
        blocks = [i for r in self.results for i in (r.items if isinstance(r, Group) else (r,))]
        if unread := [b.text.title or b.stem for b in blocks
                      if isinstance(b, FigureBlock) and not b.text.reading]:  # fmt: skip
            raise ValueError(f"figures without a reading guide: {unread}")
        if len(self.key_numbers) > MAX_KEY_NUMBERS:
            raise ValueError(f"{len(self.key_numbers)} key numbers; at most {MAX_KEY_NUMBERS}")


Page = FigurePage | ModulePage | HubPage | CoursePage

_ENV = Environment(
    loader=PackageLoader("ctx.report", "templates"),
    autoescape=select_autoescape(default=True, default_for_string=True),
    trim_blocks=True,
    lstrip_blocks=True,
)
_CODE = re.compile(r"`([^`]+)`")


def code_spans(value: str) -> Markup:
    """Page text with each backticked span as code, everything else escaped (D-045)."""
    return Markup(_CODE.sub(lambda m: f"<code>{m.group(1)}</code>", str(escape(value))))


_ENV.filters["code"] = code_spans


def plotly_js() -> str:
    r"""plotly.js as a page embeds it, its one raw U+FFFD written as the escape ``\\uFFFD``.

    The character sits in a regex literal, where JavaScript reads the escape as the same
    character; a host that takes a raw replacement character for damaged text refuses the page.
    """
    import plotly.offline  # noqa: PLC0415 - heavy; imported after the runtime profile

    library: str = plotly.offline.get_plotlyjs()
    return library.replace("\ufffd", "\\uFFFD")


def _json_for_script(value: Any) -> Markup:
    """JSON safe inside a ``<script>`` element: every ``<`` is escaped as a unicode escape.

    ``allow_nan=False``: a browser's JSON.parse rejects NaN and Infinity, which Python writes
    by default; rows go through ``_finite`` first.
    """
    return Markup(json.dumps(value, default=str, allow_nan=False).replace("<", "\\u003c"))


def _finite(row: Mapping[str, Any]) -> dict[str, Any]:
    """A row with every NaN or infinite number as null, which is what they mean on a page."""
    return {k: None if isinstance(v, float) and not math.isfinite(v) else v for k, v in row.items()}


def _figure_html(block: FigureBlock, index: int) -> Markup:
    """The plotly div and its drawing script, without plotly.js itself."""
    html = block.figure.to_html(
        include_plotlyjs=False, full_html=False, div_id=f"ctx-figure-{index}", config=PLOTLY_CONFIG
    )
    return Markup(html)


def _files(block: FigureBlock) -> list[tuple[str, str, str]]:
    """(label, theme or "", file name) of every file a reader may open beside the figure."""
    data = text.word("ui", "figures.data")
    if not block.stem:
        return [(data, "", block.data_file)]
    names = figure_files(block.stem)
    rendered = [
        (ext.upper(), theme, name)
        for name in names
        for theme in THEMES
        for ext in ("pdf", "svg")
        if name == f"{block.stem}_{theme}.{ext}"
    ]
    return [
        (data, "", block.data_file),
        (text.word("ui", "figures.columns"), "", f"{block.stem}_README.md"),
        *rendered,
    ]


def _prepare(block: FigureBlock, index: int, number: int, level: int) -> dict[str, Any]:
    """Everything the figure macro prints for one figure."""
    return {
        "html": _figure_html(block, index),
        "caption": block.caption,
        "rows": _json_for_script([_finite(r) for r in block.rows]),
        "data_file": block.data_file,
        "index": index,
        "number": number,
        "level": level,
        "text": block.text,
        "title": block.text.title,
        "tab": block.text.tab or block.text.title or text.fill("ui", "figure_n", n=number),
        "page": f"{block.stem}.html" if block.stem else "",
        "files": _files(block),
        "controls": _json_for_script(block.controls) if block.controls else None,
        "kind": block.controls["kind"] if block.controls else None,
        "defaults": _defaults(block.controls) if block.controls else None,
    }


def _defaults(controls: Mapping[str, Any]) -> dict[str, Any]:
    """What the controls show before the reader moves them: the default rule, no paint, or all.

    A filter's live line is ``ui.toml`` [filter] summary, filled here with the numbers of every
    bar and passed unfilled for the page script to fill again (``ctxFilter.summary``).
    """
    if controls["kind"] == "brush":
        return {"colours": list(zip(controls["colours"], controls["names"], strict=True)),
                "empty": controls["empty"]}  # fmt: skip
    if controls["kind"] == "filter":
        return {"filters": list(controls["filters"].items()),
                "summary": text.fill("ui", "filter.summary", **controls["start"]),
                "template": text.raw("ui", "filter.summary"),
                "none": text.word("ui", "filter.none")}  # fmt: skip
    return {
        "fc": f"{controls['lfc_min']:g}",
        "bases": list(controls["bases"].items()),
        "alpha": f"{controls['alpha']:g}",
        "presets": sorted({*controls["presets"], f"{controls['min_cutoff']:g}"}, key=float),
        "summary": controls["default_summary"],
    }


def _walk(items: Iterable[Result]) -> tuple[list[tuple[str, Any]], list[tuple[str, str]], bool]:
    """Results ready for the template, the contents entries of the groups, and any figure.

    Figures are numbered in page order across groups; each group gets an anchor.
    """
    out: list[tuple[str, Any]] = []
    outline: list[tuple[str, str]] = []
    counter = 0

    def one(item: FigureBlock | ResultTable | str, level: int) -> tuple[str, Any]:
        """One item, prepared."""
        nonlocal counter
        if isinstance(item, FigureBlock):
            counter += 1
            return "figure", _prepare(item, counter - 1, counter, level)
        return ("table", item) if isinstance(item, ResultTable) else ("text", item)

    for item in items:
        if isinstance(item, Group):
            anchor = f"group-{len(outline) + 1}"
            outline.append((anchor, item.title))
            prepared = [one(i, 4) for i in item.items]
            out.append(("group", {"anchor": anchor, "group": item, "items": prepared}))
        else:
            out.append(one(item, 3))
    return out, outline, counter > 0


def _heading(page: Page) -> dict[str, str]:
    """Kicker, headline and dataset, parsed from the page title."""
    if isinstance(page, FigurePage):
        figure = text.word("ui", "figure")
        return {"kicker": figure, "headline": page.title, "dataset": "", "suite": SUITE}
    match = TITLE.fullmatch(page.title)
    assert match is not None  # noqa: S101 - the pattern accepts any title
    if isinstance(page, HubPage | CoursePage):
        kicker = SUITE
    else:
        numbered = match["n"] is not None
        kicker = (
            text.fill("ui", "module_n", n=match["n"]) if numbered else text.word("ui", "module")
        )
    dataset = match["dataset"] or ""
    return {"kicker": kicker, "headline": match["headline"], "dataset": dataset, "suite": SUITE}


def render(page: Page) -> str:
    """Render a page to one self-contained HTML string.

    Args:
        page: A figure page, a module report, the hub or the course home page.

    Returns:
        The HTML.
    """
    context: dict[str, Any] = {
        "page": page,
        "sections": SECTIONS,
        "ui": text.load("ui"),
        "heading": _heading(page),
        # Our own constants, and CSS inside <style> must not be HTML-escaped.
        "css": {
            **{name: Markup(css_variables(theme)) for name, theme in THEMES.items()},
            "layout": Markup(css_layout()),
        },
        "layout": {"rail_from": RAIL_FROM_PX, "narrow": NARROW_PX},
        "templates": _json_for_script(plotly_templates()),
        "plotly_js": Markup(plotly_js()),
    }
    if isinstance(page, HubPage):
        context["plotly_js"] = ""  # 4.8 MB of library for a page that draws nothing
        return _ENV.get_template("hub.html.j2").render(**context)
    if isinstance(page, CoursePage):
        context |= {"plotly_js": "", "course": text.load("course")}
        return _ENV.get_template("course.html.j2").render(**context)
    if isinstance(page, FigurePage):
        context["figure"] = _prepare(page.figure, 0, 1, 3)
        return _ENV.get_template("figure.html.j2").render(**context)
    results, outline, draws = _walk(page.results)
    if not draws:
        context["plotly_js"] = ""
    versions = page.reproduce.versions
    context |= {
        "results": results,
        "outline": outline,
        "key_versions": [(n, versions[n]) for n in KEY_PACKAGES if n in versions],
    }
    return _ENV.get_template("module.html.j2").render(**context)


def write(page: Page, path: Path) -> Path:
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
