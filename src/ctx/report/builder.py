"""Render self-contained HTML pages: plotly.js inlined once, no external request, light/dark.

Every page kind shares one base template, so there is one theme toggle and one look: a figure
page (one figure with its reading guide, caption and data), a module report with the nine
sections of the report skeleton (CONVENTIONS.md, HTML reports), the course hub, the course home
page and the glossary page. A page that marks a glossary term links the glossary page and
embeds the tooltips of the terms it marks, and writing it writes the glossary too (D-048).
"""

import base64
import gzip
import hashlib
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from jinja2 import Environment, PackageLoader, pass_context, select_autoescape
from markupsafe import Markup

from ctx.report import export, glossary, marks, text
from ctx.report.figures import json_for_script, prepare
from ctx.viz.style import (
    ALONE_SCALE,
    NARROW_PX,
    RAIL_FROM_PX,
    THEMES,
    css_layout,
    css_variables,
    plotly_templates,
    zoom,
)

if TYPE_CHECKING:
    from ctx.viz.aggregates import Lead

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
MAX_TABLE_COLUMNS = 8  # CONVENTIONS.md, HTML reports
MAX_KEY_NUMBERS = 4
# The packages a reader asks about first; the full list sits behind a disclosure.
KEY_PACKAGES = ("python", "ctx", "polars", "numpy", "scipy", "pydeseq2", "plotly", "matplotlib")
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
        follows: A lead whose numbers are named aggregates over the rows the figure draws, in
            place of ``lead``: ``ctx.viz.save.figure_block`` fills it, and a figure with
            filters refills them over the rows the reader keeps.
    """

    title: str = ""
    lead: str = ""
    reading: tuple[str, ...] = ()
    tab: str = ""
    back: str = ""
    follows: "Lead | None" = None

    def __post_init__(self) -> None:
        """Refuse a lead given twice, as text and as numbers to fill."""
        if self.lead and self.follows is not None:
            raise ValueError("a lead is text or follows its rows, not both")


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
            "brush" for linked brushing, "filter" for bars filtered by column; or None. Its
            ``export`` entry draws the PDF and SVG of a changed view (``ctx.report.export``).
        below: A second plotly figure, drawn under the plot and its reading guide across the
            card: a linked view's parallel sets; or None.
    """

    figure: Any
    caption: str
    rows: Sequence[Mapping[str, Any]]
    data_file: str
    text: FigureText = FigureText()
    stem: str = ""
    controls: Mapping[str, Any] | None = None
    below: Any = None


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
class Code:
    """Text to run, shown exactly as given, with its SHA-256 and a Copy button.

    Attributes:
        title: What the block does, or the file it is.
        text: The lines, as a file holds them or a reader types them; never reflowed.
        source: The file the text is, named beside its hash; empty for typed commands.
    """

    title: str
    text: str
    source: str = ""

    @property
    def sha256(self) -> str:
        """The SHA-256 of the text's UTF-8 bytes, as ``sha256sum`` prints it for the file."""
        return hashlib.sha256(self.text.encode()).hexdigest()


@dataclass(frozen=True)
class Group:
    """Results that answer one question, under one heading; tabs for alternative views.

    A packed group's items travel gzipped in the page and are drawn when the reader reaches
    them (``templates/_page.js.j2``, packed): M9's 50 TCGA cards drawn whole made a
    29,472,416-byte page, where their figures' JSON and rows compress about eight times (C3.1).
    """

    title: str
    items: tuple["FigureBlock | ResultTable | Code | str", ...]
    lead: str = ""
    tabs: bool = False
    packed: bool = False

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
    glossary: str = ""
    guide: str = ""


@dataclass(frozen=True)
class GlossaryPage:
    """The glossary: every entry of ``glossary.yaml``, by category, each anchored by its id."""

    title: str
    groups: tuple[tuple[str, tuple[glossary.Entry, ...]], ...]


Result = Group | FigureBlock | ResultTable | Code | str


@dataclass(frozen=True)
class ModulePage:
    """A module report: every section of the skeleton is required and non-empty."""

    title: str
    purpose: str
    input: str
    method: str
    parameters: tuple[ParameterCard, ...]
    results: tuple[Result, ...]
    interpretation: str
    questions: tuple[SelfCheck, ...]
    pitfalls: tuple[str, ...]
    reproduce: Reproduction
    key_numbers: tuple[KeyNumber, ...] = ()
    terms: tuple[tuple[str, str], ...] = ()
    kicker: str = ""

    def __post_init__(self) -> None:
        """Refuse a page with an empty section; the skeleton is the teaching contract.

        Each section is the field of its key (P2P.3), every one required, Input too. Key
        numbers and terms are optional: a page may have nothing to summarise.
        """
        required = ("title", *(s for s, _ in SECTIONS))
        if empty := [name for name in required if not getattr(self, name)]:
            raise ValueError(f"empty report sections: {empty}")
        blocks = [i for r in self.results for i in (r.items if isinstance(r, Group) else (r,))]
        if unread := [b.text.title or b.stem for b in blocks
                      if isinstance(b, FigureBlock) and not b.text.reading]:  # fmt: skip
            raise ValueError(f"figures without a reading guide: {unread}")
        if len(self.key_numbers) > MAX_KEY_NUMBERS:
            raise ValueError(f"{len(self.key_numbers)} key numbers; at most {MAX_KEY_NUMBERS}")


@dataclass(frozen=True)
class EntryPage:
    """A module folder's entry page, ``index.html``: it forwards a link's gene list (GS.2).

    Pages are dated and a module is rebuilt after the pages that link it, so a link names the
    folder's entry page, which each build of the module rewrites to its newest pages.

    Attributes:
        title: The page's title.
        links: (dataset, href) of each dataset's newest overview in the folder.
    """

    title: str
    links: tuple[tuple[str, str], ...]


Page = FigurePage | ModulePage | HubPage | CoursePage | GlossaryPage | EntryPage

_ENV = Environment(
    loader=PackageLoader("ctx.report", "templates"),
    autoescape=select_autoescape(default=True, default_for_string=True),
    trim_blocks=True,
    lstrip_blocks=True,
)


def code_spans(value: str) -> Markup:
    """Page text with each backticked span as code, everything else escaped (D-045).

    A glossary mark shows its display text: only ``explain_terms`` draws the icon.
    """
    return marks.html(value)


@pass_context
def explain_terms(context: Any, value: str, section: str) -> Markup:
    """Page text with each backticked span as code and each glossary mark as a term (D-048).

    The first mention of an entry in a section carries the info icon; a later one shows its
    display text alone, unless its mark ends in ``!``.
    """
    terms: glossary.Terms = context["terms"]
    return terms.html(value, section)


def pack(html: str) -> str:
    """A packed group's HTML as the page carries it: gzip with no timestamp, then base64.

    ``mtime=0``: a timestamp in the header would make two builds of one page differ.
    """
    return base64.b64encode(gzip.compress(html.encode(), mtime=0)).decode("ascii")


PACKED = re.compile(
    r'<script type="application/octet-stream" id="ctx-packed-([0-9]+)">([^<]*)</script>'
)


def packed_bytes(html: str) -> list[int]:
    """The bytes of each packed group in a rendered page, in page order: what each one costs."""
    return [len(m[0].encode()) for m in PACKED.finditer(html)]


def unpacked(html: str) -> str:
    """A page with each packed group's HTML in place of its script, as the reader sees it drawn.

    For whatever reads a page's words or links without a browser: a test, the kit's link walk.
    """
    return PACKED.sub(lambda m: gzip.decompress(base64.b64decode(m[2])).decode(), html)


_ENV.filters["code"] = code_spans
_ENV.filters["explain_terms"] = explain_terms
_ENV.filters["pack"] = pack


def _marked(page: Page) -> list[str]:
    """Every text of a page that ``explain_terms`` renders: sections, captions and cards."""
    if isinstance(page, FigurePage):
        return [page.figure.caption]
    cards = getattr(page, "parameters", ())
    texts = [t for card in cards for t in (card.meaning, card.effect)]
    if not isinstance(page, ModulePage):
        return texts if isinstance(page, HubPage) else []
    blocks = [i for r in page.results for i in (r.items if isinstance(r, Group) else (r,))]
    return [
        *texts, page.purpose, page.input, page.method, page.interpretation, *page.pitfalls,
        *(t for check in page.questions for t in (check.prompt, check.answer)),
        *(b.caption for b in blocks if isinstance(b, FigureBlock)),
    ]  # fmt: skip


def plotly_js() -> str:
    r"""plotly.js as a page embeds it, its one raw U+FFFD written as the escape ``\\uFFFD``.

    The character sits in a regex literal, where JavaScript reads the escape as the same
    character; a host that takes a raw replacement character for damaged text refuses the page.
    """
    import plotly.offline  # noqa: PLC0415 - heavy; imported after the runtime profile

    library: str = plotly.offline.get_plotlyjs()
    return library.replace("\ufffd", "\\uFFFD")


def _walk(items: Iterable[Result]) -> tuple[list[tuple[str, Any]], list[tuple[str, str]], bool]:
    """Results ready for the template, the contents entries of the groups, and any figure.

    Figures are numbered in page order across groups; each group gets an anchor.
    """
    out: list[tuple[str, Any]] = []
    outline: list[tuple[str, str]] = []
    counter = 0
    codes = 0

    def one(item: FigureBlock | ResultTable | Code | str, level: int) -> tuple[str, Any]:
        """One item, prepared; each code block gets its own id, which its Copy button names."""
        nonlocal counter, codes
        if isinstance(item, FigureBlock):
            counter += 1
            return "figure", prepare(item, counter - 1, counter, level)
        if isinstance(item, Code):
            codes += 1
            return "code", {"id": f"ctx-code-{codes}", "code": item}
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


def _live(prepared: Iterable[Mapping[str, Any]]) -> bool:
    """Whether a page scores gene lists live.

    Such a page carries the selection and survival scripts (``_select``, ``_surv`` and
    ``_live``); other pages leave them out.
    """
    return any(f.get("kind") in ("survival", "follow") for f in prepared)


def _heading(page: Page) -> dict[str, str]:
    """Kicker, headline and dataset, parsed from the page title."""
    if isinstance(page, FigurePage):
        figure = text.word("ui", "figure")
        return {"kicker": figure, "headline": page.title, "dataset": "", "suite": SUITE}
    match = TITLE.fullmatch(page.title)
    assert match is not None  # noqa: S101 - the pattern accepts any title
    if isinstance(page, HubPage | CoursePage | GlossaryPage | EntryPage):
        kicker = SUITE
    elif page.kicker:
        kicker = page.kicker
    else:
        numbered = match["n"] is not None
        kicker = (
            text.fill("ui", "module_n", n=match["n"]) if numbered else text.word("ui", "module")
        )
    dataset = match["dataset"] or ""
    return {"kicker": kicker, "headline": match["headline"], "dataset": dataset, "suite": SUITE}


def render(page: Page, glossary_href: str = "") -> str:
    """Render a page to one self-contained HTML string.

    Args:
        page: A figure page, a module report, the hub, the course home page or the glossary.
        glossary_href: The glossary page as this page links it, which ``write`` supplies; the
            course home page carries its own. Empty, the page and its tooltips link none.

    Returns:
        The HTML.
    """
    used = glossary.used(_marked(page))
    glossary_href = page.glossary if isinstance(page, CoursePage) else glossary_href
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
        "templates": json_for_script(plotly_templates()),
        "plotly_js": Markup(plotly_js()),
        "terms": glossary.Terms(),
        "glossary_href": glossary_href,
        # Only the entries this page marks: a page carries its own tooltips, and no more.
        "tooltips": json_for_script(glossary.tooltips(used, glossary_href)) if used else "",
    }
    if isinstance(page, GlossaryPage):
        names = {e.id: e.term for _, entries in page.groups for e in entries}
        context |= {"plotly_js": "", "names": names}
        return _ENV.get_template("glossary.html.j2").render(**context)
    if isinstance(page, EntryPage):
        context["plotly_js"] = ""
        return _ENV.get_template("entry.html.j2").render(**context)
    if isinstance(page, HubPage):
        context["plotly_js"] = ""  # 4.8 MB of library for a page that draws nothing
        return _ENV.get_template("hub.html.j2").render(**context)
    if isinstance(page, CoursePage):
        context |= {"plotly_js": "", "course": text.load("course")}
        return _ENV.get_template("course.html.j2").render(**context)
    if isinstance(page, FigurePage):
        # Opened alone, at ALONE_SCALE, the theme toggle swaps templates at that scale too, and
        # at the page's own where the figure is drawn at the page's size (figures.page_layout).
        alone = {name: zoom(t, ALONE_SCALE) for name, t in plotly_templates().items()}
        figure = prepare(page.figure, 0, 1, 3, ALONE_SCALE)
        context |= {"figure": figure, "templates": json_for_script(alone),
                    "page_templates": context["templates"],
                    "export": export.assets([figure]), "live": _live([figure])}  # fmt: skip
        return _ENV.get_template("figure.html.j2").render(**context)
    results, outline, draws = _walk(page.results)
    if not draws:
        context["plotly_js"] = ""
    versions = page.reproduce.versions
    items = [i for k, it in results for i in (it["items"] if k == "group" else [(k, it)])]
    prepared = [it for k, it in items if k == "figure"]
    context |= {
        "results": results,
        "outline": outline,
        "key_versions": [(n, versions[n]) for n in KEY_PACKAGES if n in versions],
        "export": export.assets(prepared),
        "live": _live(prepared),
    }
    return _ENV.get_template("module.html.j2").render(**context)


def rendered(page: Page, path: Path) -> str:
    """The page as ``write`` would write it to ``path``, glossary link and all; nothing written.

    A builder that fits a page to ``MAX_BYTES`` measures this, then writes.
    """
    target = None
    if isinstance(page, HubPage) or glossary.used(_marked(page)):
        target = glossary.page_path(path)
    return render(page, glossary.href(path, target) if target else "")


def write(page: Page, path: Path) -> Path:
    """Render and write a page, refusing one over 15 MB before anything is written.

    The hub, and a page that marks a glossary term, link the glossary page of their ctx root
    (``glossary.page_path``), which is written first, unless the same page is there already:
    a page never links a file that does not exist (D-048). Outside a ctx root such a page keeps
    its tooltips and links no glossary.

    Args:
        page: The page.
        path: Destination ``.html`` file.

    Returns:
        The written path.

    Raises:
        ValueError: If the page or its glossary exceeds ``MAX_BYTES``.
    """
    target = None
    if isinstance(page, HubPage) or glossary.used(_marked(page)):
        target = glossary.page_path(path)
    html = _sized(rendered(page, path), path)
    if target is not None:
        book = GlossaryPage(text.word("ui", "glossary.title"), glossary.grouped())
        entries = _sized(render(book), target)
        if not (target.is_file() and target.read_bytes() == entries):
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(entries)
    path.write_bytes(html)
    return path


def _sized(html: str, path: Path) -> bytes:
    """A page's bytes, refused over ``MAX_BYTES`` before anything is written."""
    data = html.encode()
    if len(data) > MAX_BYTES:
        raise ValueError(f"{path.name}: {len(data)} bytes, over the {MAX_BYTES}-byte cap")
    return data
