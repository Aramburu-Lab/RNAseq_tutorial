"""ED.7: the core edition's lite pages (D-111), asserted by structure, never by words.

A lite page keeps each figure card (title, lead, reading guide, figure, caption, controls and
the PDF, SVG and CSV links) and every results table, and drops the sections, the Reproduce
panels, the group leads, the check tables and the notes: less text, not less to do (Rytis
2026-10-22). A full page is unchanged. Mutations planted and caught here: a section leaking
into lite (``sections`` in ``builder.render`` left whole), a reading guide or a control dropped
(``_figure_as`` emptying ``text`` or ``controls``) and a results table dropped (``_item_as``
dropping every table).
"""

from dataclasses import replace
from html.parser import HTMLParser

import pytest
from test_report import _module
from test_viz import filtered_bars

from ctx import editions
from ctx.process.blocks import page_table
from ctx.report import text
from ctx.report.builder import (
    FULL,
    SECTIONS,
    Code,
    CourseLink,
    CoursePage,
    CourseRow,
    FigurePage,
    FigureText,
    Group,
    HubLink,
    HubPage,
    HubRule,
    KeyNumber,
    Look,
    Panel,
    ParameterCard,
    ResultTable,
    look,
    render,
    results_as,
    unpacked,
)
from ctx.viz.save import figure_block

RESULTS = "Calls per contrast"
CHECK = "Flagged runs"
REPRODUCE = {"script_head": "", "summary": "Reproduce", "lead": "Lead.", "shared": "ctx x",
             "note": "Note.", "upstream_title": "Upstream", "upstream": "ctx y",
             "tabs": [{"key": "docker", "label": "Docker", "head": "docker run", "name": "a.sh",
                       "download": "Download"}]}  # fmt: skip


class Tags(HTMLParser):
    """Every start tag of a page, as (tag, attributes), and each table caption's text."""

    def __init__(self, html: str) -> None:
        super().__init__()
        self.tags: list[tuple[str, dict[str, str]]] = []
        self.captions: list[str] = []
        self._caption = False
        self.feed(unpacked(html))

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.tags.append((tag, {k: v or "" for k, v in attrs}))
        self._caption = tag == "caption"

    def handle_data(self, data: str) -> None:
        if self._caption:
            self.captions.append(data)
            self._caption = False

    def having(self, tag: str = "", css: str = "", **attrs: str) -> list[dict[str, str]]:
        """The tags of that name (any, when empty) with that class and those attributes."""
        return [a for t, a in self.tags
                if (not tag or t == tag) and (not css or css in a.get("class", "").split())
                and all(a.get(k.replace("_", "-")) == v for k, v in attrs.items())]  # fmt: skip


def _figure(number: int, controls: bool = True):
    """A filtered bar figure with its lead, reading guide, controls and Reproduce panel."""
    data = filtered_bars()
    words = FigureText(title=f"Depth {number}", lead="A lead.", reading=("Each bar is a run.",))
    block = figure_block(data, f"2026-10-22_qc_depth{number}_v1", words)
    return replace(block, controls=data.controls() if controls else None, reproduce=REPRODUCE)


def every_part() -> object:
    """A module page holding one of every part a lite page keeps or drops."""
    panel = Panel("genelist", {"kind": "genelist", "words": {}, "contrasts": [], "featured": ""})
    calls = Group("Calls", (
        _figure(1),
        ResultTable(RESULTS, ("Contrast", "Up"), (("a", "1"),), results=True),
        ResultTable(CHECK, ("Run", "Reason"), (("r1", "low"),)),
        Code("Run it", "ctx de"),
        "A note.",
        panel,
    ), lead="A group lead.")  # fmt: skip
    tabs = Group("Views", (_figure(2), _figure(3)), tabs=True)
    return _module(
        results=(calls, tabs, ResultTable("Versions", ("Tool",), (("x",),))),
        key_numbers=(KeyNumber("Runs", "6"),),
        terms=(("PCA", "A projection."),),
    )


@pytest.fixture
def core(monkeypatch):
    monkeypatch.setenv(editions.EDITION_ENV, "core")


@pytest.fixture
def full(monkeypatch):
    monkeypatch.delenv(editions.EDITION_ENV, raising=False)


def test_the_editions_say_how_their_pages_look(full):
    """core is lite with its controls, full is every section with them; each is one line."""
    assert look() == FULL == Look(lite=False, controls=True)
    core = editions.named("core")
    assert (core.page, core.controls) == ("lite", True)
    assert (editions.active().page, editions.active().controls) == ("full", True)


@pytest.mark.parametrize(
    ("line", "problem"),
    [('page = "short"', "edition.page"), ('controls = "no"', "edition.controls")],
)
def test_an_edition_with_a_bad_look_is_refused(tmp_path, line, problem):
    path = tmp_path / "odd.toml"
    path.write_text(f'[edition]\nname = "odd"\nlabel = "Odd"\nnotes = "Odd."\nmodules = "all"\n'
                    f'datasets = "all"\n{line}\n')  # fmt: skip
    with pytest.raises(editions.EditionError, match=problem):
        editions.load(path)


def test_a_core_page_keeps_only_its_figure_and_results_cards(core):
    """Catches a section, a check table, a note or a Reproduce panel leaking into lite, and a
    figure card, its lead, reading guide or controls, a caption, an export link or a results
    table dropped from it."""
    page = Tags(render(every_part()))
    assert [a["id"] for a in page.having("section")] == ["results"]
    for gone in ("toc", "keys", "group-lead", "repro", "terms"):
        assert page.having(css=gone) == [], gone
    built = every_part()
    calls, *rest = built.results
    bare = replace(calls, items=tuple(i for i in calls.items if not isinstance(i, str)))
    # The controls carry notes of their own, so the page note is shown gone by the page alone.
    assert render(built) == render(replace(built, results=(bare, *rest)))
    # The gene-list panel's own command line is a code figure too, hidden until a list is set.
    assert [a for a in page.having(css="code") if "list-command" not in a["class"]] == []
    assert len(page.having(css="reading")) == 3 and len(page.having(css="fig-lead")) == 3
    assert len(page.having("form")) >= 3 and page.having(role="tablist")
    assert page.having(css="genelist-panel")
    assert [a for a in page.having("script") if a.get("id", "").startswith("ctx-controls")]
    # The gene-list panel draws a figure of its own, with no id and a captioned class.
    figures = [a for a in page.having("figure", "fig") if "id" in a]
    assert [a["id"] for a in figures] == ["figure-1", "figure-2", "figure-3"]
    bare = [a for a in page.having("figcaption") if not a.get("class")]
    assert len(bare) == len(figures) + len(page.having("figure", "list-command"))
    assert len(page.having(css="fig-title")) == len(page.having("figure", "fig"))
    hrefs = [a["href"] for a in page.having("a") if "href" in a]
    for n in (1, 2, 3):
        for end in ("_white.pdf", "_black.pdf", "_white.svg", "_black.svg"):
            assert f"2026-10-22_qc_depth{n}_v1{end}" in hrefs
    assert sum(h.endswith(".csv") for h in hrefs) == 3
    assert page.captions == [RESULTS]


def test_the_same_page_under_full_keeps_every_part(full):
    """Catches the lite rule leaking into full: the nine sections, the forms, the Reproduce
    panels, the tabs and both tables stay, and full's results pass through untouched."""
    built = every_part()
    assert results_as(built.results, FULL) == built.results
    page = Tags(render(built))
    assert [a["id"] for a in page.having("section")] == [s for s, _ in SECTIONS]
    assert len(page.having("form")) >= 3 and len(page.having("details", "repro")) == 3
    assert page.having(role="tablist") and page.having(css="genelist-panel")
    assert page.captions[:2] == [RESULTS, CHECK]


def test_controls_go_with_one_line(tmp_path, monkeypatch):
    """A lite edition with ``controls = false`` drops the cards' forms and tabs and keeps their
    reading guides. Catches controls read from anywhere but the edition."""
    folder = tmp_path / "editions"
    folder.mkdir()
    (folder / "litectl.toml").write_text(
        '[edition]\nname = "litectl"\nlabel = "L"\nnotes = "L."\nmodules = "all"\n'
        'datasets = "all"\npage = "lite"\ncontrols = false\n')  # fmt: skip
    monkeypatch.setenv(editions.FOLDER_ENV, str(folder))
    monkeypatch.setenv(editions.EDITION_ENV, "litectl")
    page = Tags(render(every_part()))
    assert [a["id"] for a in page.having("section")] == ["results"]
    assert page.having("form") == [] and page.having(role="tablist") == []
    assert page.having(css="genelist-panel") == [] and page.having(css="repro") == []
    assert len(page.having(css="reading")) == 3


def test_a_figure_opened_alone_is_its_card_too(core):
    page = Tags(render(FigurePage("Depth", _figure(1))))
    assert page.having("form") and page.having(css="reading") and page.having(css="fig-lead")
    assert page.having(css="repro") == []
    assert len(page.having("figcaption")) == 1


def test_a_lite_hub_and_home_page_are_their_links(core):
    """The hub keeps its module list; START.html its hubs and modules. Catches the intro, key
    numbers, rules, cards or notes leaking into either."""
    hub = HubPage("Hub", "Intro.", (ParameterCard("--x", "A card.", "1", "Moves it."),),
                  (HubLink("M1 Gather", "M1 Gather", "../01_gather/a.html", "2026-10-22"),),
                  ("M2 QC",), (HubRule("r1", "A rule.", "pass", (("a > 1", "2", "pass"),)),),
                  (KeyNumber("Modules", "1 of 9"),))  # fmt: skip
    page = Tags(render(hub))
    assert [a["href"] for a in page.having("a", "module-name")] == ["../01_gather/a.html"]
    for gone in ("lede", "keys", "missing", "params", "rules"):
        assert page.having(css=gone) == [], gone
    home = CoursePage("Course", (CourseRow("M1", "Gather", (CourseLink("step", "a.html", True),)),),
                      (CourseLink("Hub", "hub.html", True),), (KeyNumber("Pages", "2"),),
                      "glossary.html", "guide.html")  # fmt: skip
    page = Tags(render(home))
    assert {a["id"] for a in page.having("section")} == {"course", "hubs", "modules"}
    assert page.having(css="lede") == [] and page.having(css="keys") == []
    hrefs = [a["href"] for a in page.having("a") if "href" in a]
    assert {"a.html", "hub.html", "glossary.html", "guide.html"} <= set(hrefs)


def test_the_guide_and_the_teacher_page_stay_whole(core):
    """Instructions are not module reports: a lite edition shows them whole."""
    page = Tags(render(replace(every_part(), instructions=True)))
    assert [a["id"] for a in page.having("section")] == [s for s, _ in SECTIONS]


def test_m3_declares_which_of_its_tables_are_results():
    """The flag is set where M3 declares a table: its measurements, never its records."""
    keys = ("alignment", "strandedness", "concordance", "reference", "folders")
    width = {k: len(text.words("process", f"tables.{k}.page.columns")) for k in keys}
    kept = {k: page_table(k, [("a",) * width[k]], ()).results for k in keys}
    assert kept == {"alignment": True, "strandedness": True, "concordance": True,
                    "reference": False, "folders": False}  # fmt: skip
