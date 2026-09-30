"""HTML pages: the round trip from one BoxData, offline, size cap, skeleton order, escaping."""

import json
import re
from datetime import date
from html import unescape
from html.parser import HTMLParser
from pathlib import Path

import plotly.graph_objects as go
import plotly.offline
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from pypdf import PdfReader
from synth import (
    counts,
    demo_band_map,
    demo_bar,
    demo_box,
    demo_bubble,
    demo_curves,
    demo_density,
    demo_flow,
    demo_heatmap,
    demo_linked,
    demo_pca,
    demo_strip,
    demo_surface,
    demo_volcano,
)

from ctx.io.naming import figure_files, lint, stem
from ctx.report.builder import (
    MAX_BYTES,
    SECTIONS,
    FigureBlock,
    FigurePage,
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    ResultTable,
    SelfCheck,
    plotly_js,
    render,
    write,
)
from ctx.viz.annotate import TEST_NAME
from ctx.viz.box import BoxData, box_data, interactive
from ctx.viz.save import write_figure
from ctx.viz.style import plotly_templates

ROOT = Path(__file__).resolve().parents[1]
HOSTILE = "</script><b>x"


DEMOS = {
    "box": lambda s: demo_box(s, planted=3, null=5),
    "scatter": demo_pca,
    "volcano": demo_volcano,
    "bar": demo_bar,
    "strip": lambda s: demo_strip(),
    "bubble": lambda s: demo_bubble(),
    "curves": lambda s: demo_curves(),
    "heatmap": lambda s: demo_heatmap(),
    "band_map": lambda s: demo_band_map(),
    "flow": lambda s: demo_flow(),
    "surface": lambda s: demo_surface(),
    "density": lambda s: demo_density(),
    "linked": lambda s: demo_linked(),
}


@pytest.fixture(scope="module", params=list(DEMOS))
def figure(request, tmp_path_factory):
    out = tmp_path_factory.mktemp("report")
    figure_stem = stem("de", f"{request.param}_demo", 1, day=date(2026, 9, 26))
    data = DEMOS[request.param](counts())
    write_figure(data, out, figure_stem, title="Demo", data_releases={})
    return out, figure_stem, data, (out / f"{figure_stem}.html").read_text()


def _tokens(data) -> list[str]:
    """What the caption must state; the renderers wrap it at different widths."""
    if not isinstance(data, BoxData):
        return data.caption.split()
    c = data.caption
    ns = {f"n = {t.n_reference} {c.groups[0]}, {t.n_comparison} {c.groups[1]}" for t in c.tests}
    return [TEST_NAME, f"BH-adjusted over {len(c.tests)} plotted features", *ns]


def _labels(data) -> list[str]:
    """The per-feature p and r a box figure prints over its boxes, in every rendering."""
    if not isinstance(data, BoxData):
        return []
    return [line for t in data.caption.tests for line in data.caption.label(t.feature).split("\n")]


def _plotly_script(html: str) -> str:
    """The figure's drawing script, with each line break read as the space it stands for."""
    (script,) = [
        s for s in re.findall(r"<script[^>]*>(.*?)</script>", html, re.S) if "Plotly.newPlot" in s
    ]
    decoded = re.sub(r"\\u([0-9a-fA-F]{4})", lambda m: chr(int(m.group(1), 16)), script)
    return decoded.replace("<br>", " ")


def _refuse(constant: str) -> None:
    """JSON.parse in a browser rejects NaN and Infinity; so does this test."""
    raise ValueError(f"{constant} is not JSON")


def test_round_trip(figure):
    out, figure_stem, data, html = figure
    assert sorted(p.name for p in out.iterdir()) == sorted(figure_files(figure_stem))
    assert lint(out) == (len(figure_files(figure_stem)), [])

    pdf = " ".join(PdfReader(out / f"{figure_stem}_white.pdf").pages[0].extract_text().split())
    svg = " ".join(re.sub(r"<[^>]+>", " ", (out / f"{figure_stem}_black.svg").read_text()).split())
    # The page prints the statistics once, under the plot (D-011); the plot carries none.
    (caption,) = re.findall(r"<figcaption>(.*?)</figcaption>", html, re.S)
    page = " ".join(unescape(caption).split())
    wanted = {
        "pdf": (pdf, _tokens(data) + _labels(data)),
        "svg": (svg, _tokens(data) + _labels(data)),
        "page caption": (page, _tokens(data)),
        "plot": (_plotly_script(html), _labels(data)),
    }
    for name, (text, tokens) in wanted.items():
        missing = [t for t in tokens if t not in text]
        assert not missing, (name, missing)
    shown = [a.text or "" for a in data.interactive().layout.annotations]
    assert not [a for a in shown if data.caption_text()[:40] in a], "caption drawn in the plot"

    csv = pl.read_csv(out / f"{figure_stem}.csv")
    assert csv.height == data.table().height
    for t in data.caption.tests if isinstance(data, BoxData) else ():
        rows = csv.filter(pl.col("feature") == t.feature)
        assert rows["p_adj"].unique().to_list() == pytest.approx([t.p_adj])
        assert rows["r"].unique().to_list() == pytest.approx([t.r])
    (embedded,) = re.findall(r'id="ctx-data-0">(.*?)</script>', html, re.S)
    rows = pl.DataFrame(json.loads(embedded, parse_constant=_refuse))  # as a browser parses it
    # The page embeds every drawn point, with the columns it calls and names points by.
    columns = list(data.page_table().columns)
    assert rows.columns == columns
    assert_frame_equal(rows, csv.select(columns).fill_nan(None), check_dtypes=False)
    for column in getattr(data, "hover", ()):
        assert f"{column}: %{{customdata" in _plotly_script(html)


def test_html_offline(figure):
    *_, html = figure
    js, raw = plotly_js(), plotly.offline.get_plotlyjs()
    assert html.count(js) == 1 and "\ufffd" not in html
    # Only the raw character changed: each one became a six-character escape. Integers, since
    # a failed == on two 4.8 MB strings sends pytest into a diff that runs for minutes.
    assert (len(js), raw.count("\ufffd")) == (len(raw) + 5 * raw.count("\ufffd"), 1)
    outside = html.replace(js, "")
    assert not re.search(r"<script[^>]+\bsrc\s*=", outside, re.IGNORECASE)
    assert not re.search(r"<link[^>]+\bhref\s*=", outside, re.IGNORECASE)
    assert not re.search(r"url\(\s*['\"]?(https?:)?//", outside, re.IGNORECASE)
    assert len(html.encode()) <= MAX_BYTES
    (templates,) = re.findall(r'id="ctx-templates">(.*?)</script>', html, re.S)
    assert json.loads(templates) == json.loads(json.dumps(plotly_templates(), default=str))
    assert html.index("prefers-color-scheme") < html.index("</head>")


def test_size_guard(tmp_path):
    page = FigurePage("big", FigureBlock(go.Figure(), "x" * (MAX_BYTES + 1), [], "a.csv"))
    with pytest.raises(ValueError, match="cap"):
        write(page, tmp_path / "big.html")
    assert not (tmp_path / "big.html").exists()


def _module(**overrides) -> ModulePage:
    fields = {
        "title": "Setup",
        "question": "Is this machine ready?",
        "input": "This machine.",
        "concept": "A concept.",
        "parameters": (ParameterCard("--root", "Where outputs go.", "cwd", "Moves outputs."),),
        "results": ("One result.",),
        "practice": "In practice.",
        "think": (SelfCheck("Why?", "Because one."), SelfCheck("How?", "Because two.")),
        "mistakes": ("A mistake.",),
        "reproduce": Reproduction(
            "ctx setup", 12.345, 678 * 2**20, 0, "laptop", "B", {"python": "3"}
        ),
    }
    return ModulePage(**{**fields, **overrides})


def test_skeleton_order():
    """The builder's sections carry CONVENTIONS' nine labels verbatim, in its order (D-045)."""
    rule = next(
        line
        for line in (ROOT / "CONVENTIONS.md").read_text().splitlines()
        if line.startswith("- Report sections")
    )
    labels = [re.sub(r"\s*\(.*", "", part).strip(" .") for part in rule.split(":", 1)[1].split("·")]
    assert labels == [label for _, label in SECTIONS]
    html = render(_module())
    positions = [html.index(f'<section id="{section}">') for section, _ in SECTIONS]
    assert positions == sorted(positions)
    with pytest.raises(ValueError, match="think"):
        _module(think=())
    # The frozen M5 and M6 pages have no Input until P2S.4: no section, no contents link.
    bare = render(_module(input=""))
    assert '<section id="input">' not in bare and 'href="#input"' not in bare


class _Elements(HTMLParser):
    """Counts start tags; script content is raw text to the parser, as to a browser."""

    def __init__(self) -> None:
        super().__init__()
        self.tags: list[str] = []

    def handle_starttag(self, tag, attrs):
        self.tags.append(tag)


def _tags(html: str) -> list[str]:
    parser = _Elements()
    parser.feed(html)
    return parser.tags


def test_escaping():
    frame = pl.DataFrame(
        {"feature": [HOSTILE] * 6, "group": ["a", "a", "a", "b", "b", "b"],
         "sample": list("uvwxyz"), "value": [1.0, 2, 3, 4, 5, 6]}
    )  # fmt: skip
    data = box_data(
        frame, columns={c: ("none", c) for c in frame.columns}, feature_label="Gene",
        value_label="Value", groups=("a", "b"),
    )  # fmt: skip
    figure_page = render(
        FigurePage("t", FigureBlock(interactive(data), HOSTILE, frame.to_dicts(), "d.csv"))
    )
    module_page = render(_module(think=(SelfCheck(HOSTILE, HOSTILE),), mistakes=(HOSTILE,)))
    for html in (figure_page, module_page):
        assert "b" not in _tags(html), "hostile text became an element"


def test_quiz_hidden():
    page = _module()
    html = render(page)
    answers = re.findall(r"<details><summary>Answer</summary><p>(.*?)</p></details>", html)
    assert len(answers) == len(page.think) >= 1
    assert answers == [check.answer for check in page.think]


# --- page structure: navigation, tabs, building blocks, copy ------------------------------


def _figure_block(title: str) -> FigureBlock:
    return FigureBlock(go.Figure(), "caption", [], "d.csv", text=FigureText(title=title))


def test_nav_and_tabs_resolve():
    group = Group("Views", (_figure_block("A"), _figure_block("B")), tabs=True)
    html = render(_module(results=("A note.", group, _figure_block("C"))))
    ids = set(re.findall(r'\bid="([^"]+)"', html))
    anchors = re.findall(r'href="#([^"]+)"', html)
    assert anchors and set(anchors) <= ids
    for attribute in ("aria-controls", "aria-labelledby"):
        values = re.findall(rf'{attribute}="([^"]+)"', html)
        assert values and set(values) <= ids, attribute
    # Figures are numbered across groups, in page order.
    assert re.findall(r'<p class="fig-num">Figure (\d+)</p>', html) == ["1", "2", "3"]


def test_building_blocks_refuse_bad_input():
    with pytest.raises(ValueError, match="empty"):
        Group("Nothing", ())
    with pytest.raises(ValueError, match="tabs"):
        Group("One view", (_figure_block("A"),), tabs=True)
    with pytest.raises(ValueError, match="columns"):
        ResultTable("Wide", tuple("abcdefghi"), ())
    with pytest.raises(ValueError, match="wrong length"):
        ResultTable("Ragged", ("a", "b"), (("1",),))
    with pytest.raises(ValueError, match="swatch"):
        ResultTable("Swatch", ("colour",), (("red",),), swatch=0)
    with pytest.raises(ValueError, match="key numbers"):
        _module(key_numbers=tuple(KeyNumber(str(i), "1") for i in range(5)))


def test_templates_use_tokens_only():
    # Colours and lengths come from ctx.viz.style; a literal in a template drifts from it.
    literal = re.compile(r"#[0-9A-Fa-f]{6}\b|rgba?\(|\b\d+(?:\.\d+)?px\b")
    for template in (ROOT / "src" / "ctx" / "report" / "templates").glob("*.j2"):
        found = literal.findall(template.read_text())
        assert not found, (template.name, found[:5])


def test_copy_has_no_dashes_or_arrows():
    # Pages are read by people: no em dash, no arrow, in any string that can reach one.
    for path in (ROOT / "src" / "ctx").rglob("*"):
        if path.suffix in {".py", ".j2", ".toml"}:
            text = path.read_text()
            assert "—" not in text and "→" not in text, path


def test_standalone_page_links_back():
    block = FigureBlock(go.Figure(), "caption", [], "d.csv", text=FigureText(back="m.html"))
    html = render(FigurePage("Alone", block))
    assert 'href="m.html"' in html


def test_hub_prints_counts_whole():
    from ctx.report.hub import _measured  # noqa: PLC0415

    assert _measured(10523.0) == "10,523"
    assert _measured(0.0123456) == "0.01235"
    assert _measured(float("nan")) == "not measured"


def test_text_runs_to_the_column():
    # Only the page is capped: a width cap on the question, paragraphs and captions wrapped
    # sentences that fit on one line at 60 % of the column (reported 2026-09-28).
    css = (ROOT / "src" / "ctx" / "report" / "templates" / "_page.css.j2").read_text()
    capped = [
        selector.strip()
        for selector, body in re.findall(r"([^{};]+)\{([^{}]*)\}", css)
        if re.search(r"max-width:(?!none)", body) and ".page" not in selector
    ]
    assert not capped, capped


def _over(fill: str, ground: str, share: float) -> str:
    """The colour of ``fill`` mixed into ``ground`` at ``share``, as CSS color-mix draws it."""
    pairs = [(int(fill[i : i + 2], 16), int(ground[i : i + 2], 16)) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(share * f + (1 - share) * g):02X}" for f, g in pairs)


def test_text_contrast():
    """Every text colour reads at 4.5:1 or more on every surface it is set on, in both themes.

    ALI's faint grey measured 3.00:1 on its own inset panel, which is why it is not used.
    """
    from ctx.viz.style import CHROME, TEXT_ON, TINT, contrast  # noqa: PLC0415

    share = float(TINT["fill"].rstrip("%")) / 100
    for theme, chrome in CHROME.items():
        for text, grounds in TEXT_ON.items():
            for ground in grounds:
                assert contrast(chrome[text], chrome[ground]) >= 4.5, (theme, text, ground)
        for text, fill in (("text", "good-fill"), ("text", "warn-fill"), ("warn", "warn-fill")):
            tinted = _over(chrome[fill], chrome["page"], share)
            assert contrast(chrome[text], tinted) >= 4.5, (theme, text, fill)


def test_each_kind_of_control_gets_its_own_form(tmp_path):
    """A linked figure gets the brush, a volcano the cutoffs; neither gets the other's form."""
    pages = {}
    for name, data in (("linked", demo_linked()), ("volcano", demo_volcano(counts()))):
        figure_stem = stem("qc", f"{name}_demo", 1, day=date(2026, 9, 30))
        write_figure(data, tmp_path, figure_stem, title="Demo", data_releases={})
        pages[name] = (tmp_path / f"{figure_stem}.html").read_text()
    linked, volcano = pages["linked"], pages["volcano"]
    assert '<form class="brush"' in linked and '<form class="calls"' not in linked
    assert len(re.findall(r'<input type="radio" name="colour"', linked)) == 7
    assert '<p class="brush-summary" aria-live="polite">Nothing painted yet.</p>' in linked
    assert '<form class="calls"' in volcano and '<form class="brush"' not in volcano
    assert '"kind": "brush"' in linked and '"kind": "calls"' in volcano


def test_hub_lists_the_course_modules_in_order():
    """The hub counts M0 to M10 in course order, each titled with its own number (D-044).

    The JupyterLite build writes a site, not a module page, so it is never counted.
    """
    from ctx.io.naming import MODULES  # noqa: PLC0415
    from ctx.report.hub import NOT_MODULES, TITLES  # noqa: PLC0415

    course = [module for module, _ in MODULES if module not in NOT_MODULES]
    assert [TITLES[m].split()[0] for m in course] == [f"M{i}" for i in range(11)]
    assert "lite_build" not in course
