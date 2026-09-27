"""HTML pages: the round trip from one BoxData, offline, size cap, skeleton order, escaping."""

import json
import re
from datetime import date
from html.parser import HTMLParser
from pathlib import Path

import plotly.graph_objects as go
import plotly.offline
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from pypdf import PdfReader
from synth import counts, demo_bar, demo_box, demo_pca, demo_volcano

from ctx.io.naming import figure_files, lint, stem
from ctx.report.builder import (
    MAX_BYTES,
    SECTIONS,
    FigureBlock,
    FigurePage,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
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
}


@pytest.fixture(scope="module", params=list(DEMOS))
def figure(request, tmp_path_factory):
    out = tmp_path_factory.mktemp("report")
    figure_stem = stem("de", f"{request.param}_demo", 1, day=date(2026, 9, 26))
    data = DEMOS[request.param](counts())
    write_figure(data, out, figure_stem, title="Demo", data_releases={})
    return out, figure_stem, data, (out / f"{figure_stem}.html").read_text()


def _tokens(data) -> list[str]:
    if not isinstance(data, BoxData):
        # Every word of the caption; the renderers wrap it at different widths.
        return data.caption.split()
    c = data.caption
    ns = {f"n = {t.n_reference} {c.groups[0]}, {t.n_comparison} {c.groups[1]}" for t in c.tests}
    labels = [line for t in c.tests for line in c.label(t.feature).split("\n")]
    return [TEST_NAME, f"BH-adjusted over {len(c.tests)} plotted features", *ns, *labels]


def _plotly_script(html: str) -> str:
    """The figure's drawing script, with each line break read as the space it stands for."""
    (script,) = [
        s for s in re.findall(r"<script[^>]*>(.*?)</script>", html, re.S) if "Plotly.newPlot" in s
    ]
    decoded = re.sub(r"\\u([0-9a-fA-F]{4})", lambda m: chr(int(m.group(1), 16)), script)
    return decoded.replace("<br>", " ")


def test_round_trip(figure):
    out, figure_stem, data, html = figure
    assert sorted(p.name for p in out.iterdir()) == sorted(figure_files(figure_stem))
    assert lint(out) == (len(figure_files(figure_stem)), [])

    pdf = " ".join(PdfReader(out / f"{figure_stem}_white.pdf").pages[0].extract_text().split())
    svg = " ".join(re.sub(r"<[^>]+>", " ", (out / f"{figure_stem}_black.svg").read_text()).split())
    for name, text in {"pdf": pdf, "svg": svg, "plotly": _plotly_script(html)}.items():
        missing = [t for t in _tokens(data) if t not in text]
        assert not missing, (name, missing)

    csv = pl.read_csv(out / f"{figure_stem}.csv")
    assert csv.height == data.table().height
    for t in data.caption.tests if isinstance(data, BoxData) else ():
        rows = csv.filter(pl.col("feature") == t.feature)
        assert rows["p_adj"].unique().to_list() == pytest.approx([t.p_adj])
        assert rows["r"].unique().to_list() == pytest.approx([t.r])
    (embedded,) = re.findall(r'id="ctx-data-0">(.*?)</script>', html, re.S)
    assert_frame_equal(pl.DataFrame(json.loads(embedded)), csv, check_dtypes=False)
    for column in data.hover:
        assert f"{column}: %{{customdata" in _plotly_script(html)


def test_html_offline(figure):
    *_, html = figure
    js = plotly.offline.get_plotlyjs()
    assert html.count(js) == 1
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
    rule = next(
        line
        for line in (ROOT / "CONVENTIONS.md").read_text().splitlines()
        if "Report skeleton" in line
    )
    parts = rule.split(":", 1)[1].split("→")
    expected = [
        ("big question", "question"), ("concept", "concept"), ("parameter cards", "parameters"),
        ("annotated results", "results"), ("in practice", "practice"), ("think-about-it", "think"),
        ("common mistakes", "mistakes"), ("reproduction", "reproduce"),
    ]  # fmt: skip
    assert len(parts) == len(SECTIONS) == len(expected)
    assert all(keyword in part for (keyword, _), part in zip(expected, parts, strict=True))
    assert [section for section, _ in SECTIONS] == [section for _, section in expected]
    html = render(_module())
    positions = [html.index(f'<section id="{section}">') for section, _ in SECTIONS]
    assert positions == sorted(positions)
    with pytest.raises(ValueError, match="think"):
        _module(think=())


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
