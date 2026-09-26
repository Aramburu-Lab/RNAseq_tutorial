"""Figure style, statistical annotation, and measurements of the written PDF and SVG files."""

import math
import re
from datetime import date
from pathlib import Path
from xml.etree import ElementTree

import matplotlib as mpl
import numpy as np
import polars as pl
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text
from pypdf import PdfReader
from pypdf.generic import ContentStream
from scipy.stats import mannwhitneyu
from synth import counts, demo_box

from ctx.io.naming import stem
from ctx.viz import box
from ctx.viz.annotate import compare_groups
from ctx.viz.save import write_figure
from ctx.viz.style import (
    DATA_COLOURS,
    FONT_CHAIN,
    INK_COLOUR,
    LINE_WIDTHS,
    MAX_CATEGORIES,
    PROJECTOR_UNSAFE,
    ROLES,
    SINGLE_COLUMN_MAX_FEATURES,
    SIZES,
    THEMES,
    WIDTHS_MM,
    apply_matplotlib,
    category_colours,
    contrast,
    css_variables,
)

ROOT = Path(__file__).resolve().parents[1]
STYLE_LITERAL = re.compile(
    r"#[0-9A-Fa-f]{6}\b|rgb\(|\b(?:fontsize|linewidth|size)\s*=\s*\d|figsize\s*=\s*\(\s*\d"
)


# --- style -------------------------------------------------------------------------------


def test_palette_contract():
    figures = next(
        line for line in (ROOT / "CONVENTIONS.md").read_text().splitlines() if "Okabe" in line
    )
    declared = re.findall(r"#[0-9A-F]{6}", figures)
    assert len(declared) == 8
    assert declared[-1] == INK_COLOUR == "#000000"
    assert tuple(declared[:-1]) == DATA_COLOURS
    assert MAX_CATEGORIES == 7
    assert set(ROLES.values()) <= set(DATA_COLOURS)
    assert category_colours(list("abc")) == dict(zip("abc", DATA_COLOURS[:3], strict=True))
    with pytest.raises(ValueError, match="at most 7"):
        category_colours(list("abcdefgh"))
    for theme in THEMES.values():
        css = css_variables(theme)
        assert all(colour in css for colour in DATA_COLOURS)
        assert f"--ink: {theme.ink};" in css


def test_projector_unsafe_is_the_lowest_contrast():
    on_white = {c: contrast(c, "#FFFFFF") for c in DATA_COLOURS}
    assert round(on_white["#F0E442"], 2) == 1.32
    assert min(on_white, key=on_white.get) == "#F0E442"
    assert {"#F0E442"} == PROJECTOR_UNSAFE


def test_no_style_literals():
    files = [p for p in (ROOT / "src" / "ctx").rglob("*.py") if p.name != "style.py"]
    assert len(files) > 5
    hits = [
        f"{p.relative_to(ROOT)}:{i}"
        for p in files
        for i, line in enumerate(p.read_text().splitlines(), 1)
        if STYLE_LITERAL.search(line)
    ]
    assert not hits


def test_context_then_sizes():
    with mpl.rc_context():
        apply_matplotlib(THEMES["white"])
        rc = mpl.rcParams
        assert rc["font.size"] == SIZES["base"]
        assert rc["axes.labelsize"] == SIZES["axis_label"]
        assert rc["xtick.labelsize"] == rc["ytick.labelsize"] == SIZES["tick"]
        assert rc["axes.linewidth"] == LINE_WIDTHS["axes"]
        assert rc["lines.linewidth"] == LINE_WIDTHS["line"]
        assert rc["patch.linewidth"] == LINE_WIDTHS["box"]


def test_font_fallback(tmp_path):
    with mpl.rc_context():
        apply_matplotlib(THEMES["white"], fonts=("NoSuchFontXYZ", "DejaVu Sans"))
        fig = Figure()
        fig.text(0.5, 0.5, "fallback")
        fig.savefig(tmp_path / "f.pdf")
        fig.savefig(tmp_path / "f.svg")
    assert _pdf_fonts(tmp_path / "f.pdf") == {"DejaVuSans"}
    assert "'DejaVu Sans'" in (tmp_path / "f.svg").read_text()


# --- statistics --------------------------------------------------------------------------


def _long(values: dict[str, tuple[list[float], list[float]]]) -> pl.DataFrame:
    rows = [
        (feature, group, f"s{i}", v)
        for feature, (ref, cmp) in values.items()
        for group, vs in (("ref", ref), ("cmp", cmp))
        for i, v in enumerate(vs)
    ]
    return pl.DataFrame(rows, schema=["feature", "group", "sample", "value"], orient="row")


def _caption(frame: pl.DataFrame):
    return compare_groups(
        frame, feature="feature", group="group", value="value", groups=("ref", "cmp")
    )


def _hand_bh(p: list[float]) -> list[float]:
    m, order = len(p), sorted(range(len(p)), key=lambda i: p[i])
    adjusted, running = [0.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, p[i] * m / rank)
        adjusted[i] = running
    return adjusted


FOUR = {
    "mid": ([1, 2, 3, 6], [4, 5, 7, 8]),
    "separated": ([1, 2, 3, 4], [5, 6, 7, 8]),
    "same": ([1, 2, 3, 4], [1, 2, 3, 4]),
    "reversed": ([5, 6, 7, 8], [1, 2, 3, 4]),
}


def test_bh_input_order():
    caption = _caption(_long(FOUR))
    raw = [t.p for t in caption.tests]
    assert [t.feature for t in caption.tests] == list(FOUR)
    assert raw[1] == pytest.approx(2 / math.comb(8, 4))
    assert raw != sorted(raw)
    assert [t.p_adj for t in caption.tests] == pytest.approx(_hand_bh(raw))
    assert max(t.p_adj for t in caption.tests) <= 1


def test_bh_after_subset():
    subset = {k: FOUR[k] for k in ("mid", "separated")}
    full = {t.feature: t.p_adj for t in _caption(_long(FOUR)).tests}
    part = _caption(_long(subset))
    assert [t.p_adj for t in part.tests] == pytest.approx(_hand_bh([t.p for t in part.tests]))
    assert [t.p_adj for t in part.tests] != pytest.approx([full[k] for k in subset])


def test_nan_dropped_before_n():
    frame = _long({"g": ([1.0, 2.0, 3.0, None], [4.0, 5.0, float("nan"), 6.0])})
    (test,) = _caption(frame).tests
    assert (test.n_reference, test.n_comparison) == (3, 3)
    assert test.p == pytest.approx(mannwhitneyu([4, 5, 6], [1, 2, 3]).pvalue)


def test_rank_biserial_truth():
    caption = _caption(_long({k: FOUR[k] for k in ("separated", "same", "reversed")}))
    assert [t.r for t in caption.tests] == [1.0, 0.0, -1.0]


def test_planted_recovered():
    synthetic = counts()
    data = demo_box(synthetic, planted=3, null=5)
    n_ref, n_cmp = 6, 6
    smallest = 2 / math.comb(n_ref + n_cmp, n_ref)
    assert smallest * len(data.caption.tests) < 0.05, "no gene could reach BH < 0.05 at this n"
    planted = dict(synthetic.truth.select("gene_id", "lfc_condition_type_a").iter_rows())
    for t in data.caption.tests:
        assert (t.n_reference, t.n_comparison) == (n_ref, n_cmp)
        if planted[t.feature] != 0:
            assert t.p_adj < 0.05
            assert np.sign(t.r) == np.sign(planted[t.feature])
        else:
            assert t.p_adj >= 0.05


# --- the written files -------------------------------------------------------------------


def _pdf_fonts(path: Path) -> set[str]:
    page = PdfReader(path).pages[0]
    fonts = [f.get_object() for f in page["/Resources"]["/Font"].values()]
    assert all(f["/Subtype"] != "/Type3" for f in fonts)
    return {str(f["/BaseFont"]).split("+")[-1].lstrip("/") for f in fonts}


COLOUR_OPS = {b"RG": "S", b"G": "S", b"rg": "F", b"g": "F"}
PAINT = {
    b"S": "s",
    b"s": "s",
    b"f": "f",
    b"F": "f",
    b"f*": "f",
    b"B": "sf",
    b"B*": "sf",
    b"b": "sf",
    b"b*": "sf",
}


def _painted(path: Path) -> tuple[set[tuple[float, ...]], set[float], set[float], bool]:
    """Colours actually painted, Tf sizes, line widths, and whether a fill covers the page."""
    reader = PdfReader(path)
    page = reader.pages[0]
    width, height = float(page.mediabox.width), float(page.mediabox.height)
    states = page["/Resources"].get("/ExtGState", {})
    state = {"S": (0.0, 0.0, 0.0), "F": (0.0, 0.0, 0.0), "w": 1.0, "CA": 1.0, "ca": 1.0}
    stack, path_points = [], []
    colours, sizes, widths, page_fill = set(), set(), set(), False
    for args, op in ContentStream(page.get_contents(), reader).operations:
        values = [float(a) for a in args if isinstance(a, int | float)]
        if op in (b"q", b"Q"):
            state = stack.pop() if op == b"Q" else (stack.append(dict(state)) or state)
        elif op in COLOUR_OPS:
            state[COLOUR_OPS[op]] = tuple(values * 3 if len(values) == 1 else values)
        elif op == b"w":
            state["w"] = values[0]
            widths.add(values[0])
        elif op == b"gs":
            gs = states[args[0]].get_object()
            state["CA"], state["ca"] = float(gs.get("/CA", 1)), float(gs.get("/ca", 1))
        elif op == b"Tf":
            sizes.add(float(args[1]))
        elif op in (b"m", b"l", b"c", b"re"):
            path_points += _points(op, values)
        elif op in (b"TJ", b"Tj", b"Do"):
            colours.add(state["F"])
        elif op in PAINT or op == b"n":
            painted, covers = _paint(PAINT.get(op, ""), state, path_points, width * height)
            colours |= painted
            page_fill |= covers
            path_points = []
    return colours, sizes, widths, page_fill


def _points(op: bytes, values: list[float]) -> list[tuple[float, float]]:
    """Corner points a path operator adds; a rectangle adds its two opposite corners."""
    if op == b"re":
        x, y, w, h = values
        return [(x, y), (x + w, y + h)]
    return list(zip(values[::2], values[1::2], strict=True))


def _paint(mode: str, state: dict, points: list, page_area: float) -> tuple[set, bool]:
    """Colours a paint operator really lays down, and whether its fill covers the page."""
    painted, covers = set(), False
    if "s" in mode and state["w"] > 0 and state["CA"] > 0:
        painted.add(state["S"])
    if "f" in mode and state["ca"] > 0:
        painted.add(state["F"])
        xs, ys = [p[0] for p in points], [p[1] for p in points]
        covers = bool(xs) and (max(xs) - min(xs)) * (max(ys) - min(ys)) >= 0.99 * page_area
    return painted, covers


def _rgb(colour: str) -> tuple[float, ...]:
    return tuple(int(colour[i : i + 2], 16) / 255 for i in (1, 3, 5))


def _same(a: tuple[float, ...], b: tuple[float, ...]) -> bool:
    return all(abs(x - y) <= 1 / 255 for x, y in zip(a, b, strict=True))


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    out = tmp_path_factory.mktemp("figure")
    figure_stem = stem("de", "box_demo", 1, day=date(2026, 9, 26))
    data = demo_box(counts(), planted=3, null=5)
    files = write_figure(data, out, figure_stem, title="Demo", data_releases={})
    return out, figure_stem, data, files


def test_caption_wrap_keeps_words_whole():
    caption = demo_box(counts(), planted=1, null=1).caption
    for width in (30, 36, 40, 73):
        assert "BH-adjusted" in caption.text(width).replace("\n", " ")


def test_static_files_measured(written):
    out, figure_stem, data, _ = written
    caption = re.sub(r"\s+", " ", data.caption.text())
    for theme in THEMES.values():
        pdf = out / f"{figure_stem}_{theme.name}.pdf"
        assert _pdf_fonts(pdf) <= {f.replace(" ", "") for f in FONT_CHAIN} | {"ArialMT"}
        colours, sizes, widths, page_fill = _painted(pdf)
        allowed = [_rgb(c) for c in (*DATA_COLOURS, theme.ink, theme.grid, theme.reference)]
        assert all(any(_same(c, a) for a in allowed) for c in colours), colours
        assert any(_same(c, _rgb(theme.ink)) for c in colours)
        assert sizes and sizes <= set(SIZES.values())
        assert widths and widths <= set(LINE_WIDTHS.values())
        assert not page_fill, "the canvas must stay transparent"
        page = PdfReader(pdf).pages[0]
        assert page["/Resources"].get("/XObject") is None or all(
            x.get_object()["/Subtype"] != "/Image" for x in page["/Resources"]["/XObject"].values()
        )
        n_features = len(data.caption.tests)
        column = "single" if n_features <= SINGLE_COLUMN_MAX_FEATURES else "double"
        assert float(page.mediabox.width) == pytest.approx(WIDTHS_MM[column] / 25.4 * 72, abs=0.5)
        text = re.sub(r"\s+", " ", page.extract_text())
        assert caption in text
        assert f"{data.value_label} ({data.unit})" in text

        svg = ElementTree.parse(out / f"{figure_stem}_{theme.name}.svg")
        families = {
            name.strip().strip("'\"")
            for element in svg.iter()
            for match in re.findall(r"font-family:\s*([^;]+)", element.get("style", ""))
            for name in match.split(",")
        }
        assert families and families <= {*FONT_CHAIN, "sans-serif"}
        assert not [e for e in svg.iter() if e.tag.endswith("image")]


@pytest.mark.parametrize(("planted", "null"), [(1, 1), (3, 5)])
def test_no_text_overlap(planted, null):
    data = demo_box(counts(), planted=planted, null=null)
    for theme in THEMES.values():
        fig = box.static(data, theme)
        renderer = FigureCanvasAgg(fig).get_renderer()
        fig.draw(renderer)
        texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
        boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in texts]
        assert len(boxes) > 2 * (planted + null)
        page = fig.bbox
        outside = [name for name, b in boxes if b.x0 < page.x0 or b.x1 > page.x1 or b.y0 < page.y0]
        assert not outside, f"text off the page: {outside}"
        clashes = [
            (a, b)
            for i, (a, box_a) in enumerate(boxes)
            for b, box_b in boxes[i + 1 :]
            if box_a.overlaps(box_b)
        ]
        assert not clashes, clashes[:5]
