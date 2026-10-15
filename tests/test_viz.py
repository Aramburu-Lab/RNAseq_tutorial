"""Figure style, statistical annotation, and measurements of the written PDF and SVG files."""

import dataclasses
import math
import re
from datetime import date
from pathlib import Path
from xml.etree import ElementTree

import matplotlib as mpl
import numpy as np
import polars as pl
import pytest
from matplotlib.figure import Figure
from pypdf import PdfReader
from pypdf.generic import ContentStream
from scipy.stats import mannwhitneyu
from synth import (
    counts,
    demo_band_map,
    demo_bar,
    demo_bars_at_scale,
    demo_box,
    demo_bubble,
    demo_curves,
    demo_density,
    demo_flow,
    demo_heatmap,
    demo_linked,
    demo_ma,
    demo_pca,
    demo_strip,
    demo_surface,
    demo_volcano,
)

from ctx.de.stats import classify, p_line
from ctx.io.naming import stem
from ctx.report.builder import FigureText
from ctx.viz import box
from ctx.viz.annotate import compare_groups, permanova
from ctx.viz.axes import axis_title
from ctx.viz.bar import BarData
from ctx.viz.contour import SurfaceData
from ctx.viz.curve import CurveData
from ctx.viz.density import DensityData
from ctx.viz.flow import FlowData
from ctx.viz.heatmap import HeatmapData
from ctx.viz.linked import LinkedData
from ctx.viz.save import drawn_texts, figure_block, layout_problems, write_figure
from ctx.viz.scatter import THIN_MAX, Calls, Line, ScatterData, bubble_pt
from ctx.viz.strip import StripData
from ctx.viz.style import (
    AXIS_TITLE_SHARE,
    CIVIDIS,
    DATA_COLOURS,
    FONT_CHAIN,
    INK_COLOUR,
    LINE_WIDTHS,
    MAX_CATEGORIES,
    PAGE_PLOT_PX,
    PLOTLY_HEIGHT_PX,
    PROJECTOR_UNSAFE,
    ROLES,
    SIZES,
    THEMES,
    WIDTHS_MM,
    apply_matplotlib,
    category_colours,
    contrast,
    css_variables,
)

ROOT = Path(__file__).resolve().parents[1]
# A plotly font or marker as a dict literal, {"size": 12}, passed the keyword form alone
# (2026-10-01); a contour's "size" is a step in data units, not a style.
STYLE_LITERAL = re.compile(
    r"#[0-9A-Fa-f]{6}\b|rgb\(|\b(?:fontsize|linewidth|size)\s*=\s*\d|figsize\s*=\s*\(\s*\d"
    r"""|(?:font|marker)["']?\s*[:=]\s*\{[^}]*["']size["']\s*:\s*\d"""
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


def test_style_literal_reads_dict_sizes():
    typed = ('"tickfont": {"size": 12}', 'legend={"font": {"size": 11}}',
             'marker={"symbol": "square", "size": 6}', "title_font=dict(size=12)")  # fmt: skip
    for line in typed:
        assert STYLE_LITERAL.search(line), line
    assert not STYLE_LITERAL.search('contours={"start": level, "end": level, "size": 1}')


def test_rc_from_style_alone():
    """Every size comes from style.py: no seaborn context resizes ticks or markers after it."""
    from ctx.viz.style import MARKER_SIZE, TICK_PT  # noqa: PLC0415

    with mpl.rc_context():
        apply_matplotlib(THEMES["white"])
        rc = mpl.rcParams
        assert rc["xtick.major.size"] == rc["ytick.major.size"] == TICK_PT["major"]
        assert rc["xtick.minor.width"] == rc["ytick.minor.width"] == LINE_WIDTHS["axes"]
        assert rc["lines.markersize"] == MARKER_SIZE
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


def filtered_bars(kind: str = "runs") -> BarData:
    """Bars the page filters, over two reference lines and a median line.

    ``runs``: one bar per run, as the M1 depth figure draws them. ``groups``: two batches
    stacked in each group, so a bar's length is a total. ``grouped``: the same side by side.
    """
    rng = np.random.default_rng(11)
    types, conditions = ["NDN", "Monocyte", "LDN"], ["steady", "pdac"]
    if kind == "runs":
        n = 40
        frame = pl.DataFrame({"bar": [f"r{i:02d}" for i in range(n)],
                              "batch": ["short"] * 22 + ["long"] * 18,
                              "reads": rng.uniform(5.0, 60.0, n),
                              "cell_type": rng.choice(types, n).tolist(),
                              "condition": rng.choice(conditions, n).tolist()})  # fmt: skip
        unit = "million reads"
    else:
        pairs = [(t, c) for t in types for c in conditions]
        frame = pl.DataFrame({"bar": [f"{t} {c}" for t, c in pairs] * 2,
                              "batch": ["short"] * 6 + ["long"] * 6,
                              "reads": rng.integers(1, 12, 12).astype(float),
                              "cell_type": [t for t, _ in pairs] * 2,
                              "condition": [c for _, c in pairs] * 2})  # fmt: skip
        unit = "samples"
    columns = {"bar": ("none", "Bar."), "batch": ("none", "Batch."), "reads": (unit, "Length."),
               "cell_type": ("none", "Cell type."),
               "condition": ("none", "Condition.")}  # fmt: skip
    return BarData(
        frame=frame, columns=columns, category="bar", value="reads", category_label="Run",
        value_label="Reads", series="batch", colours=category_colours(["short", "long"]),
        caption="Bars the page filters; no test.", stacked=kind != "grouped",
        hover=("cell_type",), width="double", lines=((30.0, "thirty"), (20.0, "twenty")),
        median="median of the bars shown", filters=("cell_type", "condition"),
    )  # fmt: skip


def _contracts():
    synthetic = counts()
    return {
        "box": demo_box(synthetic, planted=3, null=5),
        "volcano_single": demo_volcano(synthetic),
        "volcano_double": demo_volcano(synthetic, width="double"),
        "pca": demo_pca(synthetic),
        "ma_log_x": demo_ma(synthetic),
        "bar": demo_bar(synthetic),
        "bar_horizontal": demo_bar(synthetic, horizontal=True),
        "composition_20_groups": demo_bars_at_scale(horizontal=True),
        "library_sizes_210_runs": demo_bars_at_scale(horizontal=False),
        "bar_filtered_lines": filtered_bars(),
        "strip": demo_strip(),
        "strip_log_double": demo_strip(log=True, width="double"),
        "bubble": demo_bubble(),
        "curves": demo_curves(),
        "heatmap": demo_heatmap(),
        "band_map": demo_band_map(),
        "flow": demo_flow(),
        "surface": demo_surface(),
        "density": demo_density(),
        "linked": demo_linked(),
    }


KINDS = tuple(_contracts())


def static_caption(data) -> str:
    """The caption the static figure prints."""
    return data.caption.text() if isinstance(data, box.BoxData) else data.caption


def axis_titles(data) -> list[tuple[str, str]]:
    """(label, unit) of every value axis the static figure must print."""

    def xy():
        return [(data.x_label, data.columns[data.x][0]), (data.y_label, data.columns[data.y][0])]

    def value():
        return [(data.value_label, data.columns[data.value][0])]

    rules = {
        box.BoxData: lambda: [(data.value_label, data.unit)],
        ScatterData: xy,
        CurveData: xy,
        DensityData: xy,
        SurfaceData: lambda: [*xy(), (data.z_label, data.columns[data.z][0])],
        StripData: value,
        HeatmapData: lambda: [(data.column_label, "none")],
        FlowData: list,
        LinkedData: list,
    }
    return rules.get(type(data), value)()


def _unit_words(unit: str) -> list[str]:
    """The words of a unit, singular, as a reader would look for them in a title."""
    return [w.removesuffix("s") for w in re.findall(r"[a-z0-9]+", unit.lower()) if unit != "none"]


def test_axis_title_rule():
    # Literal cases, so the rule is pinned independently of the function under test.
    assert axis_title("Genes called", "genes") == "Genes called"
    assert axis_title("Samples", "samples") == "Samples"
    assert axis_title("Reads", "million reads") == "Reads (million reads)"
    assert axis_title("PC1", "VST units") == "PC1 (VST units)"
    assert axis_title("Dispersion", "dimensionless") == "Dispersion (dimensionless)"
    assert axis_title("Gene set", "none") == "Gene set"
    # A unit that is only a symbol was taken as named by any label, and dropped.
    assert axis_title("Reads kept", "%") == "Reads kept (%)"


def column_width(data) -> str:
    """The column width the figure should take: the double column, for every figure."""
    return "double" if isinstance(data, box.BoxData) else data.width


def test_every_figure_takes_the_double_column():
    """Rytis, 2026-10-01: at 85 mm every figure stood taller than wide, since a height is 85 mm
    times the type's aspect. Mutation caught: a contract default or a call site back at "single"."""
    root = Path(__file__).resolve().parents[1]
    sources = [*(root / "src").rglob("*.py"), *(root / "scripts").rglob("*.py")]
    named = [f"{p.relative_to(root)}:{n}" for p in sources if p.name != "style.py"
             for n, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1)
             if '"single"' in line]  # fmt: skip
    assert named == [], named
    contracts = [BarData, SurfaceData, CurveData, DensityData, FlowData, HeatmapData, LinkedData,
                 ScatterData, StripData]  # fmt: skip
    widths = {c.__name__: next(f.default for f in dataclasses.fields(c) if f.name == "width")
              for c in contracts}  # fmt: skip
    assert widths == {c.__name__: "double" for c in contracts}


@pytest.fixture(scope="module", params=KINDS)
def written(request, tmp_path_factory):
    out = tmp_path_factory.mktemp("figure")
    figure_stem = stem("de", f"{request.param}_demo", 1, day=date(2026, 9, 26))
    data = _contracts()[request.param]
    files = write_figure(data, out, figure_stem, title="Demo", data_releases={})
    return out, figure_stem, data, files


def test_on_page_is_points_plus_three_in_pixels():
    """A plot on the page draws text 3 pt larger than the figure, in px (96 px to 72 pt).

    The first version added 4 to a size in pt and used it as px: a 12 pt tick came out at
    16 px, which is 12 pt, not the 15 pt the rule asks for.
    """
    from ctx.viz.style import on_page  # noqa: PLC0415

    # 3 pt larger, in px, then 25 % larger (Rytis 2026-09-28 and 2026-09-30).
    assert on_page(9) == pytest.approx(20.0)
    assert on_page(10) == pytest.approx(65 / 3)
    assert "--size-base: 18px;" in css_variables(THEMES["white"])


def test_caption_wrap_keeps_words_whole():
    caption = demo_box(counts(), planted=1, null=1).caption
    for width in (30, 36, 40, 73):
        assert "BH-adjusted" in caption.text(width).replace("\n", " ")


def test_static_files_measured(written):
    out, figure_stem, data, _ = written
    caption = re.sub(r"\s+", " ", static_caption(data))
    for theme in THEMES.values():
        pdf = out / f"{figure_stem}_{theme.name}.pdf"
        assert _pdf_fonts(pdf) <= {f.replace(" ", "") for f in FONT_CHAIN} | {"ArialMT"}
        colours, sizes, widths, page_fill = _painted(pdf)
        allowed = [_rgb(c) for c in (*DATA_COLOURS, *CIVIDIS, theme.ink, theme.grid,
                                     theme.reference)]  # fmt: skip
        assert all(any(_same(c, a) for a in allowed) for c in colours), colours
        assert any(_same(c, _rgb(theme.ink)) for c in colours)
        assert sizes and sizes <= set(SIZES.values())
        assert widths and widths <= set(LINE_WIDTHS.values())
        assert not page_fill, "the canvas must stay transparent"
        page = PdfReader(pdf).pages[0]
        assert page["/Resources"].get("/XObject") is None or all(
            x.get_object()["/Subtype"] != "/Image" for x in page["/Resources"]["/XObject"].values()
        )
        width = WIDTHS_MM[column_width(data)] / 25.4 * 72
        assert float(page.mediabox.width) == pytest.approx(width, abs=0.5)
        text = re.sub(r"\s+", " ", page.extract_text())
        assert caption in text
        for label, unit in axis_titles(data):
            # The printed title is the one rule's; every word of the unit must be on the page,
            # in the title itself when the label does not already name it.
            assert axis_title(label, unit) in text, (label, unit)
            assert all(word in axis_title(label, unit).lower() for word in _unit_words(unit))

        svg = ElementTree.parse(out / f"{figure_stem}_{theme.name}.svg")
        families = {
            tuple(name.strip().strip("'\"") for name in match.split(","))
            for element in svg.iter()
            for match in re.findall(r"font-family:\s*([^;]+)", element.get("style", ""))
        }
        assert families == {(*FONT_CHAIN, "sans-serif")}, families
        assert not [e for e in svg.iter() if e.tag.endswith("image")]


OVERLAP_CASES = {
    "box_1_1": lambda s: demo_box(s, planted=1, null=1),
    "box_3_5": lambda s: demo_box(s, planted=3, null=5),
    "volcano_single": demo_volcano,
    "volcano_double": lambda s: demo_volcano(s, width="double"),
    "volcano_crowded": lambda s: demo_volcano(s, labels=10),
    "pca": demo_pca,
    "ma_log_x": demo_ma,
    "bar": demo_bar,
    "bar_horizontal": lambda s: demo_bar(s, horizontal=True),
    "composition_20_groups": lambda s: demo_bars_at_scale(horizontal=True),
    "library_sizes_210_runs": lambda s: demo_bars_at_scale(horizontal=False),
    "bar_filtered_lines": lambda s: filtered_bars(),
    "bar_filtered_lines_horizontal": lambda s: dataclasses.replace(
        filtered_bars("groups"), horizontal=True
    ),
    "strip": lambda s: demo_strip(),
    "strip_log": lambda s: demo_strip(log=True),
    "bubble": lambda s: demo_bubble(),
    "curves": lambda s: demo_curves(),
    "curves_double": lambda s: demo_curves(width="double"),
    "heatmap": lambda s: demo_heatmap(),
    "heatmap_unlabelled_rows": lambda s: demo_heatmap(n_genes=80),
    "band_map": lambda s: demo_band_map(),
    "flow": lambda s: demo_flow(),
    "surface": lambda s: demo_surface(),
    "surface_double": lambda s: demo_surface(width="double"),
    "density": lambda s: demo_density(),
    "linked": lambda s: demo_linked(),
}


@pytest.mark.parametrize("case", list(OVERLAP_CASES))
def test_no_text_overlap(case):
    data = OVERLAP_CASES[case](counts())
    for theme in THEMES.values():
        fig = data.static(theme)
        assert len(drawn_texts(fig)) > 4
        assert layout_problems(fig) == []


def test_write_refuses_clipped_text(tmp_path):
    data = demo_bar(counts())
    # Wrapped by measured width, 30 repeats took 36 lines, which tight_layout still fitted by
    # squeezing the axes to 9 % of the width; 60 take 71, wider than the page.
    long = dataclasses.replace(data, value_label="A value axis title far too long to fit " * 60)
    # matplotlib says so first: tight_layout gives up on a title it cannot fit, even at 170 mm.
    with (pytest.raises(ValueError, match="off the page"),
          pytest.warns(UserWarning, match="Tight layout not applied")):  # fmt: skip
        write_figure(long, tmp_path, stem("de", "clipped", 1, day=date(2026, 9, 27)),
                     title="t", data_releases={})  # fmt: skip
    assert not list(tmp_path.iterdir())


# --- scatter: drawn points, thinning, slider --------------------------------------------


def _points_in(fig) -> int:
    return sum(len(c.get_offsets()) for ax in fig.axes for c in ax.collections)


def _marker_traces(fig):
    return [t for t in fig.data if t.mode == "markers"]


def test_ordinals_read_as_written():
    """Mutation: the suffix taken from the last digit alone, so 11 reads "11st"."""
    from ctx.viz.scatter import _ordinal  # noqa: PLC0415

    cases = {1: "1st", 2: "2nd", 3: "3rd", 4: "4th", 11: "11th", 12: "12th", 13: "13th",
             21: "21st", 22: "22nd", 23: "23rd", 101: "101st", 111: "111th"}  # fmt: skip
    assert {n: _ordinal(n) for n in cases} == cases


def _trend(n: int = 12_796) -> ScatterData:
    """A dispersion plot of n genes with DESeq2's parametric trend through every one of them, as
    M5 draws Step-1's featured contrast (12,796 genes)."""
    rng = np.random.default_rng(7)
    mean = np.sort(10 ** rng.uniform(-0.5, 5.0, n))
    fitted = 0.04 + 1.8 / mean
    frame = pl.DataFrame({"mean": mean, "dispersion": fitted * rng.lognormal(0.0, 0.6, n),
                          "estimate": ["gene-wise"] * n, "keep": [False] * n})  # fmt: skip
    return ScatterData(
        frame=frame, columns={"mean": ("reads", "Mean."), "dispersion": ("none", "Dispersion."),
                              "estimate": ("none", "Estimate."), "keep": ("none", "Never.")},
        x="mean", y="dispersion", x_label="Mean expression", y_label="Dispersion",
        colour="estimate", colours=category_colours(["gene-wise"]), caption="Made up.",
        lines=(Line("fitted trend", tuple(mean), tuple(fitted)),), log_x=True, log_y=True,
        width="double", keep="keep",
    )  # fmt: skip


def _farthest(full: np.ndarray, kept: np.ndarray) -> float:
    """The farthest any point of ``full`` lies from the polyline through ``kept``, both (n, 2)."""
    a, b = kept[:-1], kept[1:]
    d = b - a
    rel = full[:, None, :] - a[None, :, :]
    t = np.clip((rel * d).sum(-1) / np.maximum((d * d).sum(-1), 1e-12), 0.0, 1.0)
    return float(np.sqrt(((rel - t[..., None] * d) ** 2).sum(-1)).min(axis=1).max())


def test_a_long_line_is_drawn_with_the_points_it_needs():
    """Rytis, 2026-10-01: M5's dispersion trend carried all 12,796 points onto the page. The page
    draws it through a few dozen, none of the 12,796 further than a quarter of a px from it at
    the largest size a page draws the plot (twice the widest page plot, opened alone); the
    figure files still draw every point, and a two-point line is left as it is. Mutations: the
    line drawn whole; the tolerance taken in axis units, not px; the figure files thinned too."""
    data = _trend()
    (trace,) = [t for t in data.interactive().data if t.mode == "lines"]
    line = data.lines[0]
    assert 2 < len(trace.x) < len(line.x) / 100
    assert (trace.x[0], trace.x[-1]) == (line.x[0], line.x[-1])
    # The line's own extent over the whole layout: more px per decade than the plot gives it,
    # whose axes span the points too, inside margins; so a bound on what a reader sees.
    size = np.array([max(PAGE_PLOT_PX.values()) * 2, PLOTLY_HEIGHT_PX["scatter"] * 2])
    full = np.log10(np.column_stack([line.x, line.y]))
    kept = np.log10(np.column_stack([trace.x, trace.y]))
    low, span = full.min(0), np.ptp(full, 0)
    assert _farthest((full - low) / span * size, (kept - low) / span * size) <= 0.25
    (drawn,) = [ln for ln in data.static(THEMES["white"]).axes[0].get_lines()
                if len(ln.get_xdata()) > 2]  # fmt: skip
    assert len(drawn.get_xdata()) == len(line.x)
    short = dataclasses.replace(data, lines=(Line("x = y", (1.0, 1e4), (1.0, 1e4)),))
    (two,) = [t for t in short.interactive().data if t.mode == "lines"]
    assert (two.x, two.y) == ((1.0, 1e4), (1.0, 1e4))


def test_points_drawn_equal_csv_rows():
    data = demo_volcano(counts(n_genes=6000))
    other = data.frame.filter(~pl.col("keep")).height
    assert other > THIN_MAX, "the case must thin"
    assert _points_in(data.static(THEMES["white"])) == data.table().height
    # With live calls the page draws its points from the rows it embeds (GS.3).
    assert sum(len(t.x) for t in _marker_traces(data.interactive())) == 0
    drawn = len(figure_block(data, "x").rows)
    stride = math.ceil(other / THIN_MAX)
    assert drawn == data.table().height - other + math.ceil(other / stride)
    # The caption read "every 3th" (review, 2026-10-01).
    assert stride == 3 and "every 3rd in order of" in data.caption_text()
    called = set(data.frame.filter(pl.col("keep"))["gene_id"])
    on_page = set(data.on_page()["gene_id"])
    assert called <= on_page
    # The page embeds the rows it draws, not the table's: the CSV beside it holds all.
    assert len(figure_block(data, "x").rows) == data.on_page().height < data.table().height


def test_empty_bars_refused():
    data = demo_bar(counts())
    with pytest.raises(ValueError, match="no bars"):
        dataclasses.replace(data, frame=data.frame.head(0))


@pytest.mark.parametrize(("kind", "horizontal"), [("runs", False), ("groups", True)])
def test_bar_lines_drawn_alike_in_both_renderers(kind, horizontal):
    """Reference lines, then the median line: one style each, named in both legends."""
    data = dataclasses.replace(filtered_bars(kind), horizontal=horizontal)
    totals = data.frame.group_by("bar").agg(pl.col("reads").sum())["reads"]
    names = ["thirty", "twenty", "median of the bars shown"]
    at = [30.0, 20.0, pytest.approx(float(np.median(totals)))]
    assert [(v, n) for v, n, _ in data.rules()] == list(zip(at, names, strict=True))
    ax = data.static(THEMES["white"]).axes[0]
    lines = ax.get_lines()
    assert [(line.get_xdata() if horizontal else line.get_ydata())[0] for line in lines] == at
    assert [line.get_linestyle() for line in lines] == ["--", ":", "-"]
    assert [t.get_text() for t in ax.get_legend().get_texts()][-3:] == names
    shapes = data.interactive().layout.shapes
    assert [s.name for s in shapes] == names and all(s.showlegend for s in shapes)
    assert [s.x0 if horizontal else s.y0 for s in shapes] == at
    assert [s.line.dash for s in shapes] == ["dash", "dot", "solid"]
    assert data.controls()["median_shape"] == len(shapes) - 1


def test_bar_hline_alone_as_frozen_callers_pass_it():
    # src/ctx/de/run.py and src/ctx/function/run.py pass hline= and nothing new.
    data = demo_bar(counts())
    assert data.rules() == [(3.0, "half a group", ("--", "dash"))]
    (shape,) = data.interactive().layout.shapes
    assert (shape.y0, shape.name, shape.line.dash) == (3.0, "half a group", "dash")
    assert data.controls() is None
    with pytest.raises(ValueError, match="at most 3 reference lines"):
        dataclasses.replace(data, lines=((1.0, "a"), (2.0, "b"), (4.0, "c")))
    with pytest.raises(ValueError, match="missing columns"):
        dataclasses.replace(data, filters=("donor",))


def test_bar_filter_controls():
    data = filtered_bars("groups")
    controls = data.controls()
    assert controls["filters"] == {"cell_type": ["LDN", "Monocyte", "NDN"],
                                   "condition": ["pdac", "steady"]}  # fmt: skip
    # Stacked: one bar per group, its length the total of both batches.
    assert controls["total"] == 6 and controls["traces"] == ["short", "long"]
    assert controls["labels"] == data.frame["bar"].to_list()
    assert controls["start"]["shown"] == controls["start"]["total"] == "6"
    assert controls["start"]["median"].endswith(" samples")


def _live_volcano(padj: list | None = None) -> ScatterData:
    frame = pl.DataFrame(
        {
            "gene_id": [f"g{i}" for i in range(10)],
            "lfc": [1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 1.5, 3.0, -3.0, 0.999],
            "p": [1e-4, 2e-4, 1e-3, 0.02, 0.01, 0.03, 0.2, 0.001, 0.5, 1e-5],
            "padj": padj or [0.001, 0.001, 0.01, 0.05, 0.04, 0.1, 0.2, None, float("nan"), 1e-4],
        }
    ).with_columns(y=-pl.col("p").log10(), call=classify(pl.col("padj"), pl.col("lfc")))
    # None and NaN are untested: the page embeds NaN as null.
    tested = frame.filter(pl.col("padj").is_not_nan())
    passing = tested.filter(pl.col("padj") < 0.05).height
    calls = Calls(
        lfc="lfc", p="p", adjusted={"padj": "BH-adjusted p"}, raw_words="raw p", tested="case",
        reference="control", alpha=0.05, lfc_min=1.0, max_cutoff=0.1, n_tested=tested.height,
        p_floor=1e-300, p_line=p_line("padj", 0.05, tested.height, passing),
        reading=("{n_up} up in {tested}; {line_rule}",),
        summary="{n_up} up, {n_down} down", line_rule="line at {y_line}",
        no_line_rule="no line", unit="genes",
    )  # fmt: skip
    return ScatterData(
        frame=frame, columns={c: ("none", c) for c in frame.columns}, x="lfc", y="y",
        x_label="LFC", y_label="y", colour="call", colours={"up": "#D55E00", "down": "#0072B2",
        "ns": None}, caption="c", hover=("gene_id",), calls=calls,
    )  # fmt: skip


def test_live_volcano_default_from_figure():
    """GS.3: the point traces carry no points, which the page script draws from the rows; the
    default colours are the call column's. Mutation: the traces' points left in."""
    data = _live_volcano()
    fig = data.interactive()
    controls = data.controls()
    assert [len(fig.data[i].x) for i in controls["traces"].values()] == [0, 0, 0]
    calls = data.frame["call"].to_list()
    # Written here, independently of classify: padj strictly below 0.05, |LFC| at least 1.
    rows = [
        (p, fc, call)
        for p, fc, call in data.frame.select("padj", "lfc", "call").rows()
        if p is not None and not math.isnan(p)
    ]
    assert calls.count("up") == sum(p < 0.05 and fc >= 1 for p, fc, _ in rows) == 2
    assert calls.count("down") == sum(p < 0.05 and fc <= -1 for p, fc, _ in rows) == 1
    fc_up, fc_down, p_line = (fig.data[i] for i in controls["cutoff_lines"])
    assert (fc_up.x, fc_down.x) == ((1.0, 1.0), (-1.0, -1.0))
    # The bar of the 5th of 8 tested genes: 5 x 0.05 / 8, written without p_line.
    assert p_line.y == pytest.approx((-math.log10(5 * 0.05 / 8),) * 2)
    assert data.reading() == ("2 up in case; line at 1.51",)
    assert data.summary() == "2 up, 1 down"
    assert list(data.page_table().columns) == ["gene_id", "lfc", "p", "padj"]


def test_live_volcano_refuses_what_the_page_cannot_redraw():
    data = _live_volcano()
    frame = data.frame.with_columns(shape=pl.lit("a"))
    with pytest.raises(ValueError, match="symbol"):
        dataclasses.replace(data, frame=frame, symbol="shape")
    with pytest.raises(ValueError, match="reading"):
        figure_block(data, "x", FigureText(reading=("A line the page would overwrite.",)))


def test_live_volcano_line_when_nothing_passes():
    # D-015: the first bar, 0.05 / 8, which no gene clears; not a hidden line.
    data = _live_volcano([0.5] * 7 + [None, float("nan"), 0.5])
    line = data.interactive().data[data.controls()["cutoff_lines"][2]]
    assert line.visible is not False
    assert line.y == pytest.approx((-math.log10(0.05 / 8),) * 2)
    assert data.reading() == ("0 up in case; line at 2.20",)


def test_live_volcano_hides_a_missing_line():
    # The contract's own no-line state; p_line gives it only when no gene was tested.
    data = _live_volcano()
    data = dataclasses.replace(data, calls=dataclasses.replace(data.calls, p_line=None))
    assert data.interactive().data[data.controls()["cutoff_lines"][2]].visible is False
    assert data.reading() == ("2 up in case; no line",)


def test_classify():
    frame = pl.DataFrame(
        {"padj": [0.01, 0.01, 0.05, None, float("nan"), 0.001], "lfc": [1.0, -1.0, 2, 3, 3, 0]}
    )
    called = frame.select(classify(pl.col("padj"), pl.col("lfc"))).to_series().to_list()
    assert called == ["up", "down", "ns", "ns", "ns", "ns"]


# --- PERMANOVA ------------------------------------------------------------------------------


def _clusters(seed: int = 0):
    rng = np.random.default_rng(seed)
    coords = rng.normal(size=(20, 3))
    coords[:10] += 3.0
    return coords, ["a"] * 10 + ["b"] * 10


def test_permanova_identity():
    coords, labels = _clusters()
    test = permanova(coords, labels)
    # Anderson (2001): sums of squares from squared distances, written without group means.
    d2 = ((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2)
    n, groups = len(labels), sorted(set(labels))
    ss_total = d2[np.triu_indices(n, 1)].sum() / n
    ss_within = sum(
        d2[np.ix_(idx, idx)][np.triu_indices(len(idx), 1)].sum() / len(idx)
        for idx in ([i for i, g in enumerate(labels) if g == level] for level in groups)
    )
    f = ((ss_total - ss_within) / (len(groups) - 1)) / (ss_within / (n - len(groups)))
    assert test.f == pytest.approx(f)
    assert test.r2 == pytest.approx((ss_total - ss_within) / ss_total)


def test_permanova_null_and_planted():
    coords, labels = _clusters()
    assert permanova(coords, labels).p == pytest.approx(1 / 1000)
    shuffled = list(np.random.default_rng(42).permutation(labels))
    null = np.random.default_rng(7).normal(size=(20, 3))
    assert permanova(null, shuffled).p > 0.05


def test_permanova_strata():
    coords, labels = _clusters()
    # Labels constant within each stratum: no permutation can move them, so there is no
    # p-value, and the caption says so rather than printing p = 1 as if it were a result.
    frozen = permanova(coords, labels, strata=labels)
    assert math.isnan(frozen.p) and frozen.r2 > 0
    assert "no p-value" in frozen.text("group") and "p = 1" not in frozen.text("group")
    paired = [f"s{i % 10}" for i in range(20)]
    assert permanova(coords, labels, strata=paired).p < 0.05


def test_fit_title_wraps_to_axis_length():
    from ctx.viz.annotate import fit_title  # noqa: PLC0415

    height = WIDTHS_MM["single"] * 1.1 / 25.4
    # Longer than the axis holds at the axis-title size, whatever that size is: an x is half
    # an em wide in Arial and Liberation Sans.
    fits = math.ceil(height * 72 * AXIS_TITLE_SHARE / (SIZES["axis_label"] / 2))
    title = "mono_hsct_early " + "x" * fits + " (log2 FC)"
    wrapped = fit_title(title, height)
    assert wrapped.replace("\n", " ") == title
    assert "\n" in wrapped
    assert fit_title("Samples (runs)", height) == "Samples (runs)"


def test_bubble_area_grows_linearly_with_the_size_value():
    from ctx.viz.style import BUBBLE_PT  # noqa: PLC0415

    low, mid, high = bubble_pt([0.0, 5.0, 10.0], 10.0)
    assert (low, high) == pytest.approx((BUBBLE_PT["min"], BUBBLE_PT["max"]))
    assert mid**2 - low**2 == pytest.approx((high**2 - low**2) / 2)
    data = demo_bubble()
    with pytest.raises(ValueError, match="non-negative"):
        dataclasses.replace(data, frame=data.frame.with_columns(mt_rrna_share=pl.lit(-0.1)))
    diameters = [t.marker.size for t in _marker_traces(data.interactive())]
    assert all(isinstance(d, tuple | list) and len(d) > 0 for d in diameters)


def test_a_curve_keeps_one_colour():
    data = demo_curves()
    mixed = data.frame.with_columns(
        group=pl.when(pl.col("position") == 1).then(pl.lit("degraded")).otherwise("group")
    )
    with pytest.raises(ValueError, match="more than one group"):
        dataclasses.replace(data, frame=mixed)
    traces = [t for t in data.interactive().data if t.mode == "lines" and t.legendgroup]
    assert len(traces) == data.frame["sample"].n_unique()
    assert sum(bool(t.showlegend) for t in traces) == 2


def test_cluster_order_keeps_planted_blocks_together():
    from ctx.viz.heatmap import cluster_order  # noqa: PLC0415

    rng = np.random.default_rng(0)
    pattern = rng.normal(size=20)
    rows = [pattern + rng.normal(0, 0.2, 20) if i % 2 else -pattern + rng.normal(0, 0.2, 20)
            for i in range(10)]  # fmt: skip
    labels = [f"{'a' if i % 2 else 'b'}{i}" for i in range(10)]
    order = cluster_order(np.vstack([*rows, np.ones(20)]), [*labels, "flat"])
    blocks = [name[0] for name in order[:-1]]
    assert blocks in (["a"] * 5 + ["b"] * 5, ["b"] * 5 + ["a"] * 5)
    assert order[-1] == "flat"


def test_diverging_levels_bin_at_their_edges():
    from ctx.viz.heatmap import diverging_levels  # noqa: PLC0415

    values = [-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0, None, float("nan")]
    levels, swatches = diverging_levels(values)
    assert levels == ["< -2", "-2 to -1", "-1 to -0.5", "-0.5 to 0.5", "-0.5 to 0.5",
                      "0.5 to 1", "1 to 2", "≥ 2", "≥ 2", None, None]  # fmt: skip
    assert [s.colour for s in swatches] == [ROLES["reference"]] * 3 + [None] + [
        ROLES["comparison"]] * 3  # fmt: skip
    assert [s.alpha for s in swatches[:3]] == sorted(s.alpha for s in swatches[4:])[::-1]


def test_a_flow_keeps_every_read():
    data = demo_flow()
    lost = data.links.with_columns(
        value=pl.when(pl.col("target") == "decoy").then(0.0).otherwise("value")
    )
    with pytest.raises(ValueError, match=r"differ at \['kept'\]"):
        dataclasses.replace(data, links=lost)
    backwards = pl.concat([data.links, pl.DataFrame(
        {"source": ["assigned"], "target": ["kept"], "value": [1.0]})])  # fmt: skip
    with pytest.raises(ValueError, match="later stage"):
        dataclasses.replace(data, links=backwards)
    sankey = data.interactive().data[0]
    assert sum(sankey.link.value) == data.links["value"].sum()


def test_cividis_stops_are_the_colormap_at_band_midpoints():
    stops = [mpl.colormaps["cividis"]((i + 0.5) / len(CIVIDIS))[:3] for i in range(len(CIVIDIS))]
    assert all(_same(_rgb(c), tuple(s)) for c, s in zip(CIVIDIS, stops, strict=True))


def test_hdr_levels_hold_their_mass_and_no_more_cells():
    from ctx.viz.contour import density_grid, hdr_levels  # noqa: PLC0415

    rng = np.random.default_rng(1)
    x, y = rng.normal(size=5000), rng.normal(size=5000)
    _, _, density = density_grid(x, y, (-5.0, 5.0, -5.0, 5.0))
    for mass, level in zip((0.5, 0.9), hdr_levels(density), strict=True):
        assert density[density >= level].sum() >= mass
        assert density[density > level].sum() < mass
    # Smoothing spreads the mass, so the points inside run a little above it at this n.
    data = demo_density()
    inner = (data.table()["region"] == "inner").sum() / data.table()["region"].is_not_null().sum()
    assert 0.5 <= inner <= 0.75


@pytest.mark.parametrize("title", ["Difference d (log2(count + 1))",
                                   "Count with -O, split (--fraction) (log2(count + 1))",
                                   "Called strand (log2(count + 1))"])  # fmt: skip
def test_a_vertical_title_never_breaks_inside_its_unit(title):
    """At every scatter height a wrapped y title keeps "(log2(count + 1))" on one line; split at
    its spaces it had left "1))" or "+ 1))" alone (M3's agreement plots). At 0.3 of the height
    the unit, 73.6 pt, is wider than a line, 43.7 pt; at 0.6 and 1.0 it fits one, where a
    measured wrap at every space also kept it whole.

    Mutation caught: ``fit_title`` wrapping at every space again."""
    from ctx.viz.annotate import fit_title  # noqa: PLC0415
    from ctx.viz.style import figure_inches  # noqa: PLC0415

    for width in ("single", "double"):
        _, height = figure_inches(width, "scatter")
        for scale in (0.3, 0.6, 1.0):
            lines = fit_title(title, height * scale).split("\n")
            assert "(log2(count + 1))" in "\n".join(lines), lines
            assert " ".join(lines) == title


def test_a_long_interactive_title_puts_its_unit_on_a_line_of_its_own():
    """Plotly never wraps a title: past 26 characters the page breaks before the unit, so a
    375 px phone shows it whole (M3.7); a short title and one without a unit stay as they are.

    Mutation caught: ``screen_title`` returning the title unchanged (M3's x titles ran 5 to 24
    px past the figure at 375 px)."""
    from ctx.viz.axes import screen_title  # noqa: PLC0415

    assert screen_title("Mean of both (log2(count + 1))") == "Mean of both<br>(log2(count + 1))"
    assert screen_title("Reads (%)") == "Reads (%)"
    assert screen_title("A label of thirty characters") == "A label of thirty characters"
