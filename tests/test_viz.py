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
    demo_bar,
    demo_bars_at_scale,
    demo_box,
    demo_pca,
    demo_volcano,
)

from ctx.de.stats import classify
from ctx.io.naming import stem
from ctx.report.builder import FigureText
from ctx.viz import box
from ctx.viz.annotate import compare_groups, permanova
from ctx.viz.axes import axis_title
from ctx.viz.save import drawn_texts, figure_block, layout_problems, write_figure
from ctx.viz.scatter import THIN_MAX, Calls, ScatterData
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


def _contracts():
    synthetic = counts()
    return {
        "box": demo_box(synthetic, planted=3, null=5),
        "volcano_single": demo_volcano(synthetic),
        "volcano_double": demo_volcano(synthetic, width="double"),
        "pca": demo_pca(synthetic),
        "bar": demo_bar(synthetic),
        "bar_horizontal": demo_bar(synthetic, horizontal=True),
        "composition_20_groups": demo_bars_at_scale(horizontal=True),
        "library_sizes_210_runs": demo_bars_at_scale(horizontal=False),
    }


KINDS = tuple(_contracts())


def static_caption(data) -> str:
    """The caption the static figure prints."""
    return data.caption.text() if isinstance(data, box.BoxData) else data.caption


def axis_titles(data) -> list[tuple[str, str]]:
    """(label, unit) of every value axis the static figure must print."""
    if isinstance(data, box.BoxData):
        return [(data.value_label, data.unit)]
    if isinstance(data, ScatterData):
        return [(data.x_label, data.columns[data.x][0]), (data.y_label, data.columns[data.y][0])]
    return [(data.value_label, data.columns[data.value][0])]


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


def column_width(data) -> str:
    """The column width the figure should take."""
    if isinstance(data, box.BoxData):
        many = len(data.caption.tests) > SINGLE_COLUMN_MAX_FEATURES
        return "double" if many else "single"
    return data.width


@pytest.fixture(scope="module", params=KINDS)
def written(request, tmp_path_factory):
    out = tmp_path_factory.mktemp("figure")
    figure_stem = stem("de", f"{request.param}_demo", 1, day=date(2026, 9, 26))
    data = _contracts()[request.param]
    files = write_figure(data, out, figure_stem, title="Demo", data_releases={})
    return out, figure_stem, data, files


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
            name.strip().strip("'\"")
            for element in svg.iter()
            for match in re.findall(r"font-family:\s*([^;]+)", element.get("style", ""))
            for name in match.split(",")
        }
        assert families and families <= {*FONT_CHAIN, "sans-serif"}
        assert not [e for e in svg.iter() if e.tag.endswith("image")]


OVERLAP_CASES = {
    "box_1_1": lambda s: demo_box(s, planted=1, null=1),
    "box_3_5": lambda s: demo_box(s, planted=3, null=5),
    "volcano_single": demo_volcano,
    "volcano_double": lambda s: demo_volcano(s, width="double"),
    "volcano_crowded": lambda s: demo_volcano(s, labels=10),
    "pca": demo_pca,
    "bar": demo_bar,
    "bar_horizontal": lambda s: demo_bar(s, horizontal=True),
    "composition_20_groups": lambda s: demo_bars_at_scale(horizontal=True),
    "library_sizes_210_runs": lambda s: demo_bars_at_scale(horizontal=False),
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
    long = dataclasses.replace(data, value_label="A value axis title far too long to fit " * 3)
    with pytest.raises(ValueError, match="off the page"):
        write_figure(long, tmp_path, stem("de", "clipped", 1, day=date(2026, 9, 27)),
                     title="t", data_releases={})  # fmt: skip
    assert not list(tmp_path.iterdir())


# --- scatter: drawn points, thinning, slider --------------------------------------------


def _points_in(fig) -> int:
    return sum(len(c.get_offsets()) for ax in fig.axes for c in ax.collections)


def _marker_traces(fig):
    return [t for t in fig.data if t.mode == "markers"]


def test_points_drawn_equal_csv_rows():
    data = demo_volcano(counts(n_genes=6000))
    other = data.frame.filter(~pl.col("keep")).height
    assert other > THIN_MAX, "the case must thin"
    assert _points_in(data.static(THEMES["white"])) == data.table().height
    drawn = sum(len(t.x) for t in _marker_traces(data.interactive()))
    stride = math.ceil(other / THIN_MAX)
    assert drawn == data.table().height - other + math.ceil(other / stride)
    assert f"every {stride}th" in data.caption_text()
    called = set(data.frame.filter(pl.col("keep"))["gene_id"])
    on_page = set(data.on_page()["gene_id"])
    assert called <= on_page
    # The page embeds the rows it draws, not the table's: the CSV beside it holds all.
    assert len(figure_block(data, "x").rows) == data.on_page().height < data.table().height


def test_empty_bars_refused():
    data = demo_bar(counts())
    with pytest.raises(ValueError, match="no bars"):
        dataclasses.replace(data, frame=data.frame.head(0))


def _live_volcano() -> ScatterData:
    frame = pl.DataFrame(
        {
            "gene_id": [f"g{i}" for i in range(10)],
            "lfc": [1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 1.5, 3.0, -3.0, 0.999],
            "p": [1e-4, 2e-4, 1e-3, 0.02, 0.01, 0.03, 0.2, 0.001, 0.5, 1e-5],
            "padj": [0.001, 0.001, 0.01, 0.05, 0.04, 0.1, 0.2, None, float("nan"), 0.0001],
        }
    ).with_columns(y=-pl.col("p").log10(), call=classify(pl.col("padj"), pl.col("lfc")))
    calls = Calls(
        lfc="lfc", p="p", adjusted={"padj": "BH-adjusted p"}, raw_words="raw p", tested="case",
        reference="control", alpha=0.05, lfc_min=1.0, max_cutoff=0.1, n_tested=7,
        p_floor=1e-300, p_line=0.01, reading=("{n_up} up in {tested}; {line_rule}",),
        summary="{n_up} up, {n_down} down", line_rule="line at {y_line}",
        no_line_rule="no line",
    )  # fmt: skip
    return ScatterData(
        frame=frame, columns={c: ("none", c) for c in frame.columns}, x="lfc", y="y",
        x_label="LFC", y_label="y", colour="call", colours={"up": "#D55E00", "down": "#0072B2",
        "ns": None}, caption="c", hover=("gene_id",), calls=calls,
    )  # fmt: skip


def test_live_volcano_default_from_figure():
    data = _live_volcano()
    fig = data.interactive()
    controls = data.controls()
    xs = {fig.data[i].name: fig.data[i].x for i in controls["traces"].values()}
    # Written here, independently of classify: padj strictly below 0.05, |LFC| at least 1.
    rows = [
        (p, fc)
        for p, fc in data.frame.select("padj", "lfc").rows()
        if p is not None and not math.isnan(p)
    ]
    assert len(xs["up"]) == sum(p < 0.05 and fc >= 1 for p, fc in rows) == 2
    assert len(xs["down"]) == sum(p < 0.05 and fc <= -1 for p, fc in rows) == 1
    fc_up, fc_down, p_line = (fig.data[i] for i in controls["cutoff_lines"])
    assert (fc_up.x, fc_down.x) == ((1.0, 1.0), (-1.0, -1.0))
    assert p_line.y == pytest.approx((2.0, 2.0))  # -log10 of p_line = 0.01
    assert data.reading() == ("2 up in case; line at 2.00",)
    assert data.summary() == "2 up, 1 down"
    assert list(data.page_table().columns) == ["gene_id", "lfc", "p", "padj"]


def test_live_volcano_refuses_what_the_page_cannot_redraw():
    data = _live_volcano()
    frame = data.frame.with_columns(shape=pl.lit("a"))
    with pytest.raises(ValueError, match="symbol"):
        dataclasses.replace(data, frame=frame, symbol="shape")
    with pytest.raises(ValueError, match="reading"):
        figure_block(data, "x", FigureText(reading=("A line the page would overwrite.",)))


def test_live_volcano_no_line_when_nothing_passes():
    data = _live_volcano()
    calls = dataclasses.replace(data.calls, p_line=None)
    data = dataclasses.replace(data, calls=calls)
    fig = data.interactive()
    p_line = fig.data[data.controls()["cutoff_lines"][2]]
    assert p_line.visible is False
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

    title = "mono_hsct_early (log2 FC)"
    height = WIDTHS_MM["single"] * 1.1 / 25.4
    wrapped = fit_title(title, height)
    assert wrapped.replace("\n", " ") == title
    assert "\n" in wrapped
    assert fit_title("Samples (runs)", height) == "Samples (runs)"
