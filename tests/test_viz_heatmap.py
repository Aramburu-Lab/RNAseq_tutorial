"""The page heatmap draws the static's layout (Rytis 2026-10-01, M7's paired z-scores): its
track floated far above the cells, every other row went unlabelled, the right-hand statistic
ran off the figure and the track swatches were missing from the legend. The review of that fix
(2026-10-01) added: the "nothing to show" grey stayed the white theme's on a dark page, every
row was as tall as the longest label, hover boxes ran off a phone, swatches were 6 px. The fix
wrapped every label at a phone's width on every screen; since the qc-page-fixes merge
(2026-10-07, work order item 2) a desktop keeps the static's labels and a phone swaps in the
wrapped ones (``layout.meta.narrow``), so the phone assertions read the narrow layout."""

import dataclasses
import itertools
import re
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from matplotlib.colors import to_rgba
from synth import demo_band_map, demo_heatmap

from ctx.viz.heatmap import GAP, HeatmapData, Swatch, Track
from ctx.viz.metrics import text_width, wrap_to
from ctx.viz.style import (
    BAND_COLOURS,
    BAR_LABEL_ROTATION,
    CHROME,
    HEAT,
    LEGEND,
    LINE_SPACING,
    MM_PER_INCH,
    NARROW_PX,
    PHONE_PLOT_PX,
    PHONE_PX,
    PLOTLY_HEAT_BASE_PX,
    PLOTLY_HEAT_ROW_PX,
    PLOTLY_HOVER_PX,
    PLOTLY_LEGEND_LINE_PX,
    PLOTLY_LEGEND_PAD_PX,
    PLOTLY_LEGEND_SWATCH_PX,
    PLOTLY_LEGEND_SYMBOL_PX,
    PLOTLY_MARGIN_PX,
    SIZES,
    THEMES,
    apply_matplotlib,
    css_layout,
    on_page,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = {
    "tracks_and_stats": demo_heatmap,
    "past_label_rows": lambda: demo_heatmap(n_genes=80),
    "band_map": demo_band_map,
}
# Metric names of the QC band maps that wrap on a phone.
LONG = ["Reads written by Trim Galore", "Per-base sequence content"]


def _band_map(metrics: list[str], *, drop_one: bool = False) -> HeatmapData:
    """A QC band map with these metrics as rows, labelled as the QC pages label them."""
    runs = [f"R{i:02d}" for i in range(4)]
    frame = pl.DataFrame([(r, m, "good") for r in runs for m in metrics], orient="row",
                         schema=["run", "metric", "band"]).with_columns(
        value=pl.lit(None, pl.Float64))  # fmt: skip
    if drop_one:  # a cell with no row in the table draws nothing, as in the static
        frame = frame.slice(1)
    return HeatmapData(
        frame=frame, row="metric", column="run", value="value", level="band",
        columns={"run": ("none", "Run."), "metric": ("none", "Metric."),
                 "band": ("none", "Band."), "value": ("none", "Not used.")},
        row_label="Registry metric.", column_label="Sample, as nf-core/rnaseq names it.",
        swatches=(Swatch("good", BAND_COLOURS["good"]), Swatch("no band", None)),
        caption="Made up.", row_order=tuple(metrics), column_order=tuple(runs),
    )  # fmt: skip


def _static(data: HeatmapData, theme: str = "white") -> dict[str, list]:
    """What the static figure draws: labels top down, the right column, legend, columns."""
    fig = data.static(THEMES[theme])
    ax = fig.axes[0]
    legends = [*fig.legends, *([ax.get_legend()] if ax.get_legend() else [])]
    return {
        "left": [t.get_text() for t in ax.get_yticklabels()],
        "right": [t.get_text() for a in ax.child_axes for t in a.get_yticklabels()],
        "right_title": [a.get_ylabel() for a in ax.child_axes],
        "legend": [t.get_text() for lg in legends for t in lg.get_texts()],
        "patches": [tuple(h.get_facecolor()) for lg in legends for h in lg.legend_handles],
        "columns": [t.get_text() for t in ax.get_xticklabels()],
    }


def _plain(text: str) -> str:
    """A plotly text with its line breaks read as the spaces they replaced."""
    return text.replace("<br>", " ")


def _lines(text: str) -> int:
    return text.count("<br>") + 1


def _legend(fig) -> list:
    return [t for t in fig.data if t.type == "scatter" and t.showlegend is not False]


def _rgba(css: str) -> tuple[float, ...]:
    """A CSS colour, hex or ``rgba(r,g,b,a)``, as matplotlib's 0-1 RGBA."""
    if not css.startswith("rgba("):
        return to_rgba(css)
    r, g, b, a = (float(v) for v in css.removeprefix("rgba(").removesuffix(")").split(","))
    return r / 255, g / 255, b / 255, a


def _over(rgba: tuple[float, ...], ground: str) -> np.ndarray:
    """``rgba`` laid over an opaque hex ground, as 0-255 RGB: what the reader sees."""
    base = np.array([int(ground[i : i + 2], 16) for i in (1, 3, 5)], float)
    return np.array(rgba[:3]) * 255 * rgba[3] + base * (1 - rgba[3])


def _steps(colorscale) -> list[str]:
    """The colour each integer level takes in a ``plotly_steps`` colourscale."""
    return [colour for _, colour in colorscale[::2]]


def _fits(text: str, width: float, size: float) -> bool:
    """Every line of a plotly text fits ``width`` px at ``size`` pt, as ctx.viz.metrics
    measures it, or is one word."""
    return all(text_width(line, on_page(size)) <= width or " " not in line
               for line in text.split("<br>"))  # fmt: skip


def _narrowed(fig):
    """The figure as the page script draws it narrower than ``meta.narrow.below``: its trace
    updates, then its layout's, as plotly's restyle and relayout apply them."""
    import plotly.graph_objects as go  # noqa: PLC0415

    narrow = fig.layout.meta["narrow"]
    out = go.Figure(fig)
    for update, traces in narrow["restyle"]:
        out.plotly_restyle(update, traces)
    out.plotly_relayout(narrow["relayout"])
    return out


def _same_colour(page: tuple[float, ...], white: tuple, black: tuple) -> None:
    """A page colour against the static's in both themes; one that differs is the grid grey."""
    if white == black:
        assert page == pytest.approx(white, abs=1e-3)
        return
    # The "nothing to show" level: one translucent grey, read over either theme's card.
    for theme, static in (("white", white), ("black", black)):
        seen = _over(page, CHROME[theme]["surface"])
        assert seen == pytest.approx(np.array(static[:3]) * 255, abs=4), theme


def test_rows_sit_where_the_static_draws_them():
    # Order of operations: the tracks, then the blank row, then the cells, as the static draws
    # them, each band starting where the one above ends. A categorical axis had stretched to
    # twice them, which was the band between the track and the cells.
    data = demo_heatmap()
    fig = data.interactive()
    left = _static(data)["left"]
    above = len(data.tracks) + GAP
    assert left[:above] == [t.name for t in data.tracks] + [""] * GAP
    cells, *tracks = [t for t in fig.data if t.type == "heatmap"]
    assert [list(t.y) for t in tracks] == [[i, i + 1] for i in range(len(data.tracks))]
    assert list(cells.y) == list(range(above, len(left) + 1))
    assert list(fig.layout.yaxis.range) == [len(left), 0]
    assert list(fig.layout.yaxis.tickvals) == [i + 0.5 for i, label in enumerate(left) if label]


@pytest.mark.parametrize("case", list(FIXTURES))
def test_rows_labelled_as_the_static_labels(case):
    data = FIXTURES[case]()
    fig = data.interactive()
    assert fig.layout.yaxis.tickmode == "array"  # plotly's own ticks dropped every other row
    labelled = [label for label in _static(data)["left"] if label]
    assert [_plain(t) for t in fig.layout.yaxis.ticktext] == labelled
    assert len(labelled) == len(data.tracks) + (
        len(data.rows_in_order()) if len(data.rows_in_order()) <= HEAT["label_rows"] else 0
    )


def test_a_long_label_raises_its_own_row_alone():
    # Every row had been as tall as the longest label: 3 pitches each on the QC band maps,
    # 1078 px for 9 rows at 1280 px (2026-10-01). On a phone's layout the long label wraps
    # and raises its row alone; on a desktop no label wraps and every row is one pitch.
    data = _band_map(["Q30 bases", "Per-base sequence content", "N bases"])
    fig = data.interactive()
    assert [_lines(t) for t in fig.layout.yaxis.ticktext] == [1, 1, 1]
    assert np.diff(fig.data[0].y).tolist() == [1, 1, 1]
    phone = _narrowed(fig)
    lines = [_lines(t) for t in phone.layout.yaxis.ticktext]
    # Measured, the long label takes 2 lines of half a phone's plot; at the average glyph, 3.
    assert lines == [1, 2, 1]
    assert np.diff(phone.data[0].y).tolist() == lines
    short = _narrowed(_band_map(["Q30 bases", "GC", "N bases"]).interactive())
    assert phone.layout.height - short.layout.height == pytest.approx(PLOTLY_HEAT_ROW_PX)


@pytest.mark.parametrize("case", list(FIXTURES))
def test_columns_ordered_as_the_static(case):
    data = FIXTURES[case]()
    xaxis = data.interactive().layout.xaxis
    # As categories: a column named like a number would otherwise sit at its value.
    assert (xaxis.type, xaxis.categoryorder) == ("category", "array")
    assert list(xaxis.categoryarray) == data.columns_in_order()
    shown = _static(data)["columns"]
    assert xaxis.showticklabels == bool(shown)
    assert not shown or shown == data.columns_in_order()
    assert xaxis.tickangle == -BAR_LABEL_ROTATION
    # plotly thins column labels to the width it draws at. Every column forced on overlapped
    # on a phone and on M4's 33 columns at 1280 px (93 and 32 pairs, measured 2026-10-01).
    assert xaxis.tickmode in (None, "auto")


def test_statistic_is_an_axis_inside_the_figure():
    data = demo_heatmap()
    fig = data.interactive()
    static = _static(data)
    cells = fig.data[0]
    y2 = fig.layout.yaxis2
    assert (y2.side, y2.overlaying, y2.matches, y2.automargin) == ("right", "y", "y", True)
    assert list(y2.tickvals) == [(a + b) / 2 for a, b in itertools.pairwise(cells.y)]
    assert list(y2.ticktext) == static["right"]
    assert y2.title.text == static["right_title"][0] == data.row_stats_label
    # plotly drops an axis no trace is on; annotations take no margin and ran off the figure.
    assert [t.yaxis for t in fig.data if t.yaxis == "y2"] == ["y2"]
    assert not fig.layout.annotations
    # The template's margins grow to fit the axis; a margin set here could stop them.
    assert "margin" not in fig.to_plotly_json()["layout"]


def test_no_statistic_no_right_axis():
    fig = demo_band_map().interactive()
    assert "yaxis2" not in fig.layout.to_plotly_json()
    assert not [t for t in fig.data if t.yaxis == "y2"]


@pytest.mark.parametrize("case", ["tracks_and_stats", "band_map"])
def test_legend_is_the_static_legend(case):
    data = FIXTURES[case]()
    entries = _legend(data.interactive())
    white, black = _static(data, "white"), _static(data, "black")
    assert [_plain(t.name) for t in entries] == white["legend"]
    for trace, on_white, on_black in zip(entries, white["patches"], black["patches"], strict=True):
        _same_colour(_rgba(trace.marker.color), on_white, on_black)


@pytest.mark.parametrize("case", [*FIXTURES, "missing_cell"])
def test_cells_draw_the_static_grid(case):
    data = _band_map(["Q30 bases", "N bases"], drop_one=True) if case == "missing_cell" else (
        FIXTURES[case]())  # fmt: skip
    grid, labels = data._grid()
    fig = data.interactive()
    cells, *tracks = [t for t in fig.data if t.type == "heatmap"]
    above = len(labels) - len(data.rows_in_order())
    assert list(cells.x) == data.columns_in_order()
    drawn = np.array([[np.nan if v is None else v for v in row] for row in cells.z], float)
    np.testing.assert_array_equal(drawn, grid[above:])
    white, black = data._rgba(THEMES["white"]), data._rgba(THEMES["black"])
    for k, colour in enumerate(_steps(cells.colorscale)):
        _same_colour(_rgba(colour), white[k], black[k])
    offset = len(data.swatches)
    for k, (track, trace) in enumerate(zip(data.tracks, tracks, strict=True)):
        np.testing.assert_array_equal(np.array(trace.z[0]) + offset, grid[k])
        assert _steps(trace.colorscale) == list(track.colours.values())
        offset += len(track.colours)


def test_nothing_to_show_reads_as_each_themes_grid():
    # The white theme's grid, #E6E6E6, had been the brightest cell on a dark page: 19 % of
    # M7's drawn cells and 92 % of M11's (2026-10-01).
    data = demo_band_map()
    level = [s.colour for s in data.swatches].index(None)
    shown = [_rgba(_steps(data.interactive().data[0].colorscale)[level])]
    shown.append(_rgba(_legend(data.interactive())[level].marker.color))
    for page in shown:
        for theme, other in (("white", "black"), ("black", "white")):
            seen = _over(page, CHROME[theme]["surface"])
            grid = [int(THEMES[theme].grid[i : i + 2], 16) for i in (1, 3, 5)]
            wrong = [int(THEMES[other].grid[i : i + 2], 16) for i in (1, 3, 5)]
            assert np.abs(seen - grid).max() <= 4, theme
            assert np.abs(seen - wrong).max() > 100, theme


def test_hover_names_the_cell_and_its_value():
    data = demo_heatmap()
    cells = data.interactive().data[0]
    rows, cols = data.rows_in_order(), data.columns_in_order()
    values = dict(zip(zip(data.frame["gene"], data.frame["sample"], strict=True),
                      data.frame["z"], strict=True))  # fmt: skip
    for i, j in itertools.product((0, len(rows) - 1), (0, len(cols) - 1)):
        name = f"Gene: {rows[i]} Sample: {cols[j]} z (SD): {values[rows[i], cols[j]]:.3g}"
        assert _plain(cells.text[i][j]).startswith(name)
    band = _band_map(["Q30 bases"]).interactive().data[0]
    # No value line where the map has none; a label written as a sentence loses its stop.
    assert _plain(band.text[0][0]) == (
        "Registry metric: Q30 bases Sample, as nf-core/rnaseq names it: R00 good"
    )


@pytest.mark.parametrize("case", ["tracks_and_stats", "long_labels"])
def test_hover_lines_fit_a_phone(case):
    # At 375 px a band map's hover ran from -141.9 to 401.3 px of a 317 px plot (2026-10-01).
    data = demo_heatmap() if case == "tracks_and_stats" else _band_map(LONG)
    texts = [t for trace in data.interactive().data if trace.type == "heatmap"
             for row in trace.text for t in row]  # fmt: skip
    assert texts
    assert all(_fits(t, PLOTLY_HOVER_PX, SIZES["tick"]) for t in texts)
    assert not [t for t in texts if ".:" in t]
    assert PLOTLY_HOVER_PX < PHONE_PLOT_PX


@pytest.mark.parametrize("case", ["tracks_and_stats", "band_map"])
def test_swatches_are_the_static_patch(case):
    # plotly's default 6 px square beside a 20 px label hid the palest levels (2026-10-01).
    import matplotlib as mpl  # noqa: PLC0415

    entries = _legend(FIXTURES[case]().interactive())
    assert {t.marker.size for t in entries} == {PLOTLY_LEGEND_SWATCH_PX}
    assert LEGEND["handleheight"] * on_page(SIZES["legend"]) == PLOTLY_LEGEND_SWATCH_PX
    with mpl.rc_context({"legend.handleheight": 2 * LEGEND["handleheight"]}):
        apply_matplotlib(THEMES["white"])  # the static's patch reads the same value
        assert mpl.rcParams["legend.handleheight"] == LEGEND["handleheight"]


def test_legend_flows_entry_by_entry():
    # Ungrouped, plotly sets every entry at the widest one's width: one entry a line on the
    # 567 px overview, 9 lines for M7 where the static has 3 (2026-10-01).
    fig = demo_heatmap().interactive()
    groups = [t.legendgroup for t in _legend(fig)]
    assert len(set(groups)) == len(groups)
    assert (fig.layout.legend.traceorder, fig.layout.legend.tracegroupgap) == ("grouped", 0)


def test_a_long_legend_entry_wraps_on_a_phone_only():
    """A track level as long as a tissue's full name: a desktop's legend prints it on one line,
    as the static does, and a phone's wraps it to the plot beside its swatch, its height
    growing by the lines. Mutation: the phone's legend left unwrapped, wider than its plot."""
    data = demo_heatmap()
    long = "bone marrow aspirate from the posterior iliac crest of the donor"
    cells = data.columns_in_order()
    track = Track("tissue", {c: long if i % 2 else "blood" for i, c in enumerate(cells)},
                  {"blood": "#56B4E9", long: "#E69F00"})  # fmt: skip
    data = dataclasses.replace(data, tracks=(track,))
    fig = data.interactive()
    assert [t.name for t in _legend(fig)] == _static(data)["legend"]
    assert f"tissue: {long}" in [t.name for t in _legend(fig)]
    phone = _narrowed(fig)
    width = PHONE_PLOT_PX - PLOTLY_LEGEND_SYMBOL_PX
    names = [t.name for t in _legend(phone)]
    assert all(_fits(name, width, SIZES["legend"]) for name in names)
    assert [_plain(name) for name in names] == _static(data)["legend"]
    added = sum(_lines(name) - 1 for name in names)
    assert added > 0
    assert phone.layout.height - fig.layout.height == pytest.approx(added * PLOTLY_LEGEND_LINE_PX)


def _height(data: HeatmapData, fig, title_lines: int) -> float:
    """A page heatmap's height as its parts add up: the legend's lines as a phone stacks
    them, the column labels and title below, then the rows, capped as the static caps them."""
    legend = sum(_lines(t.name) for t in _legend(fig))
    top = max(PLOTLY_MARGIN_PX["t"], legend * PLOTLY_LEGEND_LINE_PX + PLOTLY_LEGEND_PAD_PX)
    labels = max(text_width(c, on_page(SIZES["tick"])) for c in data.columns_in_order())
    title = title_lines * on_page(SIZES["axis_label"]) * LINE_SPACING
    edges = fig.data[0].y
    rows = edges[-1] - edges[0]
    if len(data.rows_in_order()) > HEAT["max_rows_drawn"]:
        rows = HEAT["max_rows_drawn"]
    expected = PLOTLY_HEAT_BASE_PX + top + (labels if fig.layout.xaxis.showticklabels else 0)
    return expected + title + (edges[0] + rows) * PLOTLY_HEAT_ROW_PX


@pytest.mark.parametrize("case", [*FIXTURES, "long_labels"])
def test_height_is_the_sum_of_its_parts(case):
    # Order of operations in ``_height``, on a desktop's layout and on a phone's. A column
    # label is measured: at the average glyph it fell 38 px short of 19-character run names
    # in capitals (2026-10-01). A phone's layout makes room for the axis title wrapped to half
    # its plot, where the page script wraps it, and exists only where some text wraps.
    data = _band_map(LONG) if case == "long_labels" else FIXTURES[case]()
    fig = data.interactive()
    title = _lines(fig.layout.xaxis.title.text)
    assert fig.layout.height == pytest.approx(_height(data, fig, title))
    assert sum(_lines(t.name) for t in _legend(fig)) >= len(_legend(fig))
    assert (fig.layout.meta is not None) == (case == "long_labels")
    if fig.layout.meta is not None:
        lines = wrap_to(data.column_label, PHONE_PLOT_PX / 2, on_page(SIZES["axis_label"]))
        phone = _narrowed(fig)
        assert phone.layout.height == pytest.approx(_height(data, phone, len(lines)))
        assert phone.layout.height > fig.layout.height


def test_row_pitch_holds_a_label_line():
    # The static's pitch per point of label, at the label's page size: 18 px held a 20 px
    # label on every other row only.
    per_point = HEAT["row_mm"] / MM_PER_INCH * 72 / SIZES["tick"]
    assert pytest.approx(per_point * on_page(SIZES["tick"])) == PLOTLY_HEAT_ROW_PX
    assert on_page(SIZES["tick"]) * LINE_SPACING <= PLOTLY_HEAT_ROW_PX


def test_a_desktop_keeps_the_static_labels_and_a_phone_wraps_them():
    """Work order item 2 (Rytis 2026-10-01): the heatmap fix wrapped every label at a phone's
    width on every screen, as plotly cannot change a text by width. A desktop now keeps the
    static's row labels, legend and axis title; narrower than they need, the page swaps in a
    phone's (``meta.narrow``): every line fits half a phone's plot, measured. Mutations: the
    narrow labels at a desktop's width; the desktop's labels wrapped for a phone; the legend
    left unwrapped on a phone."""
    data = _band_map(LONG)
    fig = data.interactive()
    static = _static(data)
    assert list(fig.layout.yaxis.ticktext) == static["left"]
    assert [t.name for t in _legend(fig)] == static["legend"]
    assert fig.layout.xaxis.title.text == data.column_label
    below = fig.layout.meta["narrow"]["below"]
    widest = max(text_width(t, on_page(SIZES["tick"])) for t in static["left"])
    # The static's labels and a line's width per upright column label, inside a desktop's plot.
    assert widest + len(data.columns_in_order()) * on_page(SIZES["tick"]) < below < 686
    phone = _narrowed(fig)
    assert all(_fits(t, PHONE_PLOT_PX / 2, SIZES["tick"]) for t in phone.layout.yaxis.ticktext)
    assert [_plain(t) for t in phone.layout.yaxis.ticktext] == static["left"]
    assert any("<br>" in t for t in phone.layout.yaxis.ticktext)
    width = PHONE_PLOT_PX - PLOTLY_LEGEND_SYMBOL_PX
    assert all(_fits(t.name, width, SIZES["legend"]) for t in _legend(phone))
    legend = fig.layout.legend
    # From the div's top left: at the plot's, a phone's labels pushed it off the right edge.
    assert (legend.xref, legend.yref, legend.x, legend.y, legend.yanchor) == (
        "container", "container", 0, 1, "top")  # fmt: skip
    assert legend.maxheight == 1


def test_phone_plot_is_what_the_stylesheet_leaves():
    # PHONE_PLOT_PX was a second copy of the phone's padding: the stylesheet's could change
    # and the labels still wrap to the old width (a 293 px plot against 317, 2026-10-01).
    css = (ROOT / "src" / "ctx" / "report" / "templates" / "_page.css.j2").read_text()
    narrow = css.split("@media (max-width:{{ layout.narrow }}px){", 1)[1].split("\n}\n", 1)[0]
    rules = dict(re.findall(r"^([^{\n]+)\{([^}]*)\}$", narrow, re.MULTILINE))
    card = re.search(r"^\.fig\{([^}]*)\}$", css, re.MULTILINE).group(1)
    values = dict(re.findall(r"(--[\w-]+): ([^;]+);", css_layout()))

    def px(rule: str, pattern: str) -> int:
        name = re.search(pattern, rule).group(1)
        return int(values[name].removesuffix("px"))

    gutter = px(rules[".page"], r"padding-left:var\((--[\w-]+)\)")
    padding = px(rules[".fig"], r"padding:var\(--[\w-]+\) var\((--[\w-]+)\)")
    border = px(card, r"border:var\((--[\w-]+)\) solid")
    assert PHONE_PX <= NARROW_PX  # the narrow rules are the phone's
    assert PHONE_PX - 2 * (gutter + padding + border) == PHONE_PLOT_PX == 317


def _sizes(node):
    """Every number under a ``size`` key, anywhere in a plotly JSON tree."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "size" and isinstance(value, int | float):
                yield value
            else:
                yield from _sizes(value)
    elif isinstance(node, list):
        for value in node:
            yield from _sizes(value)


@pytest.mark.parametrize("case", list(FIXTURES))
def test_every_size_comes_from_style(case):
    # A font size typed into the figure passed every test: the literal check read only the
    # keyword form, size=12 (2026-10-01).
    figure = FIXTURES[case]().interactive().to_plotly_json()
    figure["layout"].pop("template")
    sizes = set(_sizes(figure))
    assert sizes
    assert sizes <= {*(on_page(v) for v in SIZES.values()), PLOTLY_LEGEND_SWATCH_PX}
