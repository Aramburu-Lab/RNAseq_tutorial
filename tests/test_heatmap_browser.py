"""The page heatmap measured in headless Chromium (work order item 2, package D, 2026-10-07).

Kit v6's made-up M7 page labelled every other one of its 21 columns, as plotly thins a
categorical axis's labels to the width it draws at (D12); M12's confounding map drew its cells
0.34 high per unit of width at 1,440 px where its files draw 0.44 (D13); a phone's layout was
sized for an x title wrapped to half its plot and a right axis whose title took a whole line,
neither what the page draws (D16). The pages are built as ``test_export_browser.py`` builds
its, run through ``tests/browser.mjs``, and read off the drawn page: the boxes of the labels,
the cells' image, plotly's margins. Skipped where no Chromium is found; set ``CTX_CHROMIUM``.
"""

import itertools
from datetime import date

import numpy as np
import polars as pl
import pytest
from synth import demo_heatmap
from test_export_browser import CHROMIUM, _browse

from ctx.io.naming import stem
from ctx.report.builder import (
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    write,
)
from ctx.viz.heatmap import HeatmapData, Swatch, Track, diverging_levels
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import BAND_COLOURS, DATA_COLOURS, HEAT_ALPHAS, ROLES, THEMES

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
# Rytis, 2026-10-01: opened alone, a figure is drawn twice the page's size; stated apart from
# style.ALONE_SCALE.
TWICE = 2
# The windows a desktop reader opens a page in, a module page's plot 582 and 686 px wide; and
# the 375 px phone every page is checked at.
DESKTOPS = (1280, 1440)
PHONE = 375
# Every heatmap plot of a page as drawn: its size and margins, the cells' image, the texts of
# its labels and titles with their lines and boxes, its box and scroll hint, and its narrow
# layout. Read once plotly has drawn and the page script has fitted the plot to its box.
PLOTS = """new Promise(function (done) { setTimeout(function () {
  function box(e) { var r = e.getBoundingClientRect(); return [r.left, r.top, r.right, r.bottom]; }
  done(Array.from(document.querySelectorAll(".js-plotly-plot")).filter(function (d) {
    return d.data.some(function (t) { return t.type === "heatmap"; });
  }).map(function (d) {
    var full = d._fullLayout, x = full.xaxis, holder = d.closest(".fig-plot");
    var texts = function (s) {
      return Array.from(d.querySelectorAll(s)).filter(function (t) { return t.textContent; })
        .map(function (t) { return {text: t.textContent, box: box(t),
                                    lines: t.querySelectorAll("tspan.line").length || 1}; });
    };
    return {width: full.width, size: full._size, x: [x._offset, x._length],
            cells: box(d.querySelector("g.hm image")), columns: d.data[0].x.length,
            rows: d.data[0].z.length, edges: Array.from(d.data[0].y),
            svg: box(d.querySelector("svg.main-svg")), xticks: texts(".xtick text"),
            yticks: texts(".ytick text"), xtitle: texts(".xtitle"),
            annotations: texts(".annotation-text"), box: holder.clientWidth,
            hint: !holder.querySelector(".scroll-hint").hidden,
            meta: d.layout.meta || {}};
  }));
}, 500); })"""


def _confounding() -> HeatmapData:
    """M12's map: Cramér's V between four design columns, binned at Cohen's steps."""
    variables = ("condition", "cell_type", "batch", "sex")
    v = {("condition", "cell_type"): 0.05, ("condition", "batch"): 0.2,
         ("condition", "sex"): 0.15, ("cell_type", "batch"): 0.46, ("cell_type", "sex"): 0.03,
         ("batch", "sex"): 0.07}  # fmt: skip
    names = ["< 0.1", "0.1 to 0.3", "0.3 to 0.5", "≥ 0.5"]
    rows = []
    for a, b in itertools.product(variables, variables):
        value = None if a == b else v.get((a, b), v.get((b, a)))
        level = None if value is None else names[int(np.searchsorted((0.1, 0.3, 0.5), value,
                                                                     "right"))]  # fmt: skip
        rows.append((a, b, value, level))
    frame = pl.DataFrame(rows, schema=["row_variable", "column_variable", "v", "level"],
                         orient="row")  # fmt: skip
    shades = zip(names[1:], HEAT_ALPHAS, strict=True)
    return HeatmapData(
        frame=frame, row="row_variable", column="column_variable", value="v", level="level",
        columns={"row_variable": ("none", "Variable."), "column_variable": ("none", "Variable."),
                 "v": ("dimensionless", "Cramér's V."), "level": ("none", "Its bin.")},
        row_label="Variable", column_label="Variable",
        swatches=(Swatch(names[0], None), *(Swatch(n, ROLES["comparison"], a) for n, a in shades)),
        caption="Made-up sample sheet: Cramér's V from each pair's contingency table over 33 "
        "samples, without bias correction; the diagonal is blank. Descriptive, no test.",
        row_order=variables, column_order=variables,
    )  # fmt: skip


def _correlation() -> HeatmapData:
    """M4's map: the correlation of 33 made-up samples, two tracks above; a tall heatmap."""
    rng = np.random.default_rng(4)
    samples = [f"S{i:02d}" for i in range(1, 34)]
    profiles = rng.normal(size=(33, 200))
    profiles[:16] += rng.normal(size=200)
    rho = np.corrcoef(profiles)
    flat = [None if i == j else float(rho[i, j] * 4) for i in range(33) for j in range(33)]
    levels, swatches = diverging_levels(flat)
    frame = pl.DataFrame({"row": [s for s in samples for _ in samples], "column": samples * 33,
                          "rho": flat, "level": levels}, strict=False)  # fmt: skip
    kinds = (("cell type", ("monocyte", "neutrophil"), DATA_COLOURS[:2]),
             ("batch", ("short", "long"), (DATA_COLOURS[2], DATA_COLOURS[6])))  # fmt: skip
    tracks = tuple(Track(name, {s: values[i % 2] for i, s in enumerate(samples)},
                         dict(zip(values, colours, strict=True)))
                   for name, values, colours in kinds)  # fmt: skip
    return HeatmapData(
        frame=frame, row="row", column="column", value="rho", level="level",
        columns={"row": ("none", "Sample."), "column": ("none", "Sample."),
                 "rho": ("dimensionless", "Pearson r, times 4."), "level": ("none", "Its bin.")},
        row_label="Sample", column_label="Sample", swatches=swatches,
        caption="Made up: Pearson r of 200 values between 33 samples. Descriptive, no test.",
        row_order=tuple(samples), column_order=tuple(samples), tracks=tracks,
    )  # fmt: skip


def _phone() -> HeatmapData:
    """Band-map rows that wrap on a phone, a statistic at right and a long axis title.

    Its phone layout shows at 375 px, the axis title wrapped there by the page script.
    """
    metrics = ["Reads written by Trim Galore", "Per-base sequence content", "Q30 bases"]
    runs = [f"R{i:02d}" for i in range(4)]
    frame = pl.DataFrame([(r, m, "good") for r in runs for m in metrics], orient="row",
                         schema=["run", "metric", "band"]).with_columns(
        value=pl.lit(None, pl.Float64))  # fmt: skip
    return HeatmapData(
        frame=frame, row="metric", column="run", value="value", level="band",
        columns={"run": ("none", "Run."), "metric": ("none", "Metric."),
                 "band": ("none", "Band."), "value": ("none", "Not used.")},
        row_label="Registry metric.", column_label="Sample, as nf-core/rnaseq names it.",
        swatches=(Swatch("good", BAND_COLOURS["good"]), Swatch("no band", None)),
        caption="Made up.", row_order=tuple(metrics), column_order=tuple(runs),
        row_stats={m: 0.5 + i / 10 for i, m in enumerate(metrics)},
        row_stats_label="Share of runs", row_stats_format=".2f",
    )  # fmt: skip


MAPS = {"labels": demo_heatmap, "confounding": _confounding, "correlation": _correlation,
        "phone": _phone}  # fmt: skip


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Each map's files and Open alone page, and a module page holding them all, in order."""
    out = tmp_path_factory.mktemp("heatmaps")
    data = {name: make() for name, make in MAPS.items()}
    blocks = []
    for name, map_ in data.items():
        figure_stem = stem("qc", f"{name}_demo", 1, day=date(2026, 10, 7))
        text = FigureText(name, reading=("Each cell is one made-up value.",))
        write_figure(map_, out, figure_stem, title=name, data_releases={}, text=text)
        blocks.append(figure_block(map_, figure_stem, text))
    page = ModulePage(
        title="M2 Demo: heatmaps", question="Purpose.", concept="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Figures", tuple(blocks)),), practice="Interpretation.",
        think=(SelfCheck("A question.", "An answer."),), mistakes=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    module = write(page, out / f"{stem('qc', 'heatmaps_overview', 1, day=date(2026, 10, 7))}.html")
    alone = {name: out / f"{b.stem}.html" for name, b in zip(data, blocks, strict=True)}
    return data, module, alone


@pytest.fixture(scope="module")
def ran(built, tmp_path_factory):
    """Every page at every window it is checked at, by name: the module page's maps in order."""
    _, module, alone = built
    tasks = [{"name": f"page_{w}", "url": module.as_uri(), "width": w, "height": 900,
              "act": PLOTS} for w in (*DESKTOPS, PHONE)]  # fmt: skip
    tasks += [{"name": f"{name}_alone_{w}", "url": path.as_uri(), "width": w, "height": 900,
               "act": PLOTS}
              for name, path in alone.items() for w in (*DESKTOPS, 2560)]  # fmt: skip
    return _browse(tasks, tmp_path_factory.mktemp("downloads"))


def _plot(ran: dict, name: str, where: str, width: int) -> dict:
    """One map's plot as drawn: on the module page, or opened alone."""
    if where == "alone":
        (plot,) = ran[f"{name}_alone_{width}"]["acted"]
        return plot
    return ran[f"page_{width}"]["acted"][list(MAPS).index(name)]


def _intersect(a: list[float], b: list[float]) -> bool:
    """Whether two boxes (left, top, right, bottom) share more than an edge."""
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def _static_cells(data: HeatmapData) -> tuple[float, float]:
    """One cell of the figure files, width and height in inches, off the drawn axes."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415

    fig = data.static(THEMES["white"])
    renderer = FigureCanvasAgg(fig).get_renderer()
    fig.draw(renderer)
    box = fig.axes[0].get_window_extent(renderer)
    rows, columns = data._grid()[0].shape
    return box.width / fig.dpi / columns, box.height / fig.dpi / rows


def _page_cells(plot: dict) -> tuple[float, float]:
    """One cell of the page's plot, width and height in px, off the drawn image of the cells."""
    left, top, right, bottom = plot["cells"]
    return (right - left) / plot["columns"], (bottom - top) / plot["rows"]


@pytest.mark.parametrize("width", DESKTOPS)
@pytest.mark.parametrize("where", ["page", "alone"])
def test_every_column_the_static_labels_is_drawn_apart(ran, built, where, width):
    """D12: all 24 columns of an M7-like map labelled as the files label them, no two labels'
    boxes overlapping, on the module page and opened alone at 1,280 and 1,440 px; where the
    window is narrower than the labels need, the plot keeps their width and scrolls under the
    hint. Mutations: the column axis's ticks left to plotly (five numbers for 24 names at
    1,440 px on the page); the least width dropped (neighbouring boxes overlapped at
    1,280 px)."""
    data = built[0]["labels"]
    plot = _plot(ran, "labels", where, width)
    fig = data.static(THEMES["white"])
    labelled = [t.get_text() for t in fig.axes[0].get_xticklabels() if t.get_text()]
    assert len(labelled) == len(data.columns_in_order()) == 24
    assert [t["text"] for t in plot["xticks"]] == labelled
    boxes = sorted((t["box"] for t in plot["xticks"]), key=lambda b: b[0])
    assert not [(a, b) for a, b in itertools.pairwise(boxes) if _intersect(a, b)]
    scale = TWICE if where == "alone" else 1
    ((least, fixed),) = data.interactive().layout.meta["least_width"]
    assert plot["width"] >= scale * least + fixed - 0.5
    assert plot["hint"] == (plot["width"] > plot["box"] + 1)


@pytest.mark.parametrize("name", ["confounding", "correlation"])
@pytest.mark.parametrize(("where", "width"), [("page", 1280), ("page", 1440), ("alone", 1280),
                                              ("alone", 1440), ("alone", 2560)])  # fmt: skip
def test_cells_are_drawn_as_the_files_draw_them(ran, built, name, where, width):
    """D13: a cell's height over its width on the page, off its drawn image, equals the files',
    off their drawn axes, within 2 %, for M12's 4 x 4 map and a 33 x 33 one with two tracks, at
    every desktop width, on the module page and opened alone. Mutation: the rows at the fixed
    pitch, PLOTLY_HEAT_ROW_PX, the files' aspect unread (the 4 x 4 map's cells 0.20 at
    1,440 px against the files' 0.44)."""
    data = built[0][name]
    file_width, file_height = _static_cells(data)
    page_width, page_height = _page_cells(_plot(ran, name, where, width))
    assert page_height / page_width == pytest.approx(file_height / file_width, rel=0.02)


@pytest.mark.parametrize("name", ["confounding", "correlation"])
def test_at_a_module_pages_width_the_figure_has_no_room_to_spare(ran, name):
    """The height is the cells at the files' aspect, the labels, legend and title as plotly sets
    them at a module page's 686 px: under the lowest text the figure ends within its own bottom
    margin and plotly's padding, 16 px. Mutation: the legend budgeted one entry a line, as a
    phone stacks it (the 4 x 4 map's figure ended 89 px under its title)."""
    plot = _plot(ran, name, "page", 1440)
    texts = [*plot["xticks"], *plot["xtitle"], *plot["annotations"]]
    lowest = max(t["box"][3] for t in texts)
    assert 0 <= plot["svg"][3] - lowest <= 16


def test_rows_hold_their_labels_where_the_files_aspect_would_not(ran, built):
    """Where the files' cells would put two 20 px row labels on each other, the rows hold them:
    each band a label's box high, so no two labels' boxes overlap, at every width it is drawn
    at. The M7-like map's files draw cells 0.98 high per unit of width, which at its 716 px
    gives rows of 22.0 px for labels 22.2 px tall. Mutation: no floor (neighbouring labels'
    boxes overlapped by 0.2 px)."""
    data = built[0]["labels"]
    file_width, file_height = _static_cells(data)
    for where, width in itertools.product(("page", "alone"), DESKTOPS):
        plot = _plot(ran, "labels", where, width)
        boxes = [t["box"] for t in plot["yticks"]]
        assert len(boxes) == len(data.rows_in_order()) + len(data.tracks)
        assert not [(a, b) for a, b in itertools.pairwise(boxes) if _intersect(a, b)]
        page_width, page_height = _page_cells(plot)
        assert page_height / page_width >= file_height / file_width


def test_a_phones_layout_is_sized_as_the_page_draws_it(ran, built):
    """D16, at 375 px by device-metrics emulation: the page script wraps the x title to twice
    the room from the plot's centre to the nearer edge, and the phone's layout predicted the
    same lines; plotly's right margin, from the cells to the figure's edge, is ``_right_px``
    within 1 px. Mutations: the title wrapped to half a phone's plot (3 lines for 2); the right
    axis's title counted as a whole line on top of the template's margin (12 px wide)."""
    data = built[0]["phone"]
    plot = _plot(ran, "phone", "page", PHONE)
    narrow = plot["meta"]["narrow"]
    assert plot["edges"] == narrow["restyle"][0][0]["y"][0]  # the phone's layout is drawn
    (title,) = plot["xtitle"]
    assert title["lines"] == narrow["relayout"]["xaxis.title.text"].count("<br>") + 1 == 2
    left, length = plot["x"]
    assert plot["width"] - left - length == pytest.approx(data._right_px(), abs=1)
    assert plot["size"]["r"] == pytest.approx(data._right_px(), abs=1)
