"""M2's fastp grid: settings that keep no read are its floor, and the bands span the rest.

The grid is the one ``ctx qc --reads`` measured on the Step-1 teaching subset's grid run on
2026-10-07 (outputs/02_qc/2026-10-07_qc_Non-cancerous-blood-teaching-subset_grid_v1.csv,
1,000,000 reads in). -l 102 is longer than that run's 76-nt reads (FastQC's average length,
76.0), so its 6 cells keep none; the other 30 keep 78.3537 to 99.9985 %. Eight equal bands
from the grid's 0 put all 30 in the top two bands, 6 and 7.
"""

import json
import math
import os
import shutil
import subprocess
from datetime import date
from pathlib import Path

import polars as pl
import pytest

from ctx.io.naming import stem
from ctx.qc import Fastp
from ctx.qc.figures import flat_grid, grid_figure
from ctx.report import text
from ctx.viz.save import write_figure
from ctx.viz.style import CIVIDIS, THEMES

ROOT = Path(__file__).resolve().parents[1]
Q = (10, 15, 20, 25, 30, 35)
L = (15, 25, 36, 50, 75, 102)
# Percent of reads kept, a row per -q in Q's order, a column per -l in L's order.
KEPT = (
    (99.9985, 99.9983, 99.9385, 99.3935, 97.5137, 0.0),
    (98.9688, 98.9686, 98.9091, 98.3667, 96.4968, 0.0),
    (98.9688, 98.9686, 98.9091, 98.3667, 96.4968, 0.0),
    (98.0675, 98.0673, 98.0086, 97.4703, 95.6129, 0.0),
    (97.0184, 97.0182, 96.9596, 96.4236, 94.5781, 0.0),
    (80.3623, 80.3622, 80.3205, 79.8933, 78.3537, 0.0),
)
GRID = pl.DataFrame(
    [(q, n, round(pct * 10_000), pct) for q, row in zip(Q, KEPT, strict=True)
     for n, pct in zip(L, row, strict=True)],
    schema=["q", "l", "passed_reads", "passed_pct"], orient="row",
)  # fmt: skip
FLOOR = [pct == 0.0 for pct in GRID["passed_pct"]]
LOWEST, HIGHEST = 78.3537, 99.9985
GREY = THEMES["white"].reference


@pytest.fixture(scope="module")
def surface():
    data = grid_figure(GRID, Fastp(), "the grid run")
    assert data is not None
    return data


def test_the_bands_span_the_cells_that_keep_reads(surface):
    """Band 0 starts at the lowest share a setting keeps, 78.3537 %, not at the floor's 0.

    Mutation caught: the edges over every cell, 0 to 99.9985, which put the 30 cells in bands
    6 (5 cells) and 7 (25 cells)."""
    assert (surface.edges()[0], surface.edges()[-1]) == (LOWEST, HIGHEST)
    bands = surface.table()["band"].to_list()
    assert [b is None for b in bands] == FLOOR and sum(FLOOR) == 6
    kept = [b for b in bands if b is not None]
    assert len(kept) == 30 and (min(kept), max(kept)) == (0, len(CIVIDIS) - 1)


def test_a_grid_flat_above_its_floor_is_a_note():
    """A run whose every read has one quality keeps all of its reads or none at each setting:
    nothing lies between to map, which the old edges drew as a two-colour map of eight bands.
    """
    two = GRID.with_columns(passed_pct=pl.when(pl.col("passed_pct") > 0).then(100.0)
                            .otherwise(0.0))  # fmt: skip
    assert grid_figure(two, Fastp(), "the grid run") is None
    assert flat_grid(two, "the grid run") == text.fill(
        "qc", "notes.flat_floor", run="the grid run", kept=30, none=6, pct=100.0)  # fmt: skip
    assert "every setting keeps 0.0 %" in flat_grid(
        GRID.with_columns(passed_pct=pl.lit(0.0)), "the grid run")  # fmt: skip


def _static_fills(surface, theme) -> tuple[object, list[str]]:
    """The figure, and the colour matplotlib fills at each cell by its own rule.

    A band holds lower < z <= upper, and its lower edge too where that is the grid's minimum
    (``ContourSet._get_lowers_and_uppers``).
    """
    from matplotlib.colors import to_hex  # noqa: PLC0415 - heavy
    from matplotlib.contour import ContourSet  # noqa: PLC0415 - heavy

    fig = surface.static(theme)
    filled = [c for c in fig.axes[0].collections if isinstance(c, ContourSet) and c.filled]
    bands = [(lo, hi, to_hex(colour).upper(), path)
             for cs in filled
             for lo, hi, colour, path in zip(*cs._get_lowers_and_uppers(), cs.get_facecolor(),
                                             cs.get_paths(), strict=True)]  # fmt: skip
    # Every band of the eight is drawn with an area, and the floor is one more fill.
    assert [b[2] for b in bands if len(b[3].vertices)] == [*CIVIDIS, theme.reference]
    fills = []
    for z in GRID["passed_pct"]:
        (colour,) = [c for lo, hi, c, _ in bands if lo < z <= hi]
        fills.append(colour)
    return fig, fills


def _plotly_fills(surface) -> list[str]:
    """The colour plotly fills at each cell, from the traces the page draws.

    plotly.js 3 (contour/make_color_map.js) gives band j its colour j only when the bands'
    range over their width floors to exactly eight: one less, from an edge difference an ulp
    wide, dropped a level and coloured bands 3 to 6 one step up (read from the DOM,
    2026-10-07). A constraint trace shades the side its operation rules out.
    """
    fig = surface.interactive()
    bands = fig.data[0]
    (floor,) = [t for t in fig.data if t.type == "contour" and t.contours.type == "constraint"]
    assert math.floor((bands.zmax - bands.zmin) / bands.contours.size) == len(CIVIDIS)
    assert (floor.contours.operation, floor.contours.value) == (">", bands.contours.start)
    steps = [bands.colorscale[2 * j][1] for j in range(len(CIVIDIS))]
    fills = []
    for z in GRID["passed_pct"]:
        j = min(int((z - bands.contours.start) // bands.contours.size), len(CIVIDIS) - 1)
        fills.append(floor.fillcolor if z < floor.contours.value else steps[j])
    return fills


def _expected(fills: list[str], grey: str, palette: tuple[str, ...] = CIVIDIS) -> None:
    """The floor's 6 cells grey; the 30 others in cividis, from band 0 to band 7."""
    assert [f == grey for f in fills] == FLOOR
    kept = [f for f, on_floor in zip(fills, FLOOR, strict=True) if not on_floor]
    assert set(kept) <= set(palette)
    at = dict(zip(GRID["passed_pct"], fills, strict=True))
    assert (at[LOWEST], at[HIGHEST]) == (palette[0], palette[-1])


@pytest.mark.parametrize("theme", list(THEMES))
def test_the_static_draws_the_floor_grey_and_the_bands_over_the_rest(surface, theme):
    """Mutation caught: the edges over every cell (the 30 fill bands 6 and 7 only)."""
    from matplotlib.colors import to_hex  # noqa: PLC0415 - heavy

    fig, fills = _static_fills(surface, THEMES[theme])
    _expected(fills, THEMES[theme].reference)
    legend = fig.axes[0].get_legend()
    keys = dict(zip((t.get_text() for t in legend.get_texts()), legend.legend_handles,
                    strict=True))  # fmt: skip
    key = keys[text.word("qc", "figures.grid.floor")]
    assert to_hex(key.get_facecolor()).upper() == THEMES[theme].reference


def test_the_interactive_draws_the_same(surface):
    """Mutations caught: the edges over every cell; the plotly band width as the difference
    of two edges (7.99999999999999 bands on this grid)."""
    _expected(_plotly_fills(surface), GREY)
    names = {t.name: t for t in surface.interactive().data if t.type == "scatter"}
    assert names[text.word("qc", "figures.grid.floor")].marker.color == GREY


CHROMIUM = os.environ.get("CTX_CHROMIUM") or next(
    filter(None, map(shutil.which, ("chromium", "chromium-browser", "google-chrome"))), None
)
# The fill of every contour band on the page, by trace: (level, colour) in drawing order.
READ_FILLS = """(async () => {
  for (let i = 0; i < 200 && !document.querySelector('g.contour'); i++) {
    await new Promise((r) => setTimeout(r, 50));
  }
  const plot = document.querySelector('.js-plotly-plot');
  const fills = (g) => Array.from(g.querySelectorAll('g.contourfill path'))
    .map((p) => [p.__data__.level, p.style.fill]);
  return {
    traces: Array.from(plot.querySelectorAll('g.contour')).map(fills),
    legend: Array.from(plot.querySelectorAll('.legendtext')).map((t) => t.textContent),
  };
})()"""


def _rgb(colour: str) -> str:
    """A hex colour as the DOM writes a fill."""
    return "rgb({}, {}, {})".format(*(int(colour[i : i + 2], 16) for i in (1, 3, 5)))


@pytest.mark.skipif(CHROMIUM is None or shutil.which("node") is None,
                    reason="no Chromium or node; set CTX_CHROMIUM")  # fmt: skip
def test_the_page_draws_the_floor_grey_and_every_band(surface, tmp_path):
    """The figure page in Chromium: the fill plotly drew at each level, read from the DOM.

    Mutations caught: the edges over every cell; the band width as an edge difference."""
    figure_stem = stem("qc", "grid_demo", 1, day=date(2026, 10, 7))
    write_figure(surface, tmp_path, figure_stem, title="Grid", data_releases={})
    plan = {"chromium": CHROMIUM, "downloads": str(tmp_path),
            "tasks": [{"name": "grid", "url": (tmp_path / f"{figure_stem}.html").as_uri(),
                       "act": READ_FILLS}]}  # fmt: skip
    done = subprocess.run([shutil.which("node"), str(ROOT / "tests" / "browser.mjs")],
                          input=json.dumps(plan), capture_output=True, text=True, timeout=300,
                          check=True)  # fmt: skip
    (result,) = json.loads(done.stdout)
    assert result["errors"] == []
    bands, ((value, grey),) = result["acted"]["traces"][:2]
    palette = tuple(_rgb(c) for c in CIVIDIS)
    # Each of the eight bands drawn, at its edge, in its colour; the floor shaded grey below.
    assert [fill for _, fill in bands] == list(palette)
    assert [level for level, _ in bands] == pytest.approx(list(surface.edges()[:-1]))
    assert (value, grey) == (pytest.approx(LOWEST), _rgb(GREY))
    fills = [grey if z < value else max((lv, f) for lv, f in bands if lv <= z)[1]
             for z in GRID["passed_pct"]]  # fmt: skip
    _expected(fills, grey, palette)
    assert text.word("qc", "figures.grid.floor") in result["acted"]["legend"]
