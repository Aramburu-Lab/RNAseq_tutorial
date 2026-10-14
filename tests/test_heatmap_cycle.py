"""A held-square heatmap whose margins alternate with its height, drawn in headless Chromium.

Step-1's 8-sample laptop heatmap, held square, never finished loading: plotly's margins at
946.5 px held the cells at 938.1 px, and at 938.1 px held them at 946.5 px, and the page
script's guard stopped a repeat but not that alternation (C1, 2026-10-09). The fixture is that
heatmap's rho, cell types and batches with the samples renamed (``fixtures/heatmap_cycle.csv``).
Skipped where no Chromium is found; set ``CTX_CHROMIUM``.
"""

from datetime import date
from pathlib import Path

import polars as pl
import pytest
from test_export_browser import CHROMIUM, _browse

from ctx.io.naming import stem
from ctx.quant.figures import correlation_heatmap
from ctx.report.builder import (
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    write,
)
from ctx.viz.save import figure_block, write_figure

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
FIXTURE = Path(__file__).parent / "fixtures" / "heatmap_cycle.csv"
# Each heatmap as drawn once the page has settled: its height, the heights plotly drew in the
# 3 s after (a page that still sets heights draws more), and a cell's height over its width,
# read off the cells' image. The theme toggle is pressed first when ``toggle`` is true.
WATCH = """new Promise(function (done) {
  if (%s) document.getElementById("ctx-theme-toggle").click();
  var d = Array.from(document.querySelectorAll(".js-plotly-plot")).filter(function (p) {
    return p.data.some(function (t) { return t.type === "heatmap"; });
  })[0];
  var drawn = 0;
  setTimeout(function () {
    d.on("plotly_afterplot", function () { drawn += 1; });
    setTimeout(function () {
      var img = d.querySelector("g.hm image").getBoundingClientRect(), z = d.data[0].z;
      done({theme: document.documentElement.getAttribute("data-theme"),
            height: d._fullLayout.height, drawn: drawn,
            aspect: (img.height / z.length) / (img.width / z[0].length)});
    }, 3000);
  }, 1500);
})"""


@pytest.fixture(scope="module")
def page(tmp_path_factory) -> Path:
    """The module page holding the fixture's heatmap, as M4 builds it."""
    out = tmp_path_factory.mktemp("cycle")
    table = pl.read_csv(FIXTURE)
    rho = table.select(pl.col("^rho_.*$")).to_numpy()
    data = correlation_heatmap(rho, table.select("sample", "cell_type", "batch"), "Made-up.")
    assert data._square()
    figure_stem = stem("quant", "cycle_correlation", 1, day=date(2026, 10, 14))
    text = FigureText("cycle", reading=("Each cell is one sample pair's rho.",))
    write_figure(data, out, figure_stem, title="cycle", data_releases={}, text=text)
    built = ModulePage(
        title="M4 Demo: cycle", purpose="Purpose.", input="Input.", method="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Figures", (figure_block(data, figure_stem, text),)),),
        interpretation="Interpretation.", questions=(SelfCheck("A question.", "An answer."),),
        pitfalls=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    return write(built, out / f"{stem('quant', 'cycle_overview', 1, day=date(2026, 10, 14))}.html")


def test_a_height_cycle_settles_square_at_every_width_in_both_themes(page, tmp_path):
    """The page loads, stops setting heights, and draws its cells square within 3 %, at
    1,280 and 375 px, in both themes, and again after the window narrows and widens.

    Measured 2026-10-14: 946.5 px tall, cells 1.014 at 1,280 px and 1.017 at 375 px.
    Mutation: the guard stopping a repeat only (``heightSet !== height``, as before): the
    page never fires its load event, and the driver times out.
    """
    tasks = [{"name": f"{w}_{toggle}", "url": page.as_uri(), "width": w, "height": 900,
              "act": WATCH % ("true" if toggle else "false"),
              "resize": [{"width": other, "act": WATCH % "false"},
                         {"width": w, "act": WATCH % "false"}]}
             for w, other, toggle in ((1280, 375, False), (375, 1280, True))]  # fmt: skip
    ran = _browse(tasks, tmp_path)
    themes = set()
    for name, r in ran.items():
        assert r["errors"] == [], name
        for drawn in (r["acted"], *r["resized"]):
            themes.add(drawn["theme"])
            assert drawn["drawn"] == 0, (name, drawn)
            assert abs(drawn["aspect"] - 1) < 0.03, (name, drawn)
    assert themes == {"white", "black"}
