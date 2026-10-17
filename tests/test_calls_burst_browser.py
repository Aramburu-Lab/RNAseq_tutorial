"""P.1 in headless Chromium: typing into the live-calls box redraws once per burst, and a
change redraws at once.

The demo volcano's page; five values typed 100 ms apart. The mutation caught: the input
event wired to the redraw itself, as before P.1 (five redraws). Skipped where no Chromium is
found; set ``CTX_CHROMIUM``.
"""

from datetime import date

import pytest
from synth import counts, demo_volcano
from test_export_browser import CHROMIUM, _browse, _module

from ctx.io.naming import stem
from ctx.report.builder import FigureText, write
from ctx.viz.save import figure_block

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
BURST = """(async () => {
  const form = document.querySelector("form.calls"), fig = form.closest(".fig");
  const d = fig.querySelector(".js-plotly-plot");
  const cfg = JSON.parse(document.getElementById(form.getAttribute("data-controls")).textContent);
  const up = () => d.data[cfg.traces[cfg.levels.up]].x.length;
  let redraws = 0;
  fig.addEventListener("ctx-calls", () => { redraws += 1; });
  for (const v of ["1", "1.", "1.5", "1.52", "1.5"]) {
    form.elements.fc.value = v;
    form.elements.fc.dispatchEvent(new Event("input", {bubbles: true}));
    await new Promise(r => setTimeout(r, 100));
  }
  const during = redraws;
  await new Promise(r => setTimeout(r, 1000));
  const burst = {redraws: redraws, during: during, up: up()};
  form.elements.fc.value = "1.5";
  form.elements.fc.dispatchEvent(new Event("change", {bubbles: true}));
  return {burst: burst, change: {redraws: redraws - burst.redraws, up: up()}};
})()"""


def test_a_burst_of_keys_redraws_once_and_lands_on_the_last_value(tmp_path):
    data = demo_volcano(counts(), width="double")
    figure_stem = stem("qc", "volcano_demo", 1, day=date(2026, 10, 1))
    block = figure_block(data, figure_stem, FigureText(title="volcano"))
    overview = stem("qc", "burst_overview", 1, day=date(2026, 10, 1))
    page = write(_module([block]), tmp_path / f"{overview}.html")
    got = _browse([{"name": "burst", "url": page.as_uri(), "act": BURST}], tmp_path)["burst"]
    assert got["errors"] == []
    burst, change = got["acted"]["burst"], got["acted"]["change"]
    assert (burst["during"], burst["redraws"]) == (0, 1)  # none while typing, one after
    assert change["redraws"] == 1 and change["up"] == burst["up"]
