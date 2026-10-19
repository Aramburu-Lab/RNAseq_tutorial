"""V2 E.1e on a real M9 page in Chromium: every box's donor lines against the values table.

python box_real.py <root> <dataset> <downloads dir>: at 1,280 and 375 px in both themes, each
plot that draws a paired box's lines (a lines trace named "pair") is read back; its segments
are counted and matched, end to end, with each donor's two values in M9's own_values table for
the featured contrast; the page box is measured against the window.
"""

import sys
from pathlib import Path

import polars as pl

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402

root, dataset, downloads = Path(sys.argv[1]).absolute(), sys.argv[2], Path(sys.argv[3]).absolute()
CLOSE = 1e-9  # a value read back from plotly, against the table it was drawn from
downloads.mkdir(parents=True, exist_ok=True)
m9 = root / "outputs" / "11_genecard"
(page,) = sorted(m9.glob(f"*_genecard_{dataset}_overview_v1.html"))[-1:]
values = pl.read_parquet(sorted(m9.glob(f"*_{dataset}_own_values_v1.parquet"))[-1])
SIZES = (("1280_light", 1280, "light"), ("1280_dark", 1280, "dark"),
         ("375_light", 375, "light"), ("375_dark", 375, "dark"))  # fmt: skip
ACT = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(2000);
  /* Packed cards draw when they come near the view: bring each there, then wait for all. */
  for (const holder of document.querySelectorAll('.packed-body')) {
    holder.scrollIntoView();
    for (let i = 0; i < 50 && holder.getAttribute('data-state') !== 'drawn'; i++) await wait(100);
  }
  const plots = [];
  for (const d of document.querySelectorAll('.fig-plot')) {
    d.scrollIntoView();
    for (let i = 0; i < 30 && !d.querySelector('.js-plotly-plot .main-svg'); i++) await wait(100);
    const p = d.querySelector('.js-plotly-plot');
    if (p) plots.push(p);
  }
  window.scrollTo(0, 0); await wait(500);
  const out = [];
  for (const d of plots) {
    const line = (d.data || []).find(t => t.mode === 'lines' && t.name === 'pair');
    if (!line) continue;
    const fig = d.closest('.fig');
    const bg = getComputedStyle(d.querySelector('.main-svg') || d).backgroundColor;
    out.push({title: (fig.querySelector('.fig-title') || {}).textContent || '',
              x: Array.from(line.x), y: Array.from(line.y), units: Array.from(line.customdata),
              colour: line.line.color, ticks: (d.layout.xaxis || {}).ticktext || [],
              paths: d.querySelectorAll('.scatterlayer .trace path.js-line').length,
              wide: d.getBoundingClientRect().width});
  }
  return {plots: out, theme: document.documentElement.getAttribute('data-theme'),
          box: document.documentElement.clientWidth, scroll: document.documentElement.scrollWidth};
})()"""
tasks = [{"name": name, "url": page.as_uri(), "width": width, "height": 900, "scheme": scheme,
          "act": ACT} for name, width, scheme in SIZES]  # fmt: skip
seen = _browse(tasks, downloads, timeout=2400)
print(f"page {page.name} {page.stat().st_size:,} B")
for name, _, _ in SIZES:
    acted, errors = seen[name]["acted"], seen[name]["errors"]
    matched = segments = 0
    for plot in acted["plots"]:
        symbol = plot["title"].split(" in ")[0]
        rows = values.filter(pl.col("symbol") == symbol)
        for i in range(0, len(plot["x"]), 3):
            unit, a, b = plot["units"][i], plot["y"][i], plot["y"][i + 1]
            segments += 1
            mine = rows.filter(pl.col("donor") == unit)
            found = {g: mine.filter(pl.col("group") == g)["value"].to_list()
                     for g in ("reference", "tested")}  # fmt: skip
            matched += any(abs(r - a) < CLOSE for r in found["reference"]) and any(
                abs(t - b) < CLOSE for t in found["tested"])  # fmt: skip
    print(f"{name}: theme {acted['theme']}, {len(acted['plots'])} plots with donor lines, "
          f"{segments} lines, {matched} matched to own_values, colours "
          f"{sorted({p['colour'] for p in acted['plots']})}, box {acted['box']} scroll "
          f"{acted['scroll']}, errors {errors}")  # fmt: skip
