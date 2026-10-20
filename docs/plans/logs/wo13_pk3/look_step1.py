"""PK3: Step-1's rebuilt M6 opened with 6 genes in its address, in Chromium.

At 1,280 px light and 375 px dark: the live ORA's bars, its summary, the page's width, where
the set labels start, and a PNG of the plot.

Usage: python look_step1.py <url with fragment> <folder for the PNGs and JSON>
"""

import base64
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from test_export_browser import _browse  # noqa: E402 - the tests' own driver

ACT = """(async () => {
  await new Promise(r => setTimeout(r, 4000));
  const form = document.querySelector('form.ora'), fig = form.closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const png = await Plotly.toImage(d, {format: 'png', width: d.offsetWidth,
                                       height: d.offsetHeight});
  return {bars: d.data.filter(t => t.type === 'bar').map(t => [t.name, t.marker.color,
            Array.from(t.y), Array.from(t.x)]),
          summary: fig.querySelector('.ora-summary').textContent,
          notes: fig.querySelector('.ora-notes').textContent,
          top: fig.getBoundingClientRect().top, window: innerHeight,
          labels_left: Math.min(...Array.from(d.querySelectorAll('.ytick text'))
            .map(t => t.getBoundingClientRect().left)),
          plot_left: d.getBoundingClientRect().left,
          box: document.documentElement.clientWidth,
          scroll: document.documentElement.scrollWidth, png};
})()"""

url, out = sys.argv[1], Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
tasks = [{"name": f"{w}_{s}", "url": url, "width": w, "height": h, "scheme": s, "act": ACT}
         for w, h, s in ((1280, 900, "light"), (375, 812, "dark"))]  # fmt: skip
for name, result in _browse(tasks, out).items():
    seen = result["acted"]
    png = seen.pop("png")
    (out / f"{name}.png").write_bytes(base64.b64decode(png.split(",", 1)[1]))
    print(name, json.dumps({**seen, "errors": result["errors"]}, indent=1))
