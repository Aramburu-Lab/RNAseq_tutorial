"""Package P, GS.3: one live-calls redraw timed on a built M5 page, in SVG and in WebGL.

As built (SVG scatter), then with its point traces switched to WebGL (scattergl) in the
page, to decide whether the point count needs WebGL. Headless Chromium with --disable-gpu,
as tests/browser.mjs runs it, has no WebGL: its WebGL times measure nothing drawn (see log).
Run: python browser_gl_probe.py <page.html> [...]
"""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from test_export_browser import _browse  # noqa: E402

PROBE = """(async () => {
  await new Promise(r => setTimeout(r, 4000));
  const fig = document.querySelector('form.select').closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const form = fig.querySelector('form.calls:not(.select)');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const traces = Object.values(cfg.traces);
  async function timed() {
    const out = [];
    for (const v of ['1.5', '0.5', '2', '1.2', '1']) {
      form.elements.fc.value = v;
      const t0 = performance.now();
      form.elements.fc.dispatchEvent(new Event('input', {bubbles: true}));
      out.push(Math.round(performance.now() - t0));
      await new Promise(r => setTimeout(r, 300));
    }
    return out;
  }
  const svg = await timed();
  Plotly.restyle(d, {type: 'scattergl'}, traces);
  await new Promise(r => setTimeout(r, 1000));
  const gl = await timed();
  const drawn = traces.reduce((n, i) => n + d.data[i].x.length, 0);
  return {drawn, svg_ms: svg, gl_ms: gl, gl_type: d.data[traces[0]].type};
})()"""

for name in sys.argv[1:]:
    page = Path(name).resolve()
    task = {"name": "probe", "url": page.as_uri(), "width": 1280, "height": 900, "act": PROBE}
    with tempfile.TemporaryDirectory(prefix="p_gl_") as downloads:
        seen = _browse([task], Path(downloads))["probe"]
    print(json.dumps({page.name: {**seen["acted"], "errors": seen["errors"]}}))
