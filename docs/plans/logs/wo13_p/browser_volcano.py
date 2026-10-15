"""Package P, GS.3: the unthinned volcano on built M5 pages in headless Chromium.

A real box and lasso drag at 1,280 and 375 px in both themes (tests/browser.mjs, the picker
test's own steps), then per page: the points drawn against the rows embedded, the time one
live-calls redraw takes, and whether the page scrolls sideways.
Run: python browser_volcano.py <page.html> [...]   (CTX_CHROMIUM or a chromium on PATH)
"""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
import test_pick_browser  # noqa: E402
from test_export_browser import _browse  # noqa: E402
from test_pick_browser import BOX, LASSO, _drag  # noqa: E402

# A root built without M6 gives its picker no onward link; read none as none.
test_pick_browser.READ = test_pick_browser.READ.replace(
    "link.getAttribute('href')", "link && link.getAttribute('href')"
).replace("link: link.textContent", "link: link && link.textContent")

# Points drawn, rows embedded, and five live-calls redraws timed by typing a threshold.
MEASURE = """(async () => {
  await new Promise(r => setTimeout(r, 4000));
  const fig = document.querySelector('form.select').closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const form = fig.querySelector('form.calls:not(.select)');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const rows = JSON.parse(document.getElementById(form.getAttribute('data-rows')).textContent);
  const drawn = Object.values(cfg.traces).reduce((n, i) => n + d.data[i].x.length, 0);
  const times = [];
  for (const v of ['1.5', '0.5', '2', '1.2', '1']) {
    form.elements.fc.value = v;
    const t0 = performance.now();
    form.elements.fc.dispatchEvent(new Event('input', {bubbles: true}));
    times.push(Math.round(performance.now() - t0));
    await new Promise(r => setTimeout(r, 200));
  }
  return {drawn, rows: rows.data.length, n_tested: cfg.n_tested, type: d.data[0].type,
          redraw_ms: times, scroll: document.documentElement.scrollWidth
            <= document.documentElement.clientWidth};
})()"""

out = {}
for name in sys.argv[1:]:
    page = Path(name).resolve()
    tasks = [_drag(f"{kind}_{w}{'_dark' if dark else ''}", page, w, drag, shape, dark)
             for kind, drag, shape in (("box", "select", BOX), ("lasso", "lasso", LASSO))
             for w in (1280, 375) for dark in (False, True)]  # fmt: skip
    tasks.append({"name": "measure", "url": page.as_uri(), "width": 1280, "height": 900,
                  "act": MEASURE})  # fmt: skip
    with tempfile.TemporaryDirectory(prefix="p_gs3_") as downloads:
        results = _browse(tasks, Path(downloads))
    summary = {"measure": results["measure"]["acted"], "errors": {}}
    for task, seen in results.items():
        if seen["errors"]:
            summary["errors"][task] = seen["errors"]
        if task == "measure":
            continue
        got = seen["moused"][-1]
        inside, edge, picked = set(got["want"]), set(got["edge"]), set(got["got"])
        right = bool(inside) and inside <= picked <= inside | edge
        summary[task] = {"inside": len(inside), "edge": len(edge), "picked": len(picked),
                         "picks_what_is_inside": right, "count": got["count"],
                         "scroll": got["scroll"], "font": got["font"]}  # fmt: skip
    out[page.name] = summary
print(json.dumps(out, indent=1))
