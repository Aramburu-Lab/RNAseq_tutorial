# ruff: noqa: E501 - page scripts and mutations are quoted whole, a line each
"""Package W, P.1: what a burst of keystrokes in the live-calls box costs on a built M5 page.

Per page, in headless Chromium at 1,280 px: the points the volcano draws; one redraw (a
change event, timed as it runs); then five values typed into the threshold box as input
events 100 ms apart, counting the redraws (the figure's ``ctx-calls`` events) and timing from
the first key to the last redraw, and how long each key's input handler held the page.
Run: python p1_burst.py <page.html> [...]   (CTX_CHROMIUM or a chromium on PATH)
"""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from test_export_browser import _browse  # noqa: E402

BURST = """(async () => {
  await new Promise(r => setTimeout(r, 4000));
  const form = document.querySelector('form.calls:not(.select)');
  const fig = form.closest('.fig'), d = fig.querySelector('.js-plotly-plot');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const drawn = Object.values(cfg.traces).reduce((n, i) => n + d.data[i].x.length, 0);
  const redraws = [];
  fig.addEventListener('ctx-calls', () => redraws.push(performance.now()));
  form.elements.fc.value = '0.8';
  let t0 = performance.now();
  form.elements.fc.dispatchEvent(new Event('change', {bubbles: true}));
  const one = performance.now() - t0;
  redraws.length = 0;
  const held = [], values = ['1', '1.', '1.2', '1.25', '1.5'];
  const start = performance.now();
  for (const v of values) {
    form.elements.fc.value = v;
    t0 = performance.now();
    form.elements.fc.dispatchEvent(new Event('input', {bubbles: true}));
    held.push(Math.round(performance.now() - t0));
    await new Promise(r => setTimeout(r, 100));
  }
  await new Promise(r => setTimeout(r, 3000));
  const up = fig.querySelector('.calls-summary').textContent;
  return {drawn, one_redraw_ms: Math.round(one), keys: values.length, redraws: redraws.length,
          first_key_to_last_redraw_ms: redraws.length ? Math.round(redraws[redraws.length - 1] - start) : null,
          handler_ms: held, threshold_in_force: form.elements.fc.value, summary: up};
})()"""


def main() -> None:
    """Each page twice, so a cold first run shows apart from a warm one."""
    out = {}
    for name in sys.argv[1:]:
        page = Path(name).resolve()
        tasks = [{"name": f"burst_{k}", "url": page.as_uri(), "width": 1280, "height": 900,
                  "act": BURST} for k in range(2)]  # fmt: skip
        with tempfile.TemporaryDirectory(prefix="w_p1_") as downloads:
            results = _browse(tasks, Path(downloads))
        out[f"{page.parents[2].name}/{page.name}"] = {k: {"measured": v["acted"], "errors": v["errors"]}
                          for k, v in results.items()}  # fmt: skip
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
