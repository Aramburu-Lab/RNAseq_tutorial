# ruff: noqa: E501 - page scripts and mutations are quoted whole, a line each
"""Package W: which elements of a page reach past a 375 px window, every packed card drawn.

Run: python overflow.py <page.html> [...]
"""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from test_export_browser import _browse  # noqa: E402

WIDE = """(async () => {
  for (const b of document.querySelectorAll("[data-unpack]")) b.click();
  await new Promise(r => setTimeout(r, 6000));
  const w = document.documentElement.clientWidth, out = [];
  for (const e of document.querySelectorAll("body *")) {
    const r = e.getBoundingClientRect();
    if (r.right > w + 1 && !e.closest(".table-scroll, .fig-plot, .fig-below, pre, .toc")) {
      out.push([e.tagName, e.className.baseVal === undefined ? e.className : e.className.baseVal,
                Math.round(r.right), (e.textContent || "").slice(0, 60)]);
    }
  }
  return {width: w, scroll: document.documentElement.scrollWidth, wide: out.slice(0, 12)};
})()"""

for name in sys.argv[1:]:
    page = Path(name).absolute()
    with tempfile.TemporaryDirectory() as d:
        got = _browse([{"name": "w", "url": page.as_uri(), "width": 375, "height": 812,
                        "act": WIDE}], Path(d))["w"]  # fmt: skip
    print(page.name, json.dumps(got["acted"], indent=1))
