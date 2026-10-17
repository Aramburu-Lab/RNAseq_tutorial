"""Package wo13-e2a: which elements of a page reach past a window of the given width.

Elements inside a box that scrolls or clips sideways are left out: they widen no page.
    python browser_wide.py <page.html> <width>
"""

import json
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WT / "tests"))
from test_export_browser import _browse  # noqa: E402

ACT = """(async () => {
  await new Promise(r => setTimeout(r, 3000));
  const w = document.documentElement.clientWidth, out = [];
  const clipped = e => {
    for (let p = e.parentElement; p && p !== document.body; p = p.parentElement) {
      if (['auto', 'scroll', 'hidden', 'clip'].includes(getComputedStyle(p).overflowX)) return true;
    }
    return false;
  };
  document.querySelectorAll('body *').forEach(e => {
    const r = e.getBoundingClientRect();
    if (r.right > w + 1 && e.getClientRects().length && !clipped(e)) {
      const fig = e.closest('figure.fig');
      const css = e.className && e.className.baseVal !== undefined ? e.className.baseVal
                  : e.className;
      out.push([e.tagName, String(css || '').slice(0, 40), Math.round(r.right),
                fig ? (fig.querySelector('.fig-title') || {}).textContent
                    : (e.textContent || '').slice(0, 50)]);
    }
  });
  return {client: w, inner: window.innerWidth, scroll: document.documentElement.scrollWidth,
          wide: out.slice(0, 60)};
})()"""
page, width = Path(sys.argv[1]).resolve(), int(sys.argv[2])
found = _browse([{"name": "wide", "url": page.as_uri(), "width": width, "height": 900,
                  "act": ACT}], Path(sys.argv[1]).parent)  # fmt: skip
print(json.dumps(found["wide"]["acted"], indent=1))
