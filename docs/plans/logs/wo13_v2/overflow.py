# ruff: noqa: E501 - page scripts and reports are quoted whole
"""What runs past a page's box: python overflow.py <page> <width> [scheme] [downloads].

Loads the page, waits, and lists every element whose right edge passes the page box, widest
first, with its tag, class and the nearest figure's title.
"""

import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402

page, width = Path(sys.argv[1]).absolute(), int(sys.argv[2])
scheme = sys.argv[3] if len(sys.argv) > 3 else "light"  # noqa: PLR2004
downloads = Path(sys.argv[4] if len(sys.argv) > 4 else "/tmp").absolute()  # noqa: PLR2004, S108
ACT = """(async () => {
  await new Promise(r => setTimeout(r, 4000));
  const box = document.documentElement.clientWidth, out = [];
  for (const e of document.querySelectorAll('body *')) {
    const r = e.getBoundingClientRect();
    if (r.right > box + 1 && r.width > 0) {
      const fig = e.closest('.fig');
      out.push([Math.round(r.right), e.tagName, (e.className && e.className.baseVal !== undefined
        ? e.className.baseVal : e.className) || '', fig ? (fig.querySelector('.fig-title') || {}).textContent : '']);
    }
  }
  out.sort((a, b) => b[0] - a[0]);
  return {box: box, scroll: document.documentElement.scrollWidth, items: out.slice(0, 12)};
})()"""
task = {"name": "o", "url": page.as_uri(), "width": width, "height": 900, "scheme": scheme,
        "act": ACT}  # fmt: skip
acted = _browse([task], downloads)["o"]["acted"]
print(f"box {acted['box']} scroll {acted['scroll']}")
for item in acted["items"]:
    print(item)
