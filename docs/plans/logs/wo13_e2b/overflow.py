"""The page box against its scroll width at 375 px, and the elements reaching past the window.

python overflow.py <page.html> [act]: act, if given, runs before the measurement.
"""

import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402

page = Path(sys.argv[1]).absolute()
act = (
    (sys.argv[2] if sys.argv[2:] else "null")
    + """;
(async () => {
  await new Promise(r => setTimeout(r, 1500));
  const w = document.documentElement.clientWidth;
  const wide = Array.from(document.querySelectorAll('body *')).filter(e => {
    const b = e.getBoundingClientRect(); return b.right > w + 0.5 && b.width > 0;
  }).slice(0, 8).map(e => [e.tagName, e.className, e.id,
                             Math.round(e.getBoundingClientRect().right)]);
  return {box: w, scroll: document.documentElement.scrollWidth, wide: wide};
})()"""
)
(seen,) = _browse([{"name": "o", "url": page.as_uri(), "width": 375, "act": act}], Path("/tmp")
                  ).values()  # fmt: skip
print(page.name, seen["acted"])
