"""Package wo13-l: copy built pages into outputs/wo13_l_try with one CSS text replaced.

It tries a CSS rule in the browser before rebuilding; the copies are never shipped.

Usage: python try_css.py <old> <new> <page.html> [...]
"""

import sys
from pathlib import Path

old, new, *pages = sys.argv[1:]
out = Path(__file__).absolute().parents[4] / "outputs" / "wo13_l_try"
out.mkdir(parents=True, exist_ok=True)
for page in map(Path, pages):
    html = page.read_text(encoding="utf-8")
    assert html.count(old) == 1, (page.name, html.count(old))
    (out / page.name).write_text(html.replace(old, new), encoding="utf-8")
    print(out / page.name)
