"""wo13-z: a page's bytes by part, each byte counted in one part only.

Run: python parts.py <page.html> [...]. Parts, in the order they are claimed: packed groups
(M9's cards, gzip and base64), Reproduce panels (split into their pasted blocks, their base64
downloads as main writes them, and the rest), plotly's figure scripts, the controls JSON (views
and their options), the other JSON blocks, the library scripts (plotly.js, jsPDF), the page
script, tables, CSS; the rest is text and markup. A packed group's panels are counted inside the
packed group, and again, unpacked, on a line of their own.
"""

import base64
import gzip
import re
import sys
from pathlib import Path

PACKED = re.compile(r'<script type="application/octet-stream" id="ctx-packed-\d+">([^<]*)</script>')
REPRO = re.compile(r'<details class="repro"[^>]*>.*?<details class="repro-upstream">.*?</details>'
                   r"\n</details>", re.S)  # fmt: skip
CODE = re.compile(r'<code id="ctx-repro-\d+-[a-z]+"[^>]*>.*?</code>', re.S)
HREF = re.compile(r'href="data:text/x-shellscript;base64,[^"]*"')
PARTS = [
    ("figure data", re.compile(r'<script type="application/json" id="ctx-data-\d+">.*?</script>',
                               re.S)),
    ("plotly figure scripts", re.compile(r"<script>\s*window\.PLOTLYENV.*?</script>", re.S)),
    ("controls and views", re.compile(r'<script type="application/json" id="ctx-controls-\d+">'
                                      r".*?</script>", re.S)),
    ("other JSON", re.compile(r'<script type="application/json"[^>]*>.*?</script>', re.S)),
    ("libraries", re.compile(r"<script>/\*!? ?(?:plotly|jsPDF|\*).{0,200}?.*?</script>", re.S)),
    ("page script", re.compile(r"<script>.*?</script>", re.S)),
    ("tables", re.compile(r"<table.*?</table>", re.S)),
    ("css", re.compile(r"<style.*?</style>", re.S)),
]  # fmt: skip


def size(text: str) -> int:
    """UTF-8 bytes."""
    return len(text.encode())


def parts(html: str) -> dict[str, int]:
    """Bytes by part; each claimed span is cut out before the next part is looked for."""
    found: dict[str, int] = {"page": size(html)}
    packed = PACKED.findall(html)
    found["packed groups"] = sum(size(m[0]) for m in PACKED.finditer(html))
    inner = "".join(gzip.decompress(base64.b64decode(b)).decode() for b in packed)
    found["  panels inside packed groups, unpacked"] = sum(size(p) for p in REPRO.findall(inner))
    rest = PACKED.sub("", html)
    panels = REPRO.findall(rest)
    found["Reproduce panels"] = sum(size(p) for p in panels)
    found["  panels: pasted blocks"] = sum(size(c) for p in panels for c in CODE.findall(p))
    found["  panels: base64 downloads"] = sum(size(h) for p in panels for h in HREF.findall(p))
    found["  panels: count"] = len(panels)
    rest = REPRO.sub("", rest)
    for name, pattern in PARTS:
        found[name] = sum(size(m) for m in pattern.findall(rest))
        rest = pattern.sub("", rest)
    found["text and markup"] = size(rest)
    return found


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(Path(arg).name)
        for name, value in parts(Path(arg).read_text()).items():
            print(f"  {name:<42} {value:>12,}")
