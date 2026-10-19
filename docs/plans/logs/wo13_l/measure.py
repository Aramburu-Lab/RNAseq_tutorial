# ruff: noqa: E501 - the page script is quoted whole
"""Package wo13-l: each page's document width at 1280 and 375 px, light and dark, in Chromium.

Usage: <main env python> docs/plans/logs/wo13_l/env.py python docs/plans/logs/wo13_l/measure.py
           <out.json> <page.html> [...]

For each page and each (width, scheme) it reports the page's box (documentElement.clientWidth),
its scroll width, and up to 8 elements or text runs reaching past the box that sit in no
sideways scroller of their own (any, widest first, for an out file ending _all.json). It waits
for the page script to settle (plots drawn, heights fitted, packed cards unpacked) first.
"""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parents[4] / "tests"))
from test_export_browser import _browse  # noqa: E402

WIDE = """(async () => { const ALL = %s;
  await new Promise(r => setTimeout(r, 4000));
  /* M9's packed cards draw only when unpacked: measure them drawn. */
  const packed = document.querySelectorAll("[data-unpack]");
  for (const b of packed) b.click();
  if (packed.length) await new Promise(r => setTimeout(r, 6000));
  const w = document.documentElement.clientWidth, out = [];
  const scrolls = (e) => { for (let p = e; p && p !== document.body; p = p.parentElement) {
    const s = getComputedStyle(p).overflowX; if (s === "auto" || s === "scroll" || s === "hidden") return true; } return false; };
  const name = (e) => String(e.className.baseVal === undefined ? e.className : e.className.baseVal).slice(0, 40);
  for (const e of document.querySelectorAll("body *")) {
    const r = e.getBoundingClientRect();
    if (r.width > 0 && r.right > w + 0.5 && (ALL || !scrolls(e.parentElement))) {
      out.push([e.tagName, name(e), e.id, Math.round(r.left), Math.round(r.right), (e.textContent || "").trim().slice(0, 50)]);
    }
  }
  /* Text running past its own box, which no element's box shows. */
  const walk = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  for (let n = walk.nextNode(); n; n = walk.nextNode()) {
    if (!n.textContent.trim() || (!ALL && scrolls(n.parentElement))) continue;
    const range = document.createRange(); range.selectNodeContents(n);
    const r = range.getBoundingClientRect();
    if (r.width > 0 && r.right > w + 0.5) out.push(["#text", name(n.parentElement), "", Math.round(r.left), Math.round(r.right), n.textContent.trim().slice(0, 50)]);
  }
  if (ALL) out.sort((a, b) => b[4] - a[4]);
  return {box: w, scroll: document.documentElement.scrollWidth, wide: out.slice(0, 8)};
})()"""

out, *pages = sys.argv[1:]
WIDE = WIDE % ("true" if out.endswith("_all.json") else "false")
tasks = []
for i, name in enumerate(pages):
    url = Path(name).absolute().as_uri()
    for width in (1280, 375):
        for scheme in ("light", "dark"):
            tasks.append({"name": f"{i}|{width}|{scheme}", "url": url, "width": width,
                          "height": {1280: 900, 375: 812}[width], "scheme": scheme, "act": WIDE})  # fmt: skip
rows = []
with tempfile.TemporaryDirectory() as d:
    for start in range(0, len(tasks), 8):
        got = _browse(tasks[start : start + 8], Path(d), timeout=1800)
        for name, result in got.items():
            i, width, scheme = name.split("|")
            seen = result.get("acted") or {}
            rows.append({"page": pages[int(i)], "width": int(width), "scheme": scheme,
                         "box": seen.get("box"), "scroll": seen.get("scroll"),
                         "wide": seen.get("wide"), "errors": result.get("errors")})  # fmt: skip
            print(
                Path(pages[int(i)]).name,
                width,
                scheme,
                seen.get("box"),
                seen.get("scroll"),
                flush=True,
            )
Path(out).write_text(json.dumps(rows, indent=1) + "\n")
