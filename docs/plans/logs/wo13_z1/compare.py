# ruff: noqa: D103, E501, PLC0415 - a proof record, run once and kept as run
"""wo13-z1: an M10 page before and after packing draws the same figures and carries the same rows.

Run: python compare.py <before.html> <after.html> [--browser]. Three checks:

1. Each figure's three drawing arguments and each rows table, as text: the after page's packed
   text, inflated, is byte for byte the before page's raw text.
2. The whole page: the after page with every packed part put back raw, and its unpack script
   taken out, against the before page, line by line; every line that still differs is printed
   (the Reproduce panel's wall time and memory, which differ between any two runs).
3. With --browser, in Chromium: every plot's data and layout as plotly holds them once drawn,
   its full data's types and lengths, and every rows table as the page parses it.
"""

import base64
import difflib
import json
import re
import sys
import zlib
from pathlib import Path

TESTS = Path(__file__).resolve().parents[4] / "tests"
PACKED_PLOT = re.compile(r'(Plotly\.newPlot\(\s*"[^"]+",\s*)\.\.\.ctxUnpack\("([^"]*)"\)')
PACKED_ROWS = re.compile(
    r'<script type="application/octet-stream" (id="ctx-data-\d+")>([^<]*)</script>'
)
UNPACK = re.compile(r"<script>\n/\* A packed figure's JSON.*?</script>\n", re.S)
RAW_PLOT = re.compile(r'Plotly\.newPlot\(\s*"([^"]+)",\s*(?=\[)')
RAW_ROWS = re.compile(r'<script type="application/json" id="(ctx-data-\d+)">([^<]*)</script>')


def inflate(packed: str) -> str:
    return zlib.decompress(base64.b64decode(packed), -15).decode()


def raw_arguments(html: str) -> dict[str, str]:
    """Each figure's three arguments' text, as written, by div id."""
    decoder, out = json.JSONDecoder(), {}
    for match in RAW_PLOT.finditer(html):
        at = match.end()
        for i in range(3):
            if i:
                at = re.compile(r"\s*,\s*").match(html, at).end()
            _, at = decoder.raw_decode(html, at)
        out[match.group(1)] = html[match.end() : at]
    return out


def restored(html: str) -> str:
    """The after page with each packed part raw again and its unpack script out."""
    html = UNPACK.sub("", html, count=1)
    html = PACKED_PLOT.sub(lambda m: m.group(1) + inflate(m.group(2))[1:-1], html)
    return PACKED_ROWS.sub(
        lambda m: f'<script type="application/json" {m.group(1)}>{inflate(m.group(2))}</script>',
        html,
    )


def main(before_path: str, after_path: str, browser: bool) -> int:
    before, after = Path(before_path).read_text(), Path(after_path).read_text()
    print(
        f"before {len(before.encode()):,} B {before_path}\nafter  {len(after.encode()):,} B {after_path}"
    )
    figures_before = raw_arguments(before)
    figures_after = {
        m.group(1).split('"')[1]: inflate(m.group(2))[1:-1] for m in PACKED_PLOT.finditer(after)
    }
    rows_before = dict(RAW_ROWS.findall(before))
    rows_after = {i.split('"')[1]: inflate(p) for i, p in PACKED_ROWS.findall(after)}
    same_figures = figures_before == figures_after
    same_rows = rows_before == rows_after
    print(
        f"1. figures: {len(figures_before)} before, {len(figures_after)} after, arguments equal as text: {same_figures}"
    )
    print(
        f"   rows: {len(rows_before)} before, {len(rows_after)} after, equal as text: {same_rows}"
    )
    print(
        f"   raw figures left on the after page: {len(raw_arguments(after))}; raw rows left: {len(RAW_ROWS.findall(after))}"
    )
    back = restored(after)
    diff = [
        line
        for line in difflib.unified_diff(before.splitlines(), back.splitlines(), lineterm="", n=0)
        if not line.startswith(("---", "+++", "@@"))
    ]
    print(f"2. page restored against before: {len(diff)} differing lines")
    for line in diff:
        print(f"   {line[:220]}")
    ok = same_figures and same_rows and figures_before and rows_before
    if browser:
        ok = drawn(before_path, after_path) and ok
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


DUMP = """(async () => {
  await new Promise(r => setTimeout(r, 4000));
  const plots = Array.from(document.querySelectorAll('.js-plotly-plot')).map(d => ({
    id: d.id, data: JSON.stringify(d.data), layout: JSON.stringify(d.layout),
    full: JSON.stringify(d._fullData.map(t => [t.type, t.name, t.x && t.x.length,
                                               t.y && t.y.length, t.visible]))}));
  const rows = Array.from(document.querySelectorAll('script[id^="ctx-data-"]')).map(s =>
    JSON.stringify(s.type === 'application/octet-stream' ? ctxUnpack(s.textContent)
                                                          : JSON.parse(s.textContent)));
  return {plots, rows};
})()"""


def drawn(before_path: str, after_path: str) -> bool:
    """Every plot and rows table, read in Chromium from each page, compared."""
    import tempfile

    sys.path.insert(0, str(TESTS))
    from test_export_browser import _browse

    tasks = [
        {"name": n, "url": Path(p).resolve().as_uri(), "width": 1280, "height": 900, "act": DUMP}
        for n, p in (("before", before_path), ("after", after_path))
    ]
    with tempfile.TemporaryDirectory() as folder:
        seen = _browse(tasks, Path(folder))
    b, a = seen["before"], seen["after"]
    plots = len(b["acted"]["plots"])
    print(
        f"3. Chromium: errors before {b['errors']}, after {a['errors']}; plots {plots}, {len(a['acted']['plots'])}"
    )
    for x, y in zip(b["acted"]["plots"], a["acted"]["plots"], strict=True):
        print(
            f"   {x['id']}: data {x['data'] == y['data']} ({len(x['data']):,} chars), layout {x['layout'] == y['layout']}, full {x['full'] == y['full']}"
        )
    rows = b["acted"]["rows"] == a["acted"]["rows"]
    print(f"   rows tables {len(b['acted']['rows'])}, equal as parsed: {rows}")
    return b["errors"] == a["errors"] == [] and plots > 0 and b["acted"] == a["acted"]


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2], "--browser" in sys.argv))
