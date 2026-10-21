# ruff: noqa: D103, E501, PLC0415, T201 - a proof record, run once and kept as run
"""wo13-z5: every page built before and after packing draws the same figures, rows and plotly.js.

Run: python compare.py <before dir> <after dir> [--browser] [--timing N]. Pages are paired by
their path under each folder. Per page:

1. Bytes before and after. Each figure's three drawing arguments and each rows table as text:
   the after page's packed text, inflated, is byte for byte the before page's raw text; so is
   plotly.js, against the before page's raw library script.
   M10's before page already packed its figures (Z.1): they are unpacked first.
2. The whole page: the after page unpacked (ctx.report.figures.unpacked) against the before page
   unpacked too, with its unpack script taken out (only M10's had one), line by line, each
   match of a --same REGEX (the two scratch roots' paths) read as one; every differing line is
   printed (the Reproduce panel's wall time and memory differ between any two runs, and the
   checksum of a parquet file that records the path of the tree it was built from).
3. With --browser, in Chromium, both pages from disk (file://): every plot's data and layout as
   plotly holds them once drawn, its full data's types and lengths, every rows table as the page
   parses it, and every console error.
4. With --timing N, the five largest drawing pages N times each, before and after alternating,
   in copies with one script added first in <head>: the time from navigation start to the first
   plot drawn (a .js-plotly-plot in the DOM), to the last plot drawn at load, and to the load
   event, in ms (performance.now()).
"""

import base64
import difflib
import json
import re
import statistics
import sys
import tempfile
import zlib
from pathlib import Path

TESTS = Path(__file__).resolve().parents[4] / "tests"
PACKED_PLOT = re.compile(r'(Plotly\.newPlot\(\s*"([^"]+)",\s*)\.\.\.ctxUnpack\("([^"]*)"\)')
PACKED_ROWS = re.compile(r'<script type="application/octet-stream" id="(ctx-data-\d+)">([^<]*)<')
PACKED_LIB = re.compile(r'<script>ctxUnpack\.script\("([^"]*)"\);</script>')
UNPACK = re.compile(r"<script>\n/\* A packed figure's JSON.*?</script>\n", re.S)
RAW_LIB = re.compile(r"</style>\n<script>(.*?)</script>\n", re.S)
RAW_PLOT = re.compile(r'Plotly\.newPlot\(\s*"([^"]+)",\s*(?=\[)')
RAW_ROWS = re.compile(r'<script type="application/json" id="(ctx-data-\d+)">([^<]*)<')


def inflate(packed: str) -> str:
    return zlib.decompress(base64.b64decode(packed), -15).decode()


def raw_arguments(html: str) -> dict[str, str]:
    decoder, out = json.JSONDecoder(), {}
    for match in RAW_PLOT.finditer(html):
        at = match.end()
        for i in range(3):
            if i:
                at = re.compile(r"\s*,\s*").match(html, at).end()
            _, at = decoder.raw_decode(html, at)
        out[match.group(1)] = html[match.end() : at]
    return out


def text_checks(before: str, after: str) -> tuple[bool, str]:
    figs_b = raw_arguments(before)
    figs_a = {m[2]: inflate(m[3])[1:-1] for m in PACKED_PLOT.finditer(after)}
    rows_b = dict(RAW_ROWS.findall(before))
    rows_a = {i: inflate(p) for i, p in PACKED_ROWS.findall(after)}
    lib_b = RAW_LIB.search(before)[1] if RAW_LIB.search(before) else ""
    lib_m = PACKED_LIB.search(after)
    lib_a = inflate(lib_m[1]) if lib_m else ""
    left = len(raw_arguments(after)) + len(RAW_ROWS.findall(after))
    ok = figs_b == figs_a and rows_b == rows_a and lib_b == lib_a and left == 0
    line = (
        f"figures {len(figs_b)}/{len(figs_a)} equal {figs_b == figs_a}; rows {len(rows_b)}/"
        f"{len(rows_a)} equal {rows_b == rows_a}; plotly.js {len(lib_b.encode()):,} B equal "
        f"{lib_b == lib_a}; raw left {left}"
    )
    return ok, line


def page_diff(before: str, after: str, same: list[str]) -> list[str]:
    from ctx.report.figures import unpacked

    before, after = unpacked(UNPACK.sub("", before, count=1)), unpacked(after)
    for pattern in same:
        before, after = re.sub(pattern, "<scratch>", before), re.sub(pattern, "<scratch>", after)
    return [
        line
        for line in difflib.unified_diff(before.splitlines(), after.splitlines(), lineterm="", n=0)
        if not line.startswith(("---", "+++", "@@"))
    ]


DUMP = """(async () => {
  await new Promise(r => setTimeout(r, 2500));
  const plots = Array.from(document.querySelectorAll('.js-plotly-plot')).map(d => ({
    id: d.id, data: JSON.stringify(d.data), layout: JSON.stringify(d.layout),
    full: JSON.stringify(d._fullData.map(t => [t.type, t.name, t.x && t.x.length,
                                               t.y && t.y.length, t.visible]))}));
  const rows = Array.from(document.querySelectorAll('script[id^="ctx-data-"]')).map(s =>
    JSON.stringify(s.type === 'application/octet-stream' ? ctxUnpack(s.textContent)
                                                          : JSON.parse(s.textContent)));
  return {plots, rows, plotly: window.Plotly ? Plotly.version : null};
})()"""

PROBE = """<script>window.__z5 = {plots: 0};
new MutationObserver(function () {
  var n = document.querySelectorAll('.js-plotly-plot').length;
  if (n && !__z5.first) __z5.first = performance.now();
  if (n > __z5.plots) { __z5.plots = n; __z5.last = performance.now(); }
}).observe(document.documentElement, {subtree: true, childList: true, attributes: true,
                                      attributeFilter: ['class']});
addEventListener('load', function () { __z5.load = performance.now(); });</script>"""
TIMED = """(async () => { await new Promise(r => setTimeout(r, 300)); return window.__z5; })()"""


def browse(tasks: list[dict]) -> dict[str, dict]:
    sys.path.insert(0, str(TESTS))
    from test_export_browser import _browse

    with tempfile.TemporaryDirectory() as folder:
        return _browse(tasks, Path(folder), timeout=3600)


def drawn(pairs: list[tuple[str, Path, Path]]) -> bool:
    tasks = [
        {
            "name": f"{w}:{name}",
            "url": p.resolve().as_uri(),
            "width": 1280,
            "height": 900,
            "act": DUMP,
        }
        for name, b, a in pairs
        for w, p in (("before", b), ("after", a))
    ]
    seen, ok = browse(tasks), True
    print("\n3. Chromium, from disk")
    for name, _, _ in pairs:
        b, a = seen[f"before:{name}"], seen[f"after:{name}"]
        same = b["acted"] == a["acted"]
        plots = b["acted"]["plots"]
        each = all(x == y for x, y in zip(plots, a["acted"]["plots"], strict=False))
        good = same and each and b["errors"] == [] and a["errors"] == []
        ok = ok and good
        print(
            f"   {'ok  ' if good else 'FAIL'} {name}: plots {len(plots)}/{len(a['acted']['plots'])}, "
            f"rows {len(b['acted']['rows'])}/{len(a['acted']['rows'])}, data+layout+full+rows "
            f"equal {same}, plotly {b['acted']['plotly']}/{a['acted']['plotly']}, errors "
            f"{len(b['errors'])}/{len(a['errors'])} {(b['errors'] + a['errors'])[:2]}"
        )
    return ok


def timing(pairs: list[tuple[str, Path, Path]], n: int, scratch: Path) -> None:
    pairs = sorted(pairs, key=lambda t: -t[2].stat().st_size)[:5]
    tasks = []
    for name, b, a in pairs:
        for w, p in (("before", b), ("after", a)):
            copy = scratch / w / name
            copy.parent.mkdir(parents=True, exist_ok=True)
            html = p.read_text()
            copy.write_text(html.replace("<head>\n", "<head>\n" + PROBE + "\n", 1))
    for k in range(n):
        for name, _, _ in pairs:
            for w in ("before", "after") if k % 2 == 0 else ("after", "before"):
                tasks.append(
                    {
                        "name": f"{w}:{name}:{k}",
                        "act": TIMED,
                        "width": 1280,
                        "height": 900,
                        "url": (scratch / w / name).resolve().as_uri(),
                    }
                )
    seen = browse(tasks)
    print(f"\n4. Time from navigation start, ms, median of {n} (min-max), each page from disk")
    print(
        "   page | bytes before | bytes after | first plot before | after | last plot at load before | after | load before | after"
    )
    for name, b, a in pairs:
        cols = []
        for key in ("first", "last", "load"):
            for w in ("before", "after"):
                v = [seen[f"{w}:{name}:{k}"]["acted"][key] for k in range(n)]
                cols.append(f"{statistics.median(v):.0f} ({min(v):.0f}-{max(v):.0f})")
        print(f"   {name} | {b.stat().st_size:,} | {a.stat().st_size:,} | " + " | ".join(cols))


def main(before_dir: str, after_dir: str, browser: bool, runs: int, same: list[str]) -> int:
    from ctx.report.figures import unpacked

    bd, ad = Path(before_dir), Path(after_dir)
    names = sorted(
        {str(p.relative_to(bd)) for p in bd.rglob("*.html")}
        | {str(p.relative_to(ad)) for p in ad.rglob("*.html")}
    )
    pairs, ok, total_b, total_a = [], True, 0, 0
    print("1-2. page | bytes before | bytes after | saved | text checks | differing lines")
    for name in names:
        b, a = bd / name, ad / name
        if not (b.is_file() and a.is_file()):
            print(f"   MISSING {name}: before {b.is_file()}, after {a.is_file()}")
            ok = False
            continue
        before, after = b.read_text(), a.read_text()
        good, line = text_checks(unpacked(before), after)
        diff = page_diff(before, after, same)
        sb, sa = len(before.encode()), len(after.encode())
        total_b, total_a = total_b + sb, total_a + sa
        print(
            f"   {'ok  ' if good else 'FAIL'} {name} | {sb:,} | {sa:,} | {sb - sa:,} | {line} | {len(diff)}"
        )
        for d in diff:
            print(f"        {d[:200]}")
        ok = ok and good
        if "Plotly.newPlot" in after:
            pairs.append((name, b, a))
    print(
        f"   all {len(names)} pages: {total_b:,} B before, {total_a:,} B after; {len(pairs)} draw"
    )
    if browser:
        ok = drawn(pairs) and ok
    if runs:
        with tempfile.TemporaryDirectory() as folder:
            timing(pairs, runs, Path(folder))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    args = sys.argv[1:]
    runs = int(args[args.index("--timing") + 1]) if "--timing" in args else 0
    same = [args[i + 1] for i, a in enumerate(args) if a == "--same"]
    sys.exit(main(args[0], args[1], "--browser" in args, runs, same))
