"""wo13-z: every tab of every Reproduce panel on main's pages against this package's.

Run: python compare_panels.py. Main's page holds each tab's pasted block and its script as
base64; this package's holds step 1 per tab and the shared steps once, joined as the page's
script joins them (tests/test_figure_panel.py's panels()). The two sides' roots differ only by
their scratch folder, which is replaced before comparing.
"""

import base64
import html
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tests"))
from test_figure_panel import panels  # noqa: E402

OLD = re.compile(r'<code id="ctx-repro-(\d+)-([a-z]+)">(.*?)</code>'
                 r'.*?href="data:text/x-shellscript;base64,([^"]+)"', re.S)  # fmt: skip
PACKED = re.compile(r'<script type="application/octet-stream" id="ctx-packed-\d+">([^<]*)</script>')
SCRATCH = ROOT / "outputs" / "scratch_z"
SIDE = "scratch_z/main/"  # main's root folder, in the pasted text and inside the base64


def unpack(text: str) -> str:
    """The page with its packed groups opened."""
    import gzip  # noqa: PLC0415

    return PACKED.sub(lambda m: gzip.decompress(base64.b64decode(m[1])).decode(), text)


def old_panels(text: str, tab: str) -> list[tuple[str, str]]:
    """(pasted, script) of one tab of every panel on a page as main wrote it."""
    return [(html.unescape(p), base64.b64decode(s).decode().replace(SIDE, "scratch_z/z/"))
            for _, key, p, s in OLD.findall(text) if key == tab]  # fmt: skip


total = equal = 0
for new_page in sorted(SCRATCH.glob("z/*/outputs/[0-9][0-9]_*/*_overview_v1.html")):
    old_page = SCRATCH / "main" / new_page.relative_to(SCRATCH / "z")
    if not old_page.is_file():
        continue
    old_text = unpack(old_page.read_text()).replace("scratch_z/main/", "scratch_z/z/")
    new_text = unpack(new_page.read_text())
    for tab in ("docker", "apptainer", "pixi"):
        before, after = old_panels(old_text, tab), panels(new_text, tab)
        same = sum(a == b for a, b in zip(before, after, strict=False))
        if "--diff" in sys.argv and same != len(before):
            import difflib  # noqa: PLC0415

            a, b = next((a, b) for a, b in zip(before, after, strict=False) if a != b)
            for x, y, what in ((a[0], b[0], "pasted"), (a[1], b[1], "script")):
                print(what, "".join(difflib.unified_diff(x.splitlines(True), y.splitlines(True))))
            sys.exit(0)
        total += len(before)
        equal += same
        if same != len(before) or len(before) != len(after):
            print(f"{new_page.name} {tab}: {len(before)} main, {len(after)} Z, {same} equal")
print(f"tabs compared {total}; byte-equal in pasted block and script {equal}")
