# ruff: noqa: B905 - a proof record, run once and kept as run
"""Each differing line pair of every page, at the first differing character."""

import difflib
import re
import sys
from pathlib import Path

from ctx.report.figures import unpacked

UNPACK = re.compile(r"<script>\n/\* A packed figure's JSON.*?</script>\n", re.S)
bd, ad = Path(sys.argv[1]), Path(sys.argv[2])
pat = sys.argv[3]
for b in sorted(bd.rglob("*.html")):
    a = ad / b.relative_to(bd)
    x = (
        unpacked(UNPACK.sub("", b.read_text(), count=1))
        .replace("z5/base_before", "z5/base_after")
        .replace("z5/before/", "z5/after/")
    )
    y = unpacked(a.read_text())
    xl, yl = x.splitlines(), y.splitlines()
    sm = difflib.SequenceMatcher(None, xl, yl, autojunk=False)
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op != "replace":
            continue
        for p, q in zip(xl[i1:i2], yl[j1:j2]):
            if pat not in p:
                continue
            k = next((k for k in range(min(len(p), len(q))) if p[k] != q[k]), min(len(p), len(q)))
            print(b.relative_to(bd), len(p), len(q), k)
            print("  -", p[max(0, k - 120) : k + 120])
            print("  +", q[max(0, k - 120) : k + 120])
