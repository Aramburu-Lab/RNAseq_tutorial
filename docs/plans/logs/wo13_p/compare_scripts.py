"""Package P, E.4: the page script of each page built twice, compared byte for byte.

Run: python compare_scripts.py <outputs dir A> <outputs dir B> <relative page path> [...]
Prints, per page, the script's size and whether it is identical, and whether the whole page is.
"""

import hashlib
import re
import sys
from pathlib import Path

MARK = "/* One script for every page kind."


def page_script(html: str) -> str:
    """The script element that holds the page script, with the shared objects before it."""
    for body in re.findall(r"<script>(.*?)</script>", html, re.S):
        if MARK in body:
            return body
    raise SystemExit("no page script")


a_dir, b_dir = Path(sys.argv[1]), Path(sys.argv[2])
for name in sys.argv[3:]:
    a, b = (a_dir / name).read_text(), (b_dir / name).read_text()
    sa, sb = page_script(a), page_script(b)
    digest = hashlib.sha256(sb.encode()).hexdigest()[:16]
    print(
        f"{name}: script {len(sa.encode()):,} / {len(sb.encode()):,} B, identical {sa == sb}, "
        f"sha {digest}; page {len(a.encode()):,} / {len(b.encode()):,} B, identical {a == b}"
    )
