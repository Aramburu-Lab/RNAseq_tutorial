"""Package P: where two builds of one page differ, as up to 20 short line-level hunks.

Run: python page_diff.py <page A> <page B>
"""

import difflib
import sys
from pathlib import Path

WIDTH = 160
a, b = (Path(p).read_text().splitlines() for p in sys.argv[1:3])
hunks = [d for d in difflib.unified_diff(a, b, lineterm="", n=0) if d[:1] in "+-"]
for line in hunks[2:42]:
    print(line[:WIDTH])
