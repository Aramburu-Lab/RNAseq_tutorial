"""wo13-y item 1: every function V2 or L changed that the figures package names.

The diff runs from F1's base 37efebc to main before F1's merge, fbc9ac8^1.
Usage: python moved_names.py <diff file>. Prints, per name, the figures files naming it and
the signature lines the diff removed and added, so each call can be read against the module's.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

diff = Path(sys.argv[1]).read_text().splitlines()
signs: dict[str, list[str]] = defaultdict(list)
for line in diff:
    if m := re.match(r"^([+-])\s*def (\w+)", line):
        signs[m.group(2)].append(line)
figures = sorted(Path("src/ctx/figures").glob("*.py"))
for name, lines in sorted(signs.items()):
    users = [p.name for p in figures if re.search(rf"\b{name}\b", p.read_text())]
    if users:
        print(name, users)
        for line in lines:
            print("   ", line.strip()[:110])
