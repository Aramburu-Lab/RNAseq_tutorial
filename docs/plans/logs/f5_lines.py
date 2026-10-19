"""F.5: every line of the old dataset.py below its imports is in the package once, unchanged.

Usage: ``python docs/plans/logs/f5_lines.py <old dataset.py> src/ctx/dataset``. Compares the
non-blank lines as multisets, so a line copied twice or edited shows, and lists what the
package adds (its docstrings and imports).
"""

import sys
from collections import Counter
from pathlib import Path

old = Path(sys.argv[1]).read_text().splitlines()[14:]
new = [line for f in sorted(Path(sys.argv[2]).glob("*.py")) for line in f.read_text().splitlines()]
before = Counter(line for line in old if line.strip())
after = Counter(line for line in new if line.strip())
lost = before - after
print("old lines missing from the package:", sum(lost.values()), dict(lost))
print("lines the old file had more than once and the package has more often:",
      {k: (before[k], after[k]) for k in before if after[k] > before[k]})  # fmt: skip
print("lines only the package has:", sum((after - before).values()))
for line in sorted(after - before):
    print("  +", line)
