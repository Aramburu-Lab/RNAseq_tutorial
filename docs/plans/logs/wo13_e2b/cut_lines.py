"""Cut a line range out of a file, after checking its first and last lines (one edit, asserted).

python cut_lines.py <file> <first> <last> <first line starts with> <line after last starts with>
"""

import sys
from pathlib import Path

path, first, last, head, after = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), *sys.argv[4:6]
lines = Path(path).read_text().splitlines(keepends=True)
assert lines[first - 1].startswith(head), lines[first - 1]
assert lines[last].startswith(after), lines[last]
Path(path).write_text("".join(lines[: first - 1] + lines[last:]))
print(f"cut {last - first + 1} lines from {path}")
