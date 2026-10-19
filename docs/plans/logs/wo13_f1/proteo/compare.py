"""Compare a scratch M7 run with the real root's newest run of one day, file by file.

Every file the real run wrote under the same name (date aside) is compared byte for byte;
READMEs and HTML pages carry the date, paths and panels, so they are listed by size alone.
Files only the scratch run wrote are listed with their bytes.

Usage: python compare.py <real 07_proteo> <real day> <scratch 07_proteo>
"""

import sys
from pathlib import Path

real, day, scratch = Path(sys.argv[1]).expanduser(), sys.argv[2], Path(sys.argv[3])
new = {p.name[11:]: p for p in scratch.iterdir() if p.is_file()}
old = {p.name[11:]: p for p in real.glob(f"{day}_*") if p.is_file()}
same = differ = 0
for name in sorted(set(new) & set(old)):
    if name.endswith((".md", ".html")):
        continue
    if new[name].read_bytes() == old[name].read_bytes():
        same += 1
    else:
        differ += 1
        print("DIFFERS", name)
print(f"data, pdf and svg files: {same} byte-equal, {differ} differ")
for name in sorted(set(new) & set(old)):
    if name.endswith(".html"):
        sizes = old[name].stat().st_size, new[name].stat().st_size
        print(f"page {name}: real {sizes[0]:,} B, scratch {sizes[1]:,} B")
added = sorted(set(new) - set(old))
total = sum(new[n].stat().st_size for n in added)
for name in added:
    print(f"new {name}: {new[name].stat().st_size:,} B")
print(f"new files: {len(added)}, {total:,} B")
print("only in the real run:", sorted(set(old) - set(new)))
