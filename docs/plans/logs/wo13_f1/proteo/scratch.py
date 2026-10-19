"""A scratch root for ctx proteo, written only under its own outputs/.

The real root's data and resources are linked read only; the newest gather and M5 tables M7
reads are copied, never moved, never linked.

Usage: python scratch.py <real root> <scratch root> <dataset> <gather table>... -- <M5 what>...
"""

import shutil
import sys
from pathlib import Path

from ctx.de.run import newest

real, scratch, dataset = Path(sys.argv[1]).expanduser(), Path(sys.argv[2]), sys.argv[3]
rest = sys.argv[4:]
cut = rest.index("--")
copies = {"01_gather": rest[:cut], "05_de": rest[cut + 1 :]}
if scratch.exists():
    raise SystemExit(f"{scratch} exists; move it aside first")
(scratch / "outputs").mkdir(parents=True)
for name in ("data", "resources"):
    (scratch / name).symlink_to(real / name)
for folder, whats in copies.items():
    (scratch / "outputs" / folder).mkdir()
    for what in whats:
        table = newest(real / "outputs" / folder, dataset, what)
        for part in (table, table.with_suffix(".csv"), table.with_name(f"{table.stem}_README.md")):
            if part.is_file():
                shutil.copy2(part, scratch / "outputs" / folder / part.name)
                print("copied", part)
