"""wo13-f: which modules wrote a page for which dataset, in each real ctx root (read only)."""

import re
from pathlib import Path

PAGE = re.compile(
    r"\d{4}-\d{2}-\d{2}_(?P<module>[a-z_]+?)_(?P<dataset>[a-z0-9_]+)_(overview|hub)_v\d+\.html"
)
for root in sorted(Path.home().glob("ctx_*")):
    if not (root / "outputs").is_dir():
        continue
    found: dict[str, set[str]] = {}
    for page in (root / "outputs").glob("*/*.html"):
        if m := PAGE.fullmatch(page.name):
            found.setdefault(m["dataset"], set()).add(page.parent.name)
    print(f"== {root.name}")
    for dataset, modules in sorted(found.items()):
        print(f"  {dataset}: {' '.join(sorted(modules))}")
