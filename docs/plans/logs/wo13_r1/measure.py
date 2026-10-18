"""wo13-r1 proof: per module of a ctx root, its pages, their bytes, and each page's own Reproduce.

Wall time and peak RSS are read from the page's Reproduce block (td.wall, td.rss-self,
td.rss-children), as the module measured them; a module's pages share one run, so its overview
page's numbers are its run's. Usage: measure.py <ctx root>
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

CELL = re.compile(r'<td class="(wall|rss-self|rss-children)">([^<]*)</td>')


def main() -> None:
    """Print one row per module folder: pages, bytes, wall, peak RSS self and children."""
    root = Path(sys.argv[1]).expanduser()
    pages = defaultdict(list)
    for page in sorted(root.glob("outputs/*/*.html")):
        pages[page.parent.name].append(page)
    print("folder\tpages\tbytes\tlargest page bytes\twall\tpeak RSS self\tpeak RSS children")
    for folder, found in sorted(pages.items()):
        # The overview, or the hub, holds the module's Reproduce block; figure pages do not.
        first = next((p for p in found if "overview" in p.name or "hub" in p.name), found[0])
        cells = dict(CELL.findall(first.read_text(errors="ignore")))
        sizes = [p.stat().st_size for p in found]
        print("\t".join((folder, str(len(found)), f"{sum(sizes):,}", f"{max(sizes):,}",
                         cells.get("wall", "none"), cells.get("rss-self", "none"),
                         cells.get("rss-children", "none"))))  # fmt: skip


if __name__ == "__main__":
    main()
