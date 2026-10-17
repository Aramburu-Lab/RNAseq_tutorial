"""Package M9b: a built M9 page's bytes, cards kept and left out, and what each card part costs.

Run: python bytes.py <page.html> [<page.html> ...]
"""

import re
import sys
from pathlib import Path

from ctx.report.builder import packed_bytes, unpacked

for name in sys.argv[1:]:
    path = Path(name)
    html = path.read_text()
    packed = packed_bytes(html)
    whole = unpacked(html)
    left = re.search(r"(\d[\d,]*) of the (\d[\d,]*) listed genes have no card here", whole)
    kept_titles = re.findall(r'<div class="group" id="group-\d+">\s*<h3>([^<]*)</h3>', whole)
    views = len(re.findall(r'"kind": ?"view"', whole))
    print(
        f"{path.name}: {path.stat().st_size:,} B on disk; packed cards {len(packed)}, "
        f"median {sorted(packed)[len(packed) // 2] if packed else 0:,} B, "
        f"sum {sum(packed):,} B; groups {len(kept_titles)}; "
        f"left out {left.group(1) if left else 0} of {left.group(2) if left else '-'}; "
        f"unpacked {len(whole.encode()):,} B; view controls {views}"
    )
