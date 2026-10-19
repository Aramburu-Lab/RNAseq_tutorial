"""wo13-f1: bytes of a page, and of its Reproduce panels' details blocks.

Usage: wo13_f1_py.sh python docs/plans/logs/wo13_f1/panel_bytes.py <page.html> [...]
"""

import re
import sys
from pathlib import Path

BLOCK = re.compile(r'<details[^>]*class="[^"]*repro[^"]*".*?</details>', re.S)
for name in sys.argv[1:]:
    page = Path(name).read_text()
    blocks = BLOCK.findall(page)
    print(
        f"{Path(name).stat().st_size:>12,} B page; {len(blocks)} panels, "
        f"{sum(len(b.encode()) for b in blocks):,} B; {Path(name).name}"
    )
