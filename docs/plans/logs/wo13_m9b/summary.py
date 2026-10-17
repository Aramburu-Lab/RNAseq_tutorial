"""Package M9b: one line per task of a browser.py report.

Run: python summary.py <report.json>
"""

import json
import sys
from pathlib import Path

report = json.loads(Path(sys.argv[1]).read_text())
for name, checks in report.items():
    short = {k: v for k, v in checks.items() if k not in ("first_view", "caption", "summary")}
    print(name, short)
first = next(iter(report.values()))
print("first view:", first.get("first_view"))
print("caption:", first.get("caption"))
