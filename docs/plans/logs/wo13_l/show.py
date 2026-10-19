"""Package wo13-l: print a measure.py result, the light rows by default.

Usage: python show.py <out.json> [all]
"""

import json
import sys
from pathlib import Path

rows = json.loads(Path(sys.argv[1]).read_text())
for r in rows:
    if r["scheme"] == "light" or sys.argv[2:]:
        print(Path(r["page"]).name, r["width"], r["scheme"], r["box"], r["scroll"])
        for w in r["wide"] or []:
            print("   ", w)
