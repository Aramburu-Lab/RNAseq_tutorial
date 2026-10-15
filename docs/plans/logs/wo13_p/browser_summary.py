"""Package P: one line per page from browser_volcano.py's JSON.

Run: python browser_summary.py <json>
"""

import json
import sys
from pathlib import Path

for page, v in json.loads(Path(sys.argv[1]).read_text()).items():
    m = v["measure"]
    drags = [k for k, x in v.items() if isinstance(x, dict) and "picks_what_is_inside" in x]
    ok = [k for k in drags if v[k]["picks_what_is_inside"] and v[k]["scroll"]]
    picked = sorted({v[k]["picked"] for k in drags})
    print(
        f"{page}: drawn {m['drawn']:,}, rows {m['rows']:,}, tested {m['n_tested']:,}, "
        f"redraw ms {m['redraw_ms']}, no sideways scroll {m['scroll']}, drags right "
        f"{len(ok)} of {len(drags)} (picked {picked}), errors {v['errors'] or 'none'}"
    )
