"""wo13-f1 M10: are two figure pages' survival controls (the live panel) the same?

Usage: wo13_f1_clinical_py.sh python <this file> <a.html> <b.html>
"""

import json
import re
import sys
from pathlib import Path


def survival(path: str) -> dict:
    """The page's survival controls."""
    for body in re.findall(r'id="ctx-controls-\d+">(.*?)</script>', Path(path).read_text(), re.S):
        found = json.loads(body)
        if found["kind"] == "survival":
            return found
    raise SystemExit(f"{path}: no survival controls")


a, b = (survival(p) for p in sys.argv[1:3])
print(
    f"block characters: {len(a['block']):,} and {len(b['block']):,}; genes {a['genes']}, "
    f"{b['genes']}; screen rows {len(a['screen'])}, {len(b['screen'])}"
)
print(
    "survival controls equal"
    if a == b
    else f"DIFFER in {sorted(k for k in a.keys() | b.keys() if a.get(k) != b.get(k))}"
)
