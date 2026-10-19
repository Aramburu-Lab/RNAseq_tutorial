"""Package wo13-l: evaluate one script in a page at one width and print what it returns.

Usage: <main env python> docs/plans/logs/wo13_l/env.py python docs/plans/logs/wo13_l/probe.py
           <page.html> <width> <script.js> [dark]
"""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parents[4] / "tests"))
from test_export_browser import _browse  # noqa: E402

page, width, script, *dark = sys.argv[1:]
task = {"name": "p", "url": Path(page).absolute().as_uri(), "width": int(width),
        "height": {1280: 900, 375: 812}.get(int(width), 900), "act": Path(script).read_text(),
        "scheme": "dark" if dark else "light"}  # fmt: skip
with tempfile.TemporaryDirectory() as d:
    got = _browse([task], Path(d))["p"]
print(json.dumps(got.get("acted"), indent=1))
print("errors", got.get("errors"))
