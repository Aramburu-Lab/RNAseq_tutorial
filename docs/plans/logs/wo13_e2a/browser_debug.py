"""Package wo13-e2a: run one script in a page in headless Chromium and print what it returns.

python browser_debug.py <page.html> <script file> [width]
"""

import json
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(WT / "tests"))
from test_export_browser import _browse  # noqa: E402

page, script = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).read_text()
width = int(sys.argv[3]) if len(sys.argv) > 3 else 1280  # noqa: PLR2004
found = _browse([{"name": "d", "url": page.as_uri(), "width": width, "height": 900,
                  "act": script}], WT / "outputs" / "scratch_roots")  # fmt: skip
print(json.dumps({"acted": found["d"]["acted"], "errors": found["d"]["errors"]}, indent=1)[:6000])
