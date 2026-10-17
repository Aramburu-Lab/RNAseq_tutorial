"""Evaluate one act on a page in headless Chromium and print what it returns.

python probe.py <page.html> <act file> [width] [scheme]
"""

import json
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402

page, act = Path(sys.argv[1]).absolute(), Path(sys.argv[2]).read_text()
width = int(sys.argv[3]) if sys.argv[3:] else 1280
scheme = sys.argv[4] if sys.argv[4:] else "light"
(seen,) = _browse([{"name": "p", "url": page.as_uri(), "act": act, "width": width,
                    "scheme": scheme}], Path("/tmp")).values()  # fmt: skip
print(json.dumps({"acted": seen["acted"], "errors": seen["errors"]}, indent=1)[:6000])
