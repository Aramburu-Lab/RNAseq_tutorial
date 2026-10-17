"""Bytes of each embedded JSON script on a page, with a view's option count and label.

python page_parts.py <page.html>
"""

import json
import re
import sys
from pathlib import Path

# A script below this many bytes, and of no control kind, is left out of the listing.
SHOWN = 20_000
html = Path(sys.argv[1]).read_text()
print(f"page {len(html.encode()):,} B")
pattern = r'<script type="application/json" id="(ctx-[a-z]+-\d+)">(.*?)</script>'
for m in re.finditer(pattern, html, re.S):
    body = m.group(2)
    size = len(body.encode())
    try:
        found = json.loads(body)
    except json.JSONDecodeError:
        print(m.group(1), f"{size:,}", "not JSON")
        continue
    kind = found.get("kind") if isinstance(found, dict) else None
    extra = ""
    if kind == "view":
        extra = f"{found['label']!r} options={len(found['options'])} grid={'grid' in found}"
    if size > SHOWN or kind:
        print(m.group(1), f"{size:,}", kind or "", extra)
