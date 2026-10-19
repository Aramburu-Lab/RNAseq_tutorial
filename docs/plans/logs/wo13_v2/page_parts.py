"""The bytes each script element of a page holds, largest first: python page_parts.py <page>."""

import re
import sys
from pathlib import Path

html = Path(sys.argv[1]).read_text(encoding="utf-8")
parts = re.findall(r"<script([^>]*)>(.*?)</script>", html, flags=re.S)
sized = sorted(((len(body.encode()), attrs.strip()[:80]) for attrs, body in parts), reverse=True)
print(f"page {len(html.encode()):,} B; {len(parts)} scripts")
for size, attrs in sized[:12]:
    print(f"{size:>12,}  {attrs}")
