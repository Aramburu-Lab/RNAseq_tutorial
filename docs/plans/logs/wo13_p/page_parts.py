"""Package P: bytes of a built page by part: each large script, each plotly figure's traces.

Run: python page_parts.py <page.html> [...]
"""

import json
import re
import sys
from pathlib import Path

SMALL = 20000

for name in sys.argv[1:]:
    html = Path(name).read_text()
    print(f"{Path(name).name}: page {len(html.encode()):,} B")
    for m in re.finditer(r"<script([^>]*)>(.*?)</script>", html, re.S):
        attrs, body = m.group(1), m.group(2)
        if len(body.encode()) < SMALL:
            continue
        print(f"  script {len(body.encode()):>12,} B  {attrs.strip()[:70]}")
        plot = re.search(r'Plotly\.newPlot\(\s*"[^"]+",\s*(\[.*?\]),\s*\{"', body, re.S)
        if not plot:
            continue
        for i, t in enumerate(json.loads(plot.group(1))):
            size, points = len(json.dumps(t).encode()), len(t.get("x") or [])
            label = str(t.get("name", ""))[:30]
            print(f"    trace {i} {t.get('type', 'scatter')} {label!r}: {points:,} pts {size:,} B")
