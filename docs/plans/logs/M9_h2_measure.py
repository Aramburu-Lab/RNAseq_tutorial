"""Package M9, H.2: the online lung run's requests, cache bytes, pages and unavailable panels.

Reads ~/ctx_lung only. Requests are the cache manifest's entries by retrieval time, split at
the start of each run's log; a failed request is never cached, so each page's "unavailable"
lines are listed too.
"""

import json
import re
from collections import Counter
from pathlib import Path

from ctx.report.builder import MAX_BYTES, packed_bytes, unpacked

ROOT = Path.home() / "ctx_lung"
LOGS = Path(__file__).parent
manifest = json.loads((ROOT / "resources/cache/manifest.json").read_text())
# The LUSC run met HTTP 500 from Ensembl REST on three lookups and was run once more.
RUNS = ("luad_online", "lusc_online", "lusc_online_retry")
starts = {r: (LOGS / f"M9_h2_{r}.log").read_text().split(" ", 1)[0] for r in RUNS}
print("run starts (UTC):", starts)
for slug, nxt in zip(RUNS, [*list(starts.values())[1:], "9999"], strict=True):
    begin = starts[slug]
    mine = {
        k: v for k, v in manifest.items() if begin <= v["retrieved"].replace("+00:00", "Z") < nxt
    }
    kinds = Counter(k.split("/")[0] for k in mine)
    status = Counter(v["status"] for v in mine.values())
    print(f"{slug}: {len(mine)} requests cached, {sum(v['bytes'] for v in mine.values()):,} B "
          f"bodies; by provider {dict(kinds)}; status {dict(status)}")  # fmt: skip
early = {
    k: v for k, v in manifest.items() if v["retrieved"].replace("+00:00", "Z") < starts[RUNS[0]]
}
print(f"before the runs: {len(early)} ({', '.join(early)})")
files = [p for p in (ROOT / "resources/cache").rglob("*.gz")]
print(f"cache: {len(files)} gz files, {sum(p.stat().st_size for p in files):,} B on disk; "
      f"bodies {sum(v['bytes'] for v in manifest.values()):,} B as received")  # fmt: skip
for page in sorted((ROOT / "outputs/11_genecard").glob("2026-10-15_*_overview_v1.html")):
    html = page.read_text()
    size = len(html.encode())
    shown = unpacked(html)
    gone = re.findall(r"[^>]{0,40} unavailable for [^<]{0,160}", shown)
    print(
        f"{page.name}: {size:,} B ({size / MAX_BYTES:.1%} of {MAX_BYTES:,}); packed parts "
        f"{packed_bytes(html)}; unpacked {len(shown.encode()):,} B; unavailable {len(gone)}"
    )
    for line in gone:
        print("   ", line.strip())
    print("    nan as text:", len(re.findall(r">\s*nan\s*<", shown)))
