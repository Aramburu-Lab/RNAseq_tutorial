# ruff: noqa: D205, PLW1510 - a proof record, run once and kept as run
"""One Chromium, the legend test's four tasks on the before and after pages, alternating, three
rounds; each task's seconds from its navigate to its close, from the driver's trace.
"""

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "tests")
from test_export_browser import CHROMIUM  # noqa: E402

READ = """new Promise(r => setTimeout(r, 2500)).then(() => ({
  n: document.querySelectorAll('.js-plotly-plot').length, t: performance.now()}))"""
pages = {"before": Path(sys.argv[1]), "after": Path(sys.argv[2])}
tasks = []
for k in range(3):
    for which in ("before", "after") if k % 2 == 0 else ("after", "before"):
        for w in (375, 1280):
            for s in ("light", "dark"):
                tasks.append(
                    {
                        "name": f"{which}:{w}:{s}:{k}",
                        "url": pages[which].as_uri(),
                        "width": w,
                        "height": 900,
                        "scheme": s,
                        "act": READ,
                    }
                )
plan = {
    "chromium": CHROMIUM,
    "downloads": str(pages["after"].parent),
    "tasks": tasks,
    "trace": True,
}
done = subprocess.run(
    [shutil.which("node"), "tests/browser.mjs"],
    input=json.dumps(plan),
    capture_output=True,
    text=True,
    timeout=3000,
)
stamps = re.findall(r"trace ([\d.]+) (\S+): (\w+)", done.stderr)
start, spans = {}, {}
for t, name, step in stamps:
    if step == "navigate":
        start[name] = float(t)
    if step == "close" and name in start:
        spans[name] = float(t) - start[name]
for name, span in spans.items():
    print(f"{name} {span:.1f} s")
for which in ("before", "after"):
    v = [s for n, s in spans.items() if n.startswith(which)]
    print(which, "tasks", len(v), "total", round(sum(v), 1), "s")
print("exit", done.returncode, done.stderr[-300:] if done.returncode else "")
