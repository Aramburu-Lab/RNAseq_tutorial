"""wo13-z: the rebuilt Step-1 M2 page's Reproduce tabs in headless Chromium.

The ACT clicks every panel and keeps the last panel's tabs: each one's copy and download at
1,280 and 375 px, light and dark, against what panels() in tests/test_figure_panel.py reads.
"""

import base64
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tests"))
from test_export_browser import _browse  # noqa: E402
from test_figure_panel import panels  # noqa: E402
from test_reproduce_browser import ACT  # noqa: E402

page = next((ROOT / "outputs/scratch_z/z/tier1/outputs/02_qc").glob("*_overview_v1.html"))
served = page.read_text()
want = {tab: panels(served, tab)[-1] for tab in ("docker", "apptainer", "pixi")}
tasks = [{"name": f"{w} {s}", "url": page.as_uri(), "width": w, "height": 900, "scheme": s,
          "act": ACT} for w in (1280, 375) for s in ("light", "dark")]  # fmt: skip
for name, seen in _browse(tasks, Path(tempfile.mkdtemp())).items():
    acted, ok = seen["acted"], []
    first = {tab: got for tab, got in acted["tabs"].items()}
    for tab, got in first.items():
        script = base64.b64decode(got["href"].split(",", 1)[1]).decode()
        ok.append(got["copied"] == want[tab][0] and script == want[tab][1])
    print(name, acted["theme"], "tabs", list(first), "copy and download equal", ok,
          "shared hidden", acted["shared"], "width", acted["wide"], "of", acted["window"],
          "console errors", len(seen["errors"]), flush=True)  # fmt: skip
