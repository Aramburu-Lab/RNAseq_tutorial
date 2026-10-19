# ruff: noqa: E501 - page scripts and reports are quoted whole
"""V2 M2's band map flag view on a real page in Chromium, at 1,280 and 375 px in both themes.

python m2_real.py <root> <dataset> <downloads dir>: the view's flagged cells against the run's
flags table, its other levels counted, the page box; its PDF and SVG with no raster image.
"""

import json
import re
import sys
from pathlib import Path

import polars as pl

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402
from test_qc_views import BANDS, SWAP  # noqa: E402
from views_browser import walk  # noqa: E402

from ctx.qc import flag_map  # noqa: E402
from ctx.report import text  # noqa: E402

root, dataset, downloads = Path(sys.argv[1]).absolute(), sys.argv[2], Path(sys.argv[3]).absolute()
downloads.mkdir(parents=True, exist_ok=True)
out = root / "outputs" / "02_qc"
(page,) = sorted(out.glob(f"*_qc_{dataset}_overview_v1.html"))[-1:]
flags = pl.read_parquet(sorted(out.glob(f"*_qc_{dataset}_flags_v1.parquet"))[-1])
SIZES = (("1280_light", 1280, "light"), ("1280_dark", 1280, "dark"),
         ("375_light", 375, "light"), ("375_dark", 375, "dark"))  # fmt: skip
print(f"page {page.name} {page.stat().st_size:,} B; flags table {flags.height} rows, "
      f"{flags.select('sample', 'metric').unique().height} sample-metric pairs")  # fmt: skip
words = {k: text.word("qc", f"figures.flag_map.levels.{k}") for k in flag_map.LEVELS}
seen = walk(page, downloads / "walk", SIZES)
for name, _, _ in SIZES:
    acted = seen[name]["acted"]
    (view,) = [v for v in acted["views"] if v["options"][0]["title"] == BANDS]
    shown = view["options"][1]
    print(f"{name}: theme {acted['theme']}, options {len(view['options'])}, title "
          f"{shown['title']!r}, moved {[o['moved'] for o in view['options']]}, back unmoved "
          f"{not view['back']['moved']}, widest plot {max(o['wide'] for o in view['options']):.0f}"
          f" px, box {acted['box']} scroll {acted['scroll']}, errors {seen[name]['errors']}")  # fmt: skip
for found in re.finditer(r'id="ctx-controls-\d+">(.*?)</script>', page.read_text(), re.S):
    cfg = json.loads(found.group(1))
    if cfg.get("kind") == "view" and cfg["options"][0].get("title") == BANDS:
        break
cells = next(t for t in cfg["options"][1]["figure"]["data"] if t["type"] == "heatmap")
counts = {k: sum(c.endswith(f"<br>{w}") for row in cells["text"] for c in row)
          for k, w in words.items()}  # fmt: skip
print(f"cells by level {counts}; flagged equals the table's pairs: "
      f"{counts['flagged'] == flags.select('sample', 'metric').unique().height}")  # fmt: skip
for name, width, scheme in (("1280_light", 1280, "light"), ("375_dark", 375, "dark")):
    task = {"name": "v", "url": page.as_uri(), "width": width, "height": 900, "scheme": scheme,
            "act": SWAP % json.dumps(BANDS), "figure": "[data-probe]"}  # fmt: skip
    got = _browse([task], downloads / f"export_{name}", timeout=1200)["v"]
    svg = (downloads / f"export_{name}" / got["files"]["svg"]["guid"]).read_text()
    pdf = (downloads / f"export_{name}" / got["files"]["pdf"]["guid"]).read_bytes()
    print(f"export {name}: <image> {svg.count('<image')}, filled paths "
          f"{svg.count('class=\"js-fill\"')}, svg {len(svg):,} B, pdf {len(pdf):,} B, pdf images "
          f"{pdf.count(b'/Subtype /Image')}, line {got['acted']['line'][:60]!r}, errors "
          f"{got['errors']}")  # fmt: skip
