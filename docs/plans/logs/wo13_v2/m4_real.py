# ruff: noqa: E501 - page scripts and reports are quoted whole
"""V2 M4 F3 on a real page in Chromium: the order view, and each grouped view's files.

At 1,280 and 375 px in both themes; each view's PDF and SVG with vector cells and no raster.

python m4_real.py <root> <dataset> <downloads dir>
"""

import re
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402
from test_quant_views_browser import CELLS, ORDER, _rgb, _runs  # noqa: E402
from views_browser import walk  # noqa: E402

root, dataset, downloads = Path(sys.argv[1]).absolute(), sys.argv[2], Path(sys.argv[3]).absolute()
downloads.mkdir(parents=True, exist_ok=True)
(page,) = sorted((root / "outputs" / "04_quant").glob(f"*_quant_{dataset}_overview_v1.html"))[-1:]
SIZES = (("1280_light", 1280, "light"), ("1280_dark", 1280, "dark"),
         ("375_light", 375, "light"), ("375_dark", 375, "dark"))  # fmt: skip
print(f"page {page.name} {page.stat().st_size:,} B")
seen = walk(page, downloads / "walk", SIZES)
figure = None
for name, _, _ in SIZES:
    acted = seen[name]["acted"]
    (view,) = [v for v in acted["views"] if v["label"] == ORDER]
    figure = view["figure"]
    labels = [o["caption"].split(", grouped by ")[1].split(",")[0] if "grouped by" in o["caption"]
              else "clustering" for o in view["options"]]  # fmt: skip
    print(f"{name}: theme {acted['theme']}, options {labels}, moved "
          f"{[o['moved'] for o in view['options']]}, back unmoved {not view['back']['moved']}, "
          f"widest plot {max(o['wide'] for o in view['options']):.0f} px, box {acted['box']} "
          f"scroll {acted['scroll']}, errors {seen[name]['errors']}")  # fmt: skip
n = len(view["options"])
for option in range(1, n):
    for name, width, scheme in (("1280_light", 1280, "light"), ("375_dark", 375, "dark")):
        act = CELLS.replace("FIGURE", f"#{figure}").replace("OPTION", str(option))
        (result,) = _browse([{"name": "cells", "url": page.as_uri(), "width": width,
                              "scheme": scheme, "figure": f"#{figure}", "act": act}],
                            downloads / f"export_{option}_{name}", timeout=1200).values()  # fmt: skip
        acted = result["acted"]
        out = downloads / f"export_{option}_{name}"
        svg = (out / result["files"]["svg"]["guid"]).read_text()
        pdf = (out / result["files"]["pdf"]["guid"]).read_bytes()
        fills = re.findall(r'<path class="js-fill" d="([^"]+)" style="[^"]*?fill: ([^;]+);', svg)
        drawn: dict[str, int] = {}
        for path, colour in fills:
            drawn[colour] = drawn.get(colour, 0) + path.count("Z")
        expected: dict[str, int] = {}
        for z, colours in zip(acted["z"], acted["cells"], strict=True):
            for v, k in _runs(z).items():
                expected[_rgb(colours[v])] = expected.get(_rgb(colours[v]), 0) + k
        print(f"export {acted['label']} {name}: <image> {svg.count('<image')}, rectangles "
              f"{sum(drawn.values())} of {sum(expected.values())} runs, by colour equal "
              f"{drawn == expected}, svg {len(svg):,} B, pdf {len(pdf):,} B, pdf images "
              f"{pdf.count(b'/Subtype /Image')}, errors {result['errors']}")  # fmt: skip
