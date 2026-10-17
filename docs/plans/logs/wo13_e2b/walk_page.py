"""Walk every view option of a real page at 1,280 and 375 px in both themes; print a summary.

python walk_page.py <page.html> <downloads dir> [label ...]; with labels, each named view's
option leads, captions and legend names are printed as well.
"""

import json
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from views_browser import walk  # noqa: E402

SIZES = (("1280_light", 1280, "light"), ("1280_dark", 1280, "dark"),
         ("375_light", 375, "light"), ("375_dark", 375, "dark"))  # fmt: skip
page, downloads = Path(sys.argv[1]).absolute(), Path(sys.argv[2]).absolute()
labels = sys.argv[3:]
downloads.mkdir(parents=True, exist_ok=True)
seen = walk(page, downloads, SIZES)
for name, result in seen.items():
    acted = result["acted"]
    print(f"== {name}: theme {acted['theme']}, page box {acted['box']} px, scroll "
          f"{acted['scroll']} px, errors {result['errors']}")  # fmt: skip
    for view in acted["views"]:
        same = all([(t["name"], t["y"]) for t in d["traces"]] ==
                   [(t["name"], t["y"]) for t in b]
                   for d, b in zip(view["options"], view["built"], strict=True))  # fmt: skip
        moved = [o["moved"] for o in view["options"]]
        widest = max(o["wide"] for o in view["options"])
        back = view["back"]["traces"] == view["options"][0]["traces"]
        print(f"{view['figure']} {view['label']!r}: {len(view['options'])} options, each drawn "
              f"as built {same}, build unmoved {not moved[0]}, others moved {all(moved[1:])}, "
              f"back to build {back}, widest plot {widest:.0f} px")  # fmt: skip
        if name == "1280_light" and view["label"] in labels:
            for i, o in enumerate(view["options"][:3] if view["grid"] else view["options"]):
                print(json.dumps({"option": i, "line": o["line"], "lead": o["lead"],
                                  "caption": o["caption"],
                                  "names": [t["name"] for t in o["traces"]],
                                  "notes": o["notes"]}, ensure_ascii=False))  # fmt: skip
