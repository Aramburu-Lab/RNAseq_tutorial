"""Package W: one line per task of an m9_browser log; run python m9_summary.py <log> ...."""

import json
import sys
from pathlib import Path

for name in sys.argv[1:]:
    text = Path(name).read_text()
    d = json.loads(text[text.index("{") : text.rindex("}") + 1])
    print(f"## {Path(name).name}")
    for k, v in d.items():
        if k.startswith("move"):
            print(k, v["theme"], v["width"], "fits", v["fits"], "F3", len(v["f3"]), "bad",
                  sum(v["f3"]), "F4", len(v["f4"]), "bad", sum(v["f4"]), "tau defined",
                  v["tau_defined"], v["problems"][:2])  # fmt: skip
    print("errors", d["errors"])
    print("link", {k: d["link"][k] for k in ("sent", "tested", "carded", "page", "summary",
                                            "pages", "fits", "errors")})  # fmt: skip
    print("first rows", d["link"]["first_rows"][:2])
    print("beside", list(d["beside"].values()))
