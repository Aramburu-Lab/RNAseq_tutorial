"""Package W, E.1a: what the list panel costs on a built M9 page, and how many cards it holds.

For each page: its bytes, the panel's settings and their block (gzip, base64), the block's
JSON as plain text, the genes and contrasts it carries, and the cards drawn whole and packed.
Run: python m9_bytes.py <page.html> [...]
"""

import json
import re
import sys
from pathlib import Path

from ctx.genecard.listing import unblock
from ctx.report.builder import packed_bytes

PANEL = re.compile(r'<script type="application/json" id="ctx-panel-1">([^<]*)</script>')


def main() -> None:
    """One line of numbers per page."""
    for name in sys.argv[1:]:
        html = Path(name).read_text()
        found = PANEL.search(html)
        row = {"page": Path(name).name, "page_bytes": len(html.encode()),
               "packed_cards": len(packed_bytes(html))}  # fmt: skip
        if found:
            settings = json.loads(found[1])
            data = unblock(settings["block"])
            plain = json.dumps(data, separators=(",", ":"))
            row |= {"panel_json_bytes": len(found[1].encode()),
                    "block_bytes": len(settings["block"]), "block_plain_json_bytes": len(plain),
                    "genes": len(data["numbers"]), "contrasts": len(data["calls"]),
                    "cards_linked": len(settings["cards"])}  # fmt: skip
        print(json.dumps(row))


if __name__ == "__main__":
    main()
