"""wo13-s2: each M10 page's missing-values rows as text. Usage: missing_table.py <page>..."""

import re
import sys
from html import unescape
from pathlib import Path

for name in sys.argv[1:]:
    page = Path(name).read_text()
    start = page.index("Missing values and their reasons")
    table = page[start : page.index("</table>", start)]
    print("==", "/".join(Path(name).parts[-5:-3]))
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", table, re.S):
        found = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S)
        print(" | ".join(unescape(re.sub(r"<[^>]+>", "", c)).strip() for c in found))
