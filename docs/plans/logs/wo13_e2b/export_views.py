"""Press the PDF and SVG buttons of views set away from the build, on a real page.

Prints what each file says: its size, its caption, its PDF's page size in pt.
python export_views.py <page.html> <downloads dir> <#figure-N>:<option> ...
"""

import re
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from pypdf import PdfReader  # noqa: E402
from views_browser import exported  # noqa: E402

page, folder = Path(sys.argv[1]).absolute(), Path(sys.argv[2]).absolute()
for item in sys.argv[3:]:
    figure, option = item.rsplit(":", 1)
    out = folder / f"{figure.strip('#')}_{option}"
    out.mkdir(parents=True, exist_ok=True)
    result = exported(page, out, figure, int(option))
    svg = (out / result["files"]["svg"]["guid"]).read_text()
    caption = " ".join(re.findall(r'<g class="caption">(.*?)</g>', svg, re.S))
    words = " ".join(re.findall(r">([^<>]+)</text>", caption))
    pdf = PdfReader(out / result["files"]["pdf"]["guid"])
    box = pdf.pages[0].mediabox
    size = f"{float(box.width):.1f} x {float(box.height):.1f} pt"
    print(f"{figure} option {option}: svg {len(svg):,} B, pdf {len(pdf.pages)} page {size}, "
          f"errors {result['errors']}")  # fmt: skip
    print(f"  export caption: {result['acted']}")
    print(f"  svg caption: {words}")
