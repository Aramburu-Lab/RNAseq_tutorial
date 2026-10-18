"""R3: build the made-up study root and print its M7 page's paragraphs and tables."""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from proteo_study_root import run  # noqa: E402

root = run(Path(sys.argv[1]), sys.argv[2])
(page,) = (root / "outputs" / "07_proteo").glob("*_overview_v1.html")
html = page.read_text()
for line in re.findall(r"<p[^>]*>([^<]{40,400})</p>", html):
    print("-", line)
print(sorted(p.name for p in page.parent.glob("*.csv")))
