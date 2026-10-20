"""wo13-y item 3: M9's figure files, specs and each Reproduce tab's ``ctx figure`` lines.

Builds M9's made-up root as test_figures_m9 does (or reads one built before), to find where the
24 extra ``ctx figure`` lines come from. Usage: python m9_panels.py <root>.
"""

import html
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tests"))
from test_figure_panel import panels  # noqa: E402
from test_figures import specs  # noqa: E402
from test_figures_m9 import FOLDER, _built  # noqa: E402

from ctx.report.builder import unpacked  # noqa: E402

root = Path(sys.argv[1])
spec = _built(root) if not (root / "outputs").exists() else next(root.glob("*.toml"))
folder = root / "outputs" / FOLDER[0]
(path,) = folder.glob("*_overview_v1.html")
shown = html.unescape(unpacked(path.read_text()))
print("ctx figure lines:", len(re.findall(r"ctx figure \S+ --params", shown)))
for tab in ("docker", "apptainer", "pixi"):
    lines = Counter(re.search(r"ctx figure (\S+) --params", p).group(1)
                    for p, _ in panels(unpacked(path.read_text()), tab))  # fmt: skip
    print(tab, sum(lines.values()), dict(lines))
print("white PDFs:", len(list(folder.glob("*_white.pdf"))))
found = specs(root, FOLDER, spec)
print("specs:", len(found), sorted(s.stem[11:] for s in found))
