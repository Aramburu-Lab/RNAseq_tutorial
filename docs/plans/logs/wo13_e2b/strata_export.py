"""Press F3's PDF and SVG with the cut-offs moved; count the points and print the caption.

python strata_export.py <page.html> <downloads dir>
"""

import re
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402

page, folder = Path(sys.argv[1]).absolute(), Path(sys.argv[2]).absolute()
folder.mkdir(parents=True, exist_ok=True)
act = """(async () => {
  await new Promise(r => setTimeout(r, 1500));
  const form = document.querySelector('form.strata');
  form.elements.basis.value = 'pvalue'; form.elements.fc.value = '0.5';
  form.elements.fc.dispatchEvent(new Event('change', {bubbles: true}));
  await new Promise(r => setTimeout(r, 200));
  return form.closest('.fig').id;
})()"""
probe = _browse([{"name": "p", "url": page.as_uri(), "act": act}], folder)["p"]
(result,) = _browse([{"name": "e", "url": page.as_uri(), "act": act,
                      "figure": f"#{probe['acted']}"}], folder).values()  # fmt: skip
svg = (folder / result["files"]["svg"]["guid"]).read_text()
print(
    f"svg {len(svg):,} B; point paths {svg.count('class="point"')}; "
    f"marker paths {len(re.findall(r'<path class=\"point', svg))}; errors {result['errors']}"
)
print(" ".join(re.findall(r">([^<>]+)</text>", svg))[-600:])
