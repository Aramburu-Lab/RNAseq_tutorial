"""wo13-z1: M10 on the made-up GDC root of tests/clinical_root.py, with the src PYTHONPATH names.

The page's path is printed. Run: python fixture_page.py <empty folder>.
"""

import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parents[4] / "tests"
sys.path.insert(0, str(TESTS))

from clinical_root import ctx, make  # noqa: E402 - tests/ on the path first

root, spec = make(Path(sys.argv[1]))
done = ctx("--root", str(root), "clinical", "--dataset", "fake_brca", "--dataset-file", str(spec))
(page,) = (root / "outputs" / "12_clinical").glob("*_clinical_fake_brca_overview_v1.html")
print(page)
