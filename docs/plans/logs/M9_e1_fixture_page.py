"""Package M9, E.1: a page with the evidence view for the browser check, on the test fixture.

The lung datasets have one contrast each, so their pages carry no evidence view. This builds
the M9 tests' made-up root (NPTX2's real HPA and Ensembl responses, made-up GTEx and M5)
with a second contrast, under the scratch folder given, and prints the page's path.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tests"))

from genecard_root import CALLED, SLUG, TWINS, Stub, make, responses  # noqa: E402
from test_genecard import _cards, _register, _settings  # noqa: E402

from ctx.de.run import newest  # noqa: E402
from ctx.genecard import run as module  # noqa: E402

folder = Path(sys.argv[1]) / "root"
spec = make(folder)
_register(folder, spec)
results = newest(folder / "outputs" / "05_de", SLUG, "a_case_results")
for suffix in (".parquet", ".csv"):
    source = results.with_suffix(suffix)
    source.with_name(source.name.replace("a_case", "b_case")).write_bytes(source.read_bytes())
spec.write_text(
    spec.read_text() + '\n[[contrasts]]\nname = "b_case"\nvariable = "condition"'
    '\nreference = "control"\ntested = "case"\nn = [3, 3]\nlabel = "B case"\n'
)
named = ("NPTX2", *(t.split(".")[0] for t in TWINS), *[v[0] for v in CALLED.values()])
written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG,
                     module.Options(spec, genes=named), Stub(responses()))  # fmt: skip
print(written[0])
