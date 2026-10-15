"""The page templates under the file-length ratchet, and the page script's one definition of each
function (E.4: ``_page.js.j2`` reached 1,100 lines and was split into partials by concern)."""

import re
from collections import Counter
from pathlib import Path

from ctx.report.builder import _ENV

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src" / "ctx" / "report" / "templates"


def _fail_limit() -> int:
    """The ratchet's failure line, read as the prek hook reads it."""
    text = (ROOT / ".labtemplate.yml").read_text()
    return int(re.search(r"^\s*fail:\s*(\d+)", text, re.M).group(1))


def test_every_template_is_within_the_file_length_ratchet():
    """The prek hook checks scripts/ and bin/ only; the templates are held here. Mutation: the
    page script joined again from its partials (1,100 lines)."""
    limit = _fail_limit()
    lengths = {p.name: len(p.read_text().splitlines()) for p in TEMPLATES.iterdir()
               if p.suffix == ".j2"}  # fmt: skip
    assert lengths["_page.js.j2"] > 0
    assert {n: k for n, k in lengths.items() if k > limit} == {}


def test_the_page_script_defines_each_of_its_functions_once():
    """Each function of the page script's one closure is defined once, whichever partial holds
    it, and every partial is included. Mutations: a partial included twice; one left out."""
    rendered = _ENV.get_template("_page.js.j2").render(live=True)
    # The closure follows the shared objects (ctxGlossary, ctxPick, ...), which it includes first.
    script = rendered[rendered.index("/* One script for every page kind.") :]
    names = Counter(re.findall(r"^  function (\w+)\(", script, re.M))
    assert {n: k for n, k in names.items() if k > 1} == {}
    expected = {"guard", "plots", "template", "retheme", "wire", "at", "fitting", "exporting",
                "ctxRows", "brushing", "pairing", "picking", "oraLive", "viewing", "filtering",
                "liveCalls"}  # fmt: skip
    assert set(names) == expected
    partials = sorted(p.name for p in TEMPLATES.glob("_page_*.js.j2"))
    page = (TEMPLATES / "_page.js.j2").read_text()
    assert partials and all(page.count(f'{{% include "{p}" %}}') == 1 for p in partials)
