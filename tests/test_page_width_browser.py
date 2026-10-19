"""No module page lays out wider than a 375 px phone, in either theme (L, 2026-10-19).

Every fixture page the session builds is opened at 375 px, light and dark; the page's scroll
width must equal its box. A failure lists what reaches past the box outside any sideways
scroller: an element, or a run of text, which no element's box shows (an unbroken source path
in a reading guide widened Step-1's M2 page to 712 px). Skipped where no Chromium is found;
set ``CTX_CHROMIUM``.
"""

import dataclasses
import tempfile
from pathlib import Path

import pytest
from test_export_browser import CHROMIUM, _browse, _module

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
PHONE = 375
# The page box and its scroll width once the page script settles, and what runs past the box
# outside a sideways scroller (a table wrap, a figure's plot, the contents bar).
WIDTH = """(async () => {
  await new Promise(r => setTimeout(r, 2500));
  const w = document.documentElement.clientWidth, out = [];
  const scrolls = (e) => { for (let p = e; p && p !== document.body; p = p.parentElement) {
    if (["auto", "scroll", "hidden"].includes(getComputedStyle(p).overflowX)) return true; }
    return false; };
  for (const e of document.querySelectorAll("body *")) {
    const r = e.getBoundingClientRect();
    if (r.width > 0 && r.right > w + 0.5 && !scrolls(e.parentElement)) {
      out.push([e.tagName, Math.round(r.right), (e.textContent || "").trim().slice(0, 60)]);
    }
  }
  const walk = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  for (let n = walk.nextNode(); n; n = walk.nextNode()) {
    if (!n.textContent.trim() || scrolls(n.parentElement)) continue;
    const range = document.createRange();
    range.selectNodeContents(n);
    const r = range.getBoundingClientRect();
    if (r.width > 0 && r.right > w + 0.5) {
      out.push(["#text", Math.round(r.right), n.textContent.trim().slice(0, 60)]);
    }
  }
  return {box: w, scroll: document.documentElement.scrollWidth, wide: out.slice(0, 6)};
})()"""
# The session's module fixtures, each (root, page), and the folder its pages are written to.
FIXTURES = {"setup_run": "00_setup", "gather_run": "01_gather", "qc_run": "02_qc",
            "process_run": "03_process", "quant_run": "04_quant", "de_run": "05_de",
            "function_run": "06_function", "proteo_run": "07_proteo", "ncrna_run": "08_ncrna",
            "genecard_run": "11_genecard", "clinical_run": "12_clinical",
            "report_run": "09_report", "guide_run": "09_report"}  # fmt: skip
PAGES = ("_overview_v1", "_hub_v1", "_environment_v1", "_report_guide_v1", "_report_glossary_v1")


def _wide(pages: dict[str, Path]) -> dict[str, dict]:
    """Each page at 375 px in both themes: the pages whose scroll width passes their box."""
    tasks = [{"name": f"{name} {scheme}", "url": page.as_uri(), "width": PHONE, "height": 812,
              "scheme": scheme, "act": WIDTH}
             for name, page in pages.items() for scheme in ("light", "dark")]  # fmt: skip
    with tempfile.TemporaryDirectory() as folder:
        seen = _browse(tasks, Path(folder))
    assert len(seen) == len(tasks)
    out = {}
    for name, result in seen.items():
        acted = result["acted"]
        assert acted["box"] == PHONE, (name, acted)
        if acted["scroll"] > acted["box"]:
            out[name] = acted
    return out


def test_no_fixture_page_is_wider_than_a_phone(request):
    """Every module page the session builds, the hub, its glossary and the guide.

    Catches: the body's break-word rule removed (a long word then runs past its box)."""
    pages = {}
    for fixture, folder in FIXTURES.items():
        root = request.getfixturevalue(fixture)[0]
        for page in (root / "outputs" / folder).glob("*.html"):
            if page.stem.endswith(PAGES):
                pages[f"{fixture} {page.stem[11:]}"] = page
    # Every fixture gave a page; the report folders hold the hub, glossary and guide.
    assert {name.split()[0] for name in pages} == set(FIXTURES), sorted(pages)
    assert _wide(pages) == {}


def test_a_long_word_in_page_text_wraps_inside_its_box(tmp_path):
    """The rule itself, on a page whose reading guide and note hold one word longer than a
    phone: a source path as M2's band shading names, a gene-set name as M7's class note does.

    Catches: the body's break-word rule removed (the fixture pages may hold no word that long)."""
    import plotly.graph_objects as go  # noqa: PLC0415 - heavy; this test alone draws one

    from ctx.report.builder import FigureBlock, FigureText, write  # noqa: PLC0415

    path = "Graphs/QualityBoxPlot.java(github.com/s-andrews/FastQC/blob/v0.12.1/uk/ac/babraham)"
    sets = "GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT_AND_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX"
    reading = (f"Shading, display only: {path}",)
    block = FigureBlock(go.Figure(go.Scatter(x=[1, 2], y=[1, 2])), "A caption.", (), "",
                        FigureText("A figure", reading=reading))  # fmt: skip
    page = _module([block])
    page = dataclasses.replace(page, results=(f"Gene classes from {sets}.", *page.results))
    out = write(page, tmp_path / "2026-10-19_long_words_v1.html")
    assert _wide({"long": out}) == {}
