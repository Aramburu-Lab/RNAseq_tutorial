"""V.4: a box view of 11, 15 and 20 features on the page and in its PDF and SVG, in Chromium.

On a module page from 11 features the box draws a phone's layout, its names turned with their
p and r beneath and the labels over the boxes hidden, and scrolls sideways at 1,280 and 375 px
(2026-10-21). The view's export was that layout redrawn at the files' 482 pt, in the one-row
size; it is now the files' layout, every label in the files' rows (``Stack``), at the files'
size, with its legend clear of the highest label. Skipped where no Chromium is found; set
``CTX_CHROMIUM``.
"""

from datetime import date

import pytest
from test_box_stack import boxes
from test_export_browser import CHROMIUM, _browse

from ctx.io.naming import stem
from ctx.report.builder import (
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    write,
)
from ctx.report.lookups import lookup
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import STACK_GAP_PT, THEMES

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
COUNTS = (11, 15, 20)
# The p of the view the page opens on; the view shown and exported prints 3.7e-05, wider.
BUILT_P = 1.2e-03
# Sets the view, then every drawn text of the plot: tick names, shown labels, legend entries.
SHOWN = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.view');
  form.elements.view.value = 'AT';
  form.elements.view.dispatchEvent(new Event('change', {bubbles: true}));
  await wait(400);
  const box = document.querySelector('.fig-plot');
  const rect = e => {
    const b = e.getBoundingClientRect(); return [b.left, b.top, b.right, b.bottom];
  };
  const shown = e => e.getBoundingClientRect().width > 0;
  return {ticks: Array.from(box.querySelectorAll('.xtick text')).map(rect),
          notes: Array.from(box.querySelectorAll('.annotation')).filter(shown).map(rect),
          legend: Array.from(box.querySelectorAll('.legend .legendtext')).map(rect)};
})()"""
# Every text of an exported SVG, its legend's box and its labels' boxes, in the file's own pt.
# The boxes are plotly's, written at the export: a text's own box follows the face the viewer
# finds for the SVG's font names, 0.25 pt taller here than the page's.
FILE = """(function () {
  var svg = document.documentElement, view = svg.viewBox.baseVal;
  var origin = svg.getBoundingClientRect(), k = view.width / origin.width;
  function box(e) {
    var r = e.getBoundingClientRect();
    return [(r.left - origin.left) * k, (r.right - origin.left) * k,
            (r.top - origin.top) * k, (r.bottom - origin.top) * k];
  }
  return {page: [view.width, view.height], legend: box(svg.querySelector('.legend rect.bg')),
          notes: Array.from(svg.querySelectorAll('.annotation rect.bg')).map(box),
          texts: Array.from(svg.querySelectorAll('text')).filter(function (t) {
            return t.textContent.trim(); }).map(function (t) {
            return [t.textContent].concat(box(t)); })};
})()"""


def _apart(boxes: list[list[float]]) -> bool:
    """No two [left, right, top, bottom] boxes overlap."""
    return all(a[1] <= b[0] or b[1] <= a[0] or a[3] <= b[2] or b[3] <= a[2]
               for i, a in enumerate(boxes) for b in boxes[i + 1 :])  # fmt: skip


def _page(out, n: int):
    """A module page holding one box view: the build at ``BUILT_P``, the other at 3.7e-05."""
    day = date(2026, 10, 21)
    built, other = boxes(n, p=BUILT_P), boxes(n)
    figure_stem = stem("de", f"box{n}_demo", 1, day=day)
    text = FigureText("Boxes", reading=("Each box is one gene in one group.",))
    write_figure(built, out, figure_stem, title="Box", data_releases={}, text=text)
    block = lookup(figure_block(built, figure_stem, text),
                   [("built", built, text), ("other", other, text)], "p")  # fmt: skip
    page = ModulePage(
        title="M5 Demo: box", purpose="Purpose.", input="Input.", method="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Figures", (block,)),),
        interpretation="Interpretation.", questions=(SelfCheck("A question.", "An answer."),),
        pitfalls=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    return write(page, out / f"{stem('de', f'box{n}_overview', 1, day=day)}.html"), other


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    """Each page at 1,280 and 375 px on the other view, its export at both, each file measured."""
    out, downloads = tmp_path_factory.mktemp("box"), tmp_path_factory.mktemp("downloads")
    tasks, contracts = [], {}
    for n in COUNTS:
        page, contracts[n] = _page(out, n)
        tasks += [{"name": f"{n}_{w}", "url": page.as_uri(), "width": w, "height": 900,
                   "act": SHOWN.replace("AT", "1"), "figure": "#figure-1"}
                  for w in (1280, 375)]  # fmt: skip
    seen = _browse(tasks, downloads)
    files = []
    for name, result in seen.items():
        copy = out / f"{name}.svg"
        copy.write_bytes((downloads / result["files"]["svg"]["guid"]).read_bytes())
        files.append({"name": name, "url": copy.as_uri(), "act": FILE})
    measured = {k: r["acted"] for k, r in _browse(files, downloads).items()}
    return seen, measured, contracts


@pytest.mark.parametrize("n", COUNTS)
@pytest.mark.parametrize("width", [1280, 375])
def test_the_page_draws_every_name_and_shown_label_apart(ran, n, width):
    """The page as it stood: every name drawn, no two shown texts over each other, no error. At
    both widths the box draws its phone layout, every label under its name."""
    seen = ran[0][f"{n}_{width}"]
    assert seen["errors"] == []
    drawn = seen["acted"]
    assert len(drawn["ticks"]) == n and drawn["notes"] == []
    assert _apart([[a, c, b, d] for a, b, c, d in drawn["ticks"] + drawn["legend"]])


@pytest.mark.parametrize("n", COUNTS)
@pytest.mark.parametrize("width", [1280, 375])
def test_the_export_draws_the_files_rows_at_the_files_size(ran, n, width):
    """The exported SVG: the files' size, every label and name drawn, no two texts over each
    other or off the page, the legend at least ``STACK_GAP_PT`` over the highest label and less
    than a pt more. Mutations caught: ``inches()`` left at the one-row size (265 pt, where
    the files draw 284 at 20); the export left without the files' layout (no label drawn);
    the rows at the files' 20.6 pt pitch, where plotly's two-line labels stand 22 pt tall;
    the legend left where plotly puts it (its box 4.86 to 33.86 pt, the labels from 30 at
    11)."""
    _, measured, contracts = ran
    svg = measured[f"{n}_{width}"]
    drawn = contracts[n].static(THEMES["white"]).get_size_inches()
    assert svg["page"] == pytest.approx([x * 72 for x in drawn], abs=0.01)
    texts = svg["texts"]
    notes = [t for t in texts if t[0].startswith("p=")]
    assert len(notes) == n and len([t for t in texts if t[0].startswith("GENE")]) == n
    assert _apart([t[1:] for t in texts])
    width_pt, height_pt = svg["page"]
    assert all(t[1] >= 0 and t[2] <= width_pt and t[3] >= 0 and t[4] <= height_pt for t in texts)
    gap = min(b[2] for b in svg["notes"]) - svg["legend"][3]
    # Lifted to the drawing's top at each count: plotly's rows are taller than the files'.
    assert svg["legend"][2] == pytest.approx(0, abs=0.01)
    assert STACK_GAP_PT - 0.01 <= gap < STACK_GAP_PT + 1
