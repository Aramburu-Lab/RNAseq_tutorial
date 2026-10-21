"""V.4b: a plain view whose options differ in feature count, on the page and in its files.

Each option carries its own least width and export settings (``views.viewed``): the page kept
the first option's least width, and with 4 then 11 features the 11 turned names overlapped
10 times at 1,280 px; and every option exported at the first option's size and labels. A box
view of 4, 11 and 20 features at 1,280 and 375 px in both themes. Skipped where no Chromium is
found; set ``CTX_CHROMIUM``.
"""

from datetime import date

import polars as pl
import pytest
from test_box_stack import boxes
from test_box_stack_browser import FILE, _apart
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
from ctx.report.export import file_caption, plain
from ctx.report.views import viewed
from ctx.viz.box import box_data
from ctx.viz.save import figure_block, write_figure
from ctx.viz.style import THEMES

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
COUNTS = (4, 11, 20)
# A last option of the build's 4 features from six samples a group: its own caption, "n = 6
# control, 6 case", at the build's size.
SAME = "same"
WIDTHS = (1280, 375)
SCHEMES = ("light", "dark")
# Sets option AT, then every drawn text of the plot, and what the fit's ctxBuilt sets back.
SHOWN = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.view');
  form.elements.view.value = 'AT';
  form.elements.view.dispatchEvent(new Event('change', {bubbles: true}));
  await wait(600);
  const box = document.querySelector('.fig-plot'), div = box.querySelector('.js-plotly-plot');
  const rect = e => {
    const b = e.getBoundingClientRect(); return [b.left, b.right, b.top, b.bottom];
  };
  const shown = e => e.getBoundingClientRect().width > 0;
  const built = div.ctxBuilt();
  return {ticks: Array.from(box.querySelectorAll('.xtick text')).map(rect),
          notes: Array.from(box.querySelectorAll('.annotation')).filter(shown).map(rect),
          legend: Array.from(box.querySelectorAll('.legend .legendtext')).map(rect),
          built: built && (built.relayout['xaxis.ticktext'] || []).length};
})()"""


def _six(data):
    """The same boxes from the first six samples of each group."""
    frame = data.frame.filter(
        ~pl.col("sample").str.ends_with("6") & ~pl.col("sample").str.ends_with("7")
    )
    return box_data(frame, columns=data.columns, feature_label="Gene", value_label="Expression",
                    groups=("control", "case"))  # fmt: skip


def _page(out):
    """A module page holding one plain box view of 4, 11 and 20 features, the 4 as built, and
    of 4 from six samples a group."""
    day = date(2026, 10, 21)
    contracts = {**{n: boxes(n) for n in COUNTS}, SAME: _six(boxes(4))}
    figure_stem = stem("de", "boxview_demo", 1, day=day)
    text = FigureText("Boxes", reading=("Each box is one gene in one group.",))
    write_figure(contracts[4], out, figure_stem, title="Box", data_releases={}, text=text)
    block = viewed(figure_block(contracts[4], figure_stem, text),
                   [(str(n), contracts[n], text) for n in [*COUNTS, SAME]], "Genes")  # fmt: skip
    page = ModulePage(
        title="M5 Demo: box view", purpose="Purpose.", input="Input.", method="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Figures", (block,)),),
        interpretation="Interpretation.", questions=(SelfCheck("A question.", "An answer."),),
        pitfalls=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    return write(page, out / f"{stem('de', 'boxview_overview', 1, day=day)}.html"), contracts


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    """Each option at each width and theme, each other option's export, its SVG measured. The
    first option's PDF and SVG are the figure files themselves (V.4's tests)."""
    out, downloads = tmp_path_factory.mktemp("view"), tmp_path_factory.mktemp("downloads")
    page, contracts = _page(out)
    tasks = [{"name": f"{n}_{w}_{s}", "url": page.as_uri(), "width": w, "height": 900,
              "scheme": s, "act": SHOWN.replace("AT", str(i)),
              **({"figure": "#figure-1"} if i else {})}
             for i, n in enumerate([*COUNTS, SAME]) for w in WIDTHS for s in SCHEMES]  # fmt: skip
    seen = _browse(tasks, downloads)
    files = []
    for name, result in seen.items():
        if not result["files"]:
            continue
        copy = out / f"{name}.svg"
        copy.write_bytes((downloads / result["files"]["svg"]["guid"]).read_bytes())
        files.append({"name": name, "url": copy.as_uri(), "act": FILE})
    measured = {k: r["acted"] for k, r in _browse(files, downloads).items()}
    return seen, measured, contracts


@pytest.mark.parametrize("n", COUNTS)
@pytest.mark.parametrize("width", WIDTHS)
@pytest.mark.parametrize("scheme", SCHEMES)
def test_the_shown_option_draws_every_name_apart(ran, n, width, scheme):
    """Every name of the shown option drawn, no two shown texts over each other, no error; a
    phone layout's ``ctxBuilt`` sets back the shown option's n names. Mutation caught: the page
    left at the first option's least width, where 11 turned names overlapped 10 times."""
    seen = ran[0][f"{n}_{width}_{scheme}"]
    assert seen["errors"] == []
    drawn = seen["acted"]
    assert len(drawn["ticks"]) == n
    assert _apart(drawn["ticks"] + drawn["notes"] + drawn["legend"])
    assert drawn["built"] in (None, n)


@pytest.mark.parametrize("n", COUNTS[1:])
@pytest.mark.parametrize("width", WIDTHS)
@pytest.mark.parametrize("scheme", SCHEMES)
def test_each_option_exports_at_its_own_size(ran, n, width, scheme):
    """The exported SVG of each option: that option's files' size, its n names and n labels,
    no two texts over each other. Mutation caught: every option exported with the first
    option's settings, 4 features' size and labels."""
    _, measured, contracts = ran
    svg = measured[f"{n}_{width}_{scheme}"]
    drawn = contracts[n].static(THEMES["white"]).get_size_inches()
    assert svg["page"] == pytest.approx([x * 72 for x in drawn], abs=0.01)
    texts = svg["texts"]
    assert len([t for t in texts if t[0].startswith("p=")]) == n
    assert len([t for t in texts if t[0].startswith("GENE")]) == n
    assert _apart([t[1:] for t in texts])


@pytest.mark.parametrize("width", WIDTHS)
@pytest.mark.parametrize("scheme", SCHEMES)
def test_an_option_of_the_builds_size_exports_its_own_caption(ran, width, scheme):
    """An option of the build's 4 features from six samples a group: its SVG at the build's
    size, with that option's caption, not the build's. Mutation caught: an option of the
    build's size left on the build's export settings, which printed "n = 8 control"."""
    _, measured, contracts = ran
    svg = measured[f"{SAME}_{width}_{scheme}"]
    drawn = contracts[SAME].static(THEMES["white"]).get_size_inches()
    assert svg["page"] == pytest.approx([x * 72 for x in drawn], abs=0.01)
    words = " ".join(t[0] for t in svg["texts"])
    own, built = (plain(file_caption(contracts[k].caption)) for k in (SAME, 4))
    assert own != built and own in words and built not in words
