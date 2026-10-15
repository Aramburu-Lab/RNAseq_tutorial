"""Packed groups and entry pages in Chromium (C3.1, M9.8, GS.2).

A packed group's figures are drawn when the group comes near the view, wired as the page was
at load: in the dark theme, with their glossary terms; an entry page forwards a link's gene
list to the module's page. Skipped where no Chromium is found; set ``CTX_CHROMIUM``.
"""

from dataclasses import replace
from datetime import date

import pytest
from synth import demo_strip
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
from ctx.report.entry import write_entry
from ctx.viz.save import figure_block, write_figure

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
DAY = date(2026, 10, 14)
FRAGMENT = "#v=1&src=tcga_brca:tumour_normal&n=2&g=-75624,106236&end"
# Scroll to the last packed group, wait, and read what each group drew; then turn the theme
# to dark and back, so a part drawn while dark is checked separately.
READ = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(800);
  const before = Array.from(document.querySelectorAll('.packed-body'))
    .map(h => h.getAttribute('data-state'));
  document.querySelector('.packed-body:last-of-type').scrollIntoView();
  await wait(2500);
  const groups = Array.from(document.querySelectorAll('.group')).map(g => ({
    state: (g.querySelector('.packed-body') || {getAttribute: () => 'whole'})
      .getAttribute('data-state'),
    plots: g.querySelectorAll('.js-plotly-plot').length,
    drawn: Array.from(g.querySelectorAll('.js-plotly-plot')).filter(d => d._fullLayout).length,
    terms: g.querySelectorAll('button.info[data-term]').length,
    font: Array.from(g.querySelectorAll('.js-plotly-plot')).map(d => d._fullLayout.font.color)}));
  const term = document.querySelector('.packed-body button.info[data-term]');
  if (term) term.click();
  const tip = document.getElementById('ctx-tip');
  return {before, groups, tip: !!tip && !tip.hidden && term.parentNode.contains(tip)};
})()"""


def _page(out, packed: bool):
    """A module page of three groups, each one strip with its own glossary term in its caption:
    only a term's first mention in a section carries the icon (D-048)."""
    blocks = []
    for i in range(3):
        term = ("read", "fragment", "adapter")[i]
        data = replace(demo_strip(), caption=f"Each [[{term}]] counted; no test.")
        figure_stem = stem("genecard", f"demo_{i}", 1, day=DAY)
        text = FigureText(f"Strip {i}", reading=("Each point is one made-up sample.",))
        write_figure(data, out, figure_stem, title=f"Strip {i}", data_releases={}, text=text)
        blocks.append(figure_block(data, figure_stem, text))
    groups = tuple(Group(f"Card {i}", (b,), packed=packed and i > 0) for i, b in enumerate(blocks))
    built = ModulePage(
        title="M9 Gene cards: Demo", purpose="Purpose.", input="Input.", method="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),), results=groups,
        interpretation="Interpretation.", questions=(SelfCheck("A question.", "An answer."),),
        pitfalls=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    return write(built, out / f"{stem('genecard', 'demo_overview', 1, day=DAY)}.html")


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    out = tmp_path_factory.mktemp("packed") / "outputs" / "11_genecard"
    out.mkdir(parents=True)
    page = _page(out, packed=True)
    entry = write_entry(out, "genecard")
    dark = "localStorage.setItem('ctx-theme', 'black');"
    tasks = [
        {"name": "white", "url": page.as_uri(), "width": 1280, "height": 900, "act": READ},
        {"name": "phone", "url": page.as_uri(), "width": 375, "height": 800, "act": READ},
        # The theme is read before first paint: set it, then load the page again.
        {"name": "dark", "url": page.as_uri(), "width": 1280, "height": 900,
         "act": f"{dark} location.reload(); null", "resize": [{"width": 1280, "act": READ}]},
        {"name": "entry", "url": entry.as_uri() + FRAGMENT, "width": 1280, "height": 900,
         "act": "new Promise(r => setTimeout(r, 1500)).then(() => "
                "({href: location.href, length: history.length}))"},
    ]  # fmt: skip
    return page, _browse(tasks, tmp_path_factory.mktemp("downloads"))


@pytest.mark.parametrize("name", ["white", "phone"])
def test_packed_groups_draw_when_reached_and_wire_their_terms(ran, name):
    """Before the scroll only the first card is drawn; after it every card is, with its
    glossary term's tooltip opening.

    Mutations caught: every packed part drawn at load (the IntersectionObserver gone and
    unpack called on each); the plot scripts not run anew (no plot drawn); the tooltip
    binding skipped for a packed part."""
    _, results = ran
    seen = results[name]
    assert seen["errors"] == []
    acted = seen["acted"]
    assert acted["before"][-1] is None
    assert [g["state"] for g in acted["groups"]] == ["whole", "drawn", "drawn"]
    assert all(g["plots"] == g["drawn"] == 1 and g["terms"] == 1 for g in acted["groups"]), acted
    assert acted["tip"]


def test_a_part_drawn_in_the_dark_theme_takes_it(ran):
    """Mutation: the theme relayout left to the load, so a card drawn later is light."""
    _, results = ran
    (acted,) = results["dark"]["resized"]
    fonts = {f for g in acted["groups"] for f in g["font"]}
    assert len(fonts) == 1, fonts
    white = {f for g in results["white"]["acted"]["groups"] for f in g["font"]}
    assert fonts != white


def test_the_entry_page_forwards_the_list(ran):
    """Mutations caught: the fragment dropped on the way; a push in place of a replace, which
    leaves the entry page in Back's way."""
    page, results = ran
    acted = results["entry"]["acted"]
    assert acted["href"] == page.as_uri() + FRAGMENT
    assert acted["length"] == 2  # noqa: PLR2004 - about:blank, then this page in its place
