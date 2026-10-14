"""The ``view`` kind (D-051) on M5's p-value histogram and M6's GSEA bars, by contrast
(work-order items 6 and 7): in Chromium on the pipeline fixture's pages, and its contract.

Each test names the mutation it catches.
"""

import pytest
from test_export_browser import CHROMIUM, _browse

from ctx.report.builder import FigureText
from ctx.report.figures import prepare
from ctx.report.views import viewed

ACT = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(2000);
  return await Promise.all(Array.from(document.querySelectorAll('form.view')).map(async form => {
    const fig = form.closest('.fig'), d = fig.querySelector('.js-plotly-plot');
    const words = s => (fig.querySelector(s) || {}).textContent;
    const read = () => ({y: JSON.stringify(d.data.map(t => t.y)), caption: words('figcaption'),
      title: words('.fig-title'), line: words('.view-summary'),
      moved: !!fig.querySelector('a[data-export][role=button]'),
      options: form.elements.view.options.length});
    const out = {start: read()};
    form.elements.view.value = '1'; form.dispatchEvent(new Event('change'));
    await wait(300); out.other = read();
    form.elements.view.value = '0'; form.dispatchEvent(new Event('change'));
    await wait(300); out.back = read();
    return out;
  }));
})()"""


@pytest.fixture(scope="module")
def ran(pipeline, tmp_path_factory):
    root, _, _ = pipeline
    (m5,) = (root / "outputs" / "05_de").glob("*_de_desynth_overview_v1.html")
    (m6,) = (root / "outputs" / "06_function").glob("*_function_desynth_overview_v1.html")
    tasks = [{"name": n, "url": p.as_uri(), "width": 1280, "height": 900, "act": ACT}
             for n, p in (("m5", m5), ("m6", m6))]  # fmt: skip
    return _browse(tasks, tmp_path_factory.mktemp("downloads"))


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
@pytest.mark.parametrize("page", ["m5", "m6"])
def test_a_view_swaps_the_figure_and_its_words_and_comes_back(ran, page):
    """Mutations caught: the figure left as built; the caption not swapped; the build's words
    not restored; the export switch left off."""
    seen = ran[page]
    assert seen["errors"] == []
    (view,) = seen["acted"]
    start, other, back = view["start"], view["other"], view["back"]
    assert start["options"] > 1 and not start["moved"] and start["line"] == ""
    assert other["y"] != start["y"] and other["caption"] != start["caption"]
    assert other["moved"] and other["line"].startswith("This view shows ")
    assert back["caption"] == start["caption"] and back["title"] == start["title"]
    assert back["y"] == start["y"] and not back["moved"]


def test_a_view_needs_its_build_first_and_is_dropped_alone():
    """Mutations caught: one option accepted; another figure taken as the build; the select
    kept on a page drawn at twice the size its options were drawn at."""
    from dataclasses import replace  # noqa: PLC0415

    from synth import counts, demo_bar  # noqa: PLC0415

    from ctx.viz.save import figure_block  # noqa: PLC0415

    data = demo_bar(counts())
    other = replace(demo_bar(counts(seed=7)), caption="Another made-up count; no test.")
    words = FigureText("Bars", reading=("Each bar is made up.",))
    block = figure_block(data, "2026-10-14_de_view_demo_v1", words)
    with pytest.raises(ValueError, match="two or more"):
        viewed(block, [("a", data, words)], "Which")
    with pytest.raises(ValueError, match="as built"):
        viewed(block, [("b", other, words), ("a", data, words)], "Which")
    shown = viewed(block, [("a", data, words), ("b", other, words)], "Which")
    assert prepare(shown, 0, 1, 3)["kind"] == "view"
    assert prepare(shown, 0, 1, 3, scale=2.0)["kind"] is None
