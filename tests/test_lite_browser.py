"""ED.7 in Chromium: a core edition's lite page draws every figure, scrolls no sideways and
opens its PDF and SVG, at 1,280 and 375 px, light and dark (D-111).

A lite page carries no controls, so no figure is redrawn in the page and its PDF and SVG links
name the files the build wrote, in the theme shown; Chromium opens each one the links name.
Skipped where no Chromium is found; set ``CTX_CHROMIUM``.
"""

from dataclasses import replace
from datetime import date

import pytest
from test_export_browser import CHROMIUM, _browse
from test_lite import REPRODUCE
from test_report import _module
from test_viz import filtered_bars

from ctx import editions
from ctx.io.naming import stem
from ctx.report.builder import FigureText, Group, ResultTable, write
from ctx.viz.box import box_data
from ctx.viz.save import figure_block, write_figure

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
SIZES = {"wide_light": (1280, "light"), "wide_dark": (1280, "dark"),
         "phone_light": (375, "light"), "phone_dark": (375, "dark")}  # fmt: skip
DAY = date(2026, 10, 22)
# A file a link names, opened: its type as Chromium read it, and an SVG's drawn text.
OPENED = """({type: document.contentType,
           texts: document.querySelectorAll('text').length})"""
# What the page drew, after Plotly has had time, and the PDF and SVG links a reader sees.
LOOK = """new Promise(r => setTimeout(r, 2500)).then(() => {
  const shown = Array.from(document.querySelectorAll('.fig-files a'))
    .filter(a => a.offsetParent !== null && /[.](pdf|svg)$/.test(a.getAttribute('href')));
  return {links: shown.map(a => a.getAttribute('href')),
          figures: document.querySelectorAll('figure.fig').length,
          drawn: document.querySelectorAll('figure.fig .fig-plot .js-plotly-plot').length,
          forms: document.querySelectorAll('form').length,
          sections: Array.from(document.querySelectorAll('section')).map(s => s.id),
          page: document.documentElement.clientWidth,
          scroll: document.documentElement.scrollWidth,
          theme: document.documentElement.getAttribute('data-theme')};
})"""


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    """A core page of a bar figure that carried controls and a box figure, at four looks."""
    out = tmp_path_factory.mktemp("lite")
    bars = filtered_bars()
    bars_stem = stem("qc", "lite_depth", 1, day=DAY)
    words = FigureText("Depth per run", lead="A lead.", reading=("Each bar is one run.",))
    write_figure(bars, out, bars_stem, title="Depth", data_releases={}, text=words)
    from test_box_width import kidney_like  # noqa: PLC0415

    columns = {"value": ("log2 CPM", "x"), "feature": ("none", "Gene."),
               "group": ("none", "Group."), "sample": ("none", "Sample.")}  # fmt: skip
    box = box_data(kidney_like(), columns=columns, feature_label="Gene",
                   value_label="Expression", groups=("a", "b"))  # fmt: skip
    box_stem = stem("qc", "lite_box", 1, day=DAY)
    box_words = FigureText("Four genes by group", reading=("Each box is one gene.",))
    write_figure(box, out, box_stem, title="Box", data_releases={}, text=box_words)
    blocks = (replace(figure_block(bars, bars_stem, words), controls=bars.controls(),
                      reproduce=REPRODUCE), figure_block(box, box_stem, box_words))  # fmt: skip
    table = ResultTable("Calls per contrast", ("Contrast", "Up"), (("a", "1"),), results=True)
    page = _module(title="M2 QC: lite", results=(Group("Figures", (*blocks, table)),))
    mp = pytest.MonkeyPatch()
    mp.setenv(editions.EDITION_ENV, "core")
    try:
        module = write(page, out / f"{stem('qc', 'lite_overview', 1, day=DAY)}.html")
    finally:
        mp.undo()
    tasks = [
        {
            "name": name,
            "url": module.as_uri(),
            "width": width,
            "height": 900,
            "scheme": scheme,
            "act": LOOK,
        }
        for name, (width, scheme) in SIZES.items()
    ]
    files = [f"{s}_{t}.{e}" for s in (bars_stem, box_stem) for t in ("white", "black")
             for e in ("pdf", "svg")]  # fmt: skip
    tasks += [{"name": f, "url": (out / f).as_uri(), "width": 1280, "height": 900,
               "act": OPENED} for f in files]  # fmt: skip
    return _browse(tasks, tmp_path_factory.mktemp("downloads")), (bars_stem, box_stem), files


@pytest.mark.parametrize("name", list(SIZES))
def test_a_lite_page_draws_its_figures_and_links_their_files(ran, name):
    """Both figures drawn, no form, the results section alone, no sideways scroll, no console
    error, and each figure's PDF and SVG linked in the theme the page shows."""
    results, stems, _ = ran
    done = results[name]
    seen = done["acted"]
    assert done["errors"] == []
    assert (seen["figures"], seen["drawn"], seen["forms"]) == (2, 2, 0)
    assert seen["sections"] == ["results"]
    assert seen["scroll"] <= seen["page"]
    theme = "black" if SIZES[name][1] == "dark" else "white"
    assert seen["theme"] == theme
    assert seen["links"] == [f"{s}_{theme}.{e}" for s in stems for e in ("pdf", "svg")]


def test_every_linked_pdf_and_svg_opens(ran):
    """Chromium reads each file the links name as a PDF or an SVG, each SVG with its text."""
    results, _, files = ran
    for name in files:
        opened = results[name]
        assert opened["errors"] == [], name
        kind = "application/pdf" if name.endswith(".pdf") else "image/svg+xml"
        assert opened["acted"]["type"] == kind, name
        assert name.endswith(".pdf") or opened["acted"]["texts"] > 0, name
