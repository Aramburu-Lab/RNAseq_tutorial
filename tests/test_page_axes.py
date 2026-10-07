"""A page's axes: its y title has the figure files' lines, from the one ``fit_title`` call (item
2, D10 and D15), and a strip shows every stratum's label.

``axes.screen_title`` broke every page title over 26 characters before its unit, a rule for
x titles on a 375 px phone; the page script undid it for x titles alone, and the page's own PDF
and SVG of a moved view copy the page's y title, so they carried that break too.
"""

import dataclasses

import pytest
from synth import (
    counts,
    demo_bar,
    demo_bars_at_scale,
    demo_box,
    demo_bubble,
    demo_curves,
    demo_density,
    demo_ma,
    demo_pca,
    demo_strip,
    demo_volcano,
)

from ctx.viz.axes import screen_title
from ctx.viz.style import THEMES

# A title fit_title and screen_title break apart: after "share" in the files, before "(%)" by
# the phone rule.
LONG, UNIT = "Mitochondrial rRNA share of reads", "%"


def _long(data):
    """The demo with LONG as its y title, in its unit."""
    if hasattr(data, "unit") and hasattr(data, "groups"):  # a box: its unit is a field
        return dataclasses.replace(data, value_label=LONG, unit=UNIT)
    field, column = next((f, getattr(data, c)) for f, c in (("y_label", "y"),
                                                            ("value_label", "value"))
                         if hasattr(data, f))  # fmt: skip
    meaning = data.columns[column][1]
    return dataclasses.replace(data, **{field: LONG},
                               columns={**data.columns, column: (UNIT, meaning)})  # fmt: skip


DEMOS = {
    "volcano": lambda: demo_volcano(counts()),
    "ma": lambda: demo_ma(counts()),
    "pca": lambda: demo_pca(counts()),
    "bubble": demo_bubble,
    "bar": lambda: demo_bar(counts()),
    "bar across": lambda: demo_bar(counts(), horizontal=True),
    "bars at scale": lambda: demo_bars_at_scale(False),
    "box": lambda: demo_box(counts()),
    "curves": demo_curves,
    "density": demo_density,
    "strip": demo_strip,
    "strip, log": lambda: demo_strip(log=True),
}


@pytest.mark.parametrize("long", [False, True], ids=["as built", "long title"])
@pytest.mark.parametrize("name", list(DEMOS))
def test_the_page_y_title_is_the_files(name, long):
    """Line by line, the page's y title is the white figure file's y label: a horizontal bar's
    category title as well as every value title. Mutation: screen_title applied to a y title,
    which breaks the long title before its unit."""
    data = DEMOS[name]()
    data = _long(data) if long else data
    files = data.static(THEMES["white"]).axes[0].get_ylabel().split("\n")
    page = data.interactive().layout.yaxis.title.text.split("<br>")
    assert page == files
    if long and not getattr(data, "horizontal", False):
        assert screen_title(" ".join(files)).split("<br>") != files, "the case must tell"


@pytest.mark.parametrize("log", [False, True])
def test_a_strip_shows_every_stratum_label(log):
    """A strip's strata are categories on a linear axis, which plotly gives "hide past div": it
    hid a label that, turned, ended 12 px inside the plot. Measured in Chrome on Step-1's 15 M2
    strip pages: 9 hid their last label at 1,280 px, 12 once the y titles took two lines, 0
    with "allow", none of whose labels then ran past the plot. Mutation: the setting dropped."""
    axis = demo_strip(log=log).interactive().layout.xaxis
    assert axis.ticklabeloverflow == "allow"
    assert list(axis.ticktext) == ["PB / NDN", "PB / Mono", "BM / NDN"]
