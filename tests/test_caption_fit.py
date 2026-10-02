"""A caption wraps inside its figure's width, at both column widths (2026-09-30).

Step-1's running-score caption, a gene-set name in capitals and its statistics, ran 0.07 in
off a double-column page: the wrap budgeted the full width, not the width inside the inset.
"""

import inspect

import pytest

from ctx.viz import annotate
from ctx.viz.save import layout_problems
from ctx.viz.style import CAPTION_BOTTOM, CAPTION_INSET, THEMES, WIDTHS_MM, apply_matplotlib

CAPTION = (
    "HALLMARK_INTERFERON_ALPHA_RESPONSE in ndn_hsct_early: ES = 0.74, NES = 3.18, nominal "
    "p = 0.0018,\nFDR = 0 (no null score as extreme); 92 members of 15,109 ranked genes; "
    "1000 gene-set\npermutations, seed 42."
)


@pytest.mark.parametrize("column", sorted(WIDTHS_MM))
def test_caption_stays_on_the_page(column):
    from matplotlib.figure import Figure  # noqa: PLC0415 - heavy

    apply_matplotlib(THEMES["white"])
    width = WIDTHS_MM[column] / 25.4
    fig = Figure(figsize=(width, 3.0))
    fig.subplots().plot([0, 1], [0, 1])
    annotate.place_caption(fig, CAPTION, width, 3.0)
    assert layout_problems(fig) == []


def test_inset_is_the_one_place_caption_draws_at():
    # The insets are named once in style; a page's own PDF of a changed view reads them too.
    assert "fig.text(CAPTION_INSET, CAPTION_BOTTOM, " in inspect.getsource(annotate.place_caption)
    assert (CAPTION_INSET, CAPTION_BOTTOM) == (0.02, 0.01)
