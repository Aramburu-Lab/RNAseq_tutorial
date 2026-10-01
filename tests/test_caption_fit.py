"""A caption wraps inside its figure's width, at both column widths (2026-09-30).

Step-1's running-score caption, a gene-set name in capitals and its statistics, ran 0.07 in
off a double-column page: the wrap budgeted the full width, not the width inside the inset.
"""

import inspect

import pytest

from ctx.viz import annotate
from ctx.viz.save import layout_problems
from ctx.viz.style import CAPTION_INSET, SIZES, THEMES, WIDTHS_MM, apply_matplotlib

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
    assert "fig.text(CAPTION_INSET, " in inspect.getsource(annotate.place_caption)


@pytest.mark.parametrize("column", sorted(WIDTHS_MM))
def test_a_caption_line_ends_only_where_the_next_word_would_not_fit(column):
    """No gap of a word or more at the right: an average glyph width wrapped Rytis's 85 mm
    example in 5 lines, the widest 195 of 231 pt; measured, it takes 4 (2026-10-01).

    Mutation caught: the average-glyph wrap, which ends every line early."""
    from matplotlib.font_manager import FontProperties  # noqa: PLC0415 - heavy
    from matplotlib.textpath import text_to_path  # noqa: PLC0415 - heavy

    apply_matplotlib(THEMES["white"])
    usable = WIDTHS_MM[column] / 25.4 * 72 * (1 - 2 * CAPTION_INSET)
    lines = annotate.wrap_to_width(CAPTION, usable, SIZES["caption"]).split("\n")
    width = lambda line: text_to_path.get_text_width_height_descent(  # noqa: E731
        line, FontProperties(size=SIZES["caption"]), False
    )[0]
    assert " ".join(lines) == " ".join(CAPTION.split())
    assert all(width(line) <= usable for line in lines)
    for line, after in zip(lines, lines[1:], strict=False):
        assert width(f"{line} {after.split()[0]}") > usable, line
