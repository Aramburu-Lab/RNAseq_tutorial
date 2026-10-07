"""A caption wraps inside its figure's width, at both column widths (2026-09-30).

Step-1's running-score caption, a gene-set name in capitals and its statistics, ran 0.07 in
off a double-column page: the wrap budgeted the full width, not the width inside the inset.
"""

import inspect

import pytest

from ctx.viz import annotate
from ctx.viz.metrics import text_width
from ctx.viz.save import layout_problems
from ctx.viz.style import (
    CAPTION_BOTTOM,
    CAPTION_INSET,
    SIZES,
    THEMES,
    WIDTHS_MM,
    apply_matplotlib,
)

CAPTION = (
    "HALLMARK_INTERFERON_ALPHA_RESPONSE in ndn_hsct_early: ES = 0.74, NES = 3.18, nominal "
    "p = 0.0018,\nFDR = 0 (no null score as extreme); 92 members of 15,109 ranked genes; "
    "1000 gene-set\npermutations, seed 42."
)
# Step-1's PCA caption, one of the course's 99 (2026-10-07): at 85 mm the advances without
# their kerning ended its third line a word early, "variable" still fitting as matplotlib sets it.
KERNED = (
    "PERMANOVA on Euclidean distances for condition: R² = 0.40, pseudo-F = 8.28, p = 0.001 "
    "(999 permutations; n = 55 samples in 5 groups). Blind VST, the 500 most variable genes "
    "after the transform; PC1 33%, PC2 9% of their variance; n = 55 samples."
)
# How far matplotlib's set width of a caption line may lie from ctx.viz.metrics' measure: on
# the course's captions it set every line 0.0004 to 0.038 pt narrower at 170 mm, and from
# 0.020 pt narrower to 0.0016 pt wider at 85 mm (FreeType's fixed-point rounding).
SET_SLACK_PT = 0.05


def _drawn(caption: str, column: str) -> tuple[list[str], float, float]:
    """The caption's lines as place_caption draws them, the room inside the insets, and the
    figure's caption size."""
    from matplotlib.figure import Figure  # noqa: PLC0415 - heavy

    apply_matplotlib(THEMES["white"])
    width = WIDTHS_MM[column] / 25.4
    fig = Figure(figsize=(width, 3.0))
    fig.subplots().plot([0, 1], [0, 1])
    annotate.place_caption(fig, caption, width, 3.0)
    (text,) = fig.texts
    return text.get_text().split("\n"), width * 72 * (1 - 2 * CAPTION_INSET), text.get_fontsize()


def _set_width(line: str, size: float) -> float:
    """A line's width as matplotlib sets it in the PDF, in pt."""
    from matplotlib.font_manager import FontProperties  # noqa: PLC0415 - heavy
    from matplotlib.textpath import text_to_path  # noqa: PLC0415 - heavy

    return text_to_path.get_text_width_height_descent(line, FontProperties(size=size), False)[0]


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


@pytest.mark.parametrize("caption", [CAPTION, KERNED], ids=["running_score", "pca"])
@pytest.mark.parametrize("column", sorted(WIDTHS_MM))
def test_a_caption_line_ends_only_where_the_next_word_would_not_fit(column, caption):
    """No gap of a word or more at the right, and no line past the inset, as matplotlib sets
    the lines the figure draws: an average glyph width wrapped Rytis's 85 mm example in 5
    lines, the widest 195 of 231 pt; measured, it takes 4 (2026-10-01).

    Mutations caught: the average-glyph wrap, which ends every line early; the measure without
    its kerning, which ended the PCA caption's third line a word early at 85 mm."""
    lines, usable, size = _drawn(caption, column)
    assert size == SIZES["caption"]
    assert " ".join(lines) == " ".join(caption.split())
    assert all(_set_width(line, size) <= usable for line in lines)
    for line, after in zip(lines, lines[1:], strict=False):
        assert _set_width(f"{line} {after.split()[0]}", size) > usable, line


def test_a_machine_without_arial_sets_the_files_in_the_face_measured():
    """Without Arial the files are set in Liberation Sans, the face ctx.viz.metrics measures,
    from the vendored file where the machine has none. Mutation: the vendored faces never
    given to matplotlib, which on a machine without Liberation Sans fell back to DejaVu Sans,
    wider than the measure."""
    from matplotlib.font_manager import FontProperties, findfont  # noqa: PLC0415 - heavy

    apply_matplotlib(THEMES["white"], fonts=("NoSuchFontXYZ", "Liberation Sans", "DejaVu Sans"))
    try:
        face = findfont(FontProperties(size=SIZES["caption"]))
        assert face.endswith("LiberationSans-Regular.ttf"), face
    finally:
        apply_matplotlib(THEMES["white"])


@pytest.mark.parametrize("caption", [CAPTION, KERNED], ids=["running_score", "pca"])
def test_the_measure_is_the_width_matplotlib_sets(caption):
    """ctx.viz.metrics, which the figure files, the page's wrapped labels and a page's own PDF
    of a changed view all wrap by, against matplotlib's set width of the same text, word by
    word up to the widest caption line, 170 mm inside the insets. Mutation caught: the kerning
    left out, up to 4.33 pt wide on the course's captions."""
    apply_matplotlib(THEMES["white"])
    widest = WIDTHS_MM["double"] / 25.4 * 72 * (1 - 2 * CAPTION_INSET)
    words, checked = caption.split(), 0
    for end in range(1, len(words) + 1):
        line = " ".join(words[:end])
        if text_width(line, SIZES["caption"]) > widest:
            break
        measured, drawn = text_width(line, SIZES["caption"]), _set_width(line, SIZES["caption"])
        assert measured == pytest.approx(drawn, abs=SET_SLACK_PT), line
        checked += 1
    assert checked >= 10
