"""P2S.8: a vertical title wraps to its own axis's drawn height, in every contract.

``AXIS_TITLE_SHARE`` gave every figure a scatter's share of its height, longer than some
strips' and bars' axes: "Share of the top gene (fraction)" ran 2.5 pt past its axis (package E,
2026-10-07). Each line is measured as matplotlib sets it (``test_caption_fit._set_width``),
against the axes' height as laid out in the file, never through ``ctx.viz.metrics``.
"""

import pytest
from synth import demo_surface
from test_caption_fit import _set_width
from test_page_axes import DEMOS, _long

from ctx.viz.metrics import words
from ctx.viz.style import AXIS_TITLE_SHARE, SIZES, THEMES

CASES = {**DEMOS, "surface": demo_surface}


def _height(fig, ax) -> float:
    """The axes' height in pt, read off matplotlib's own layout."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415 - heavy

    renderer = FigureCanvasAgg(fig).get_renderer()
    return ax.get_window_extent(renderer).height * 72 / fig.dpi


def _titles(fig, data):
    """(axes, its vertical title's lines) of the figure: the y axis, and a colour bar's."""
    found = [] if getattr(data, "horizontal", False) else [(fig.axes[0], fig.axes[0].get_ylabel())]
    found += [(ax, ax.get_ylabel()) for ax in fig.axes[1:] if ax.get_ylabel()]
    return [(ax, label.split("\n")) for ax, label in found]


@pytest.mark.parametrize("long", [False, True], ids=["as built", "long title"])
@pytest.mark.parametrize("name", list(CASES))
def test_a_vertical_title_fits_its_own_drawn_axis(name, long):
    """Every line fits the axis it runs along, and ends only where its next word would cross
    it. Mutation: ``refit_y_title`` left out of a contract's static (the share's estimate
    kept): the long title ends early on axes taller than the estimate, or runs past shorter
    ones."""
    data = CASES[name]()
    data = _long(data) if long and name != "surface" else data
    fig = data.static(THEMES["white"])
    size = SIZES["axis_label"]
    for ax, lines in _titles(fig, data):
        room = _height(fig, ax)
        for line in lines:
            assert _set_width(line, size) <= room + 0.01 or len(words(line)) == 1, (line, room)
        for line, after in zip(lines, lines[1:], strict=False):
            assert _set_width(f"{line} {words(after)[0]}", size) > room, (line, after, room)


def test_the_share_was_not_the_axis():
    """The rule changed something: on some contract the drawn axis differs from the share of
    the figure's height by more than a line's worth of a word, which is why a title ran past
    its axis before (stated apart from the code: 0.55 of the height)."""
    gaps = {}
    for name, make in CASES.items():
        data = make()
        fig = data.static(THEMES["white"])
        height = fig.get_size_inches()[1] * 72
        gaps[name] = _height(fig, fig.axes[0]) - 0.55 * height
    assert AXIS_TITLE_SHARE == 0.55
    assert any(abs(g) > 10 for g in gaps.values()), gaps
