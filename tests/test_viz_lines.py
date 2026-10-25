"""One line per sample: a unique colour and dash each, three times a plot line, legible on both
themes (Rytis 2026-09-30, the base-quality curves)."""

import polars as pl
import pytest

from ctx.viz.curve import CurveData
from ctx.viz.style import (
    LINE_COLOURS,
    LINE_WIDTHS,
    MAX_DISTINCT,
    PAGE_SAMPLE_LINE_SCALE,
    THEMES,
    distinct_styles,
)


def _curves(n: int) -> CurveData:
    frame = pl.DataFrame({
        "run": [f"r{i:02d}" for i in range(n) for _ in range(3)],
        "position": [1, 2, 3] * n,
        "quality": [30.0 + i for i in range(n) for _ in range(3)],
    })  # fmt: skip
    columns = {"run": ("none", "Run."), "position": ("nt", "Base."), "quality": ("Phred", "Q.")}
    return CurveData(frame=frame, columns=columns, series="run", x="position", y="quality",
                     x_label="Position", y_label="Quality", colour="run", colours={},
                     caption="Made up.", distinct=True)  # fmt: skip


def _luminance(hexa: str) -> float:
    rgb = [int(hexa[i : i + 2], 16) / 255 for i in (1, 3, 5)]
    linear = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def test_every_run_has_its_own_style_in_both_renderers():
    data = _curves(9)
    lines = data.static(THEMES["white"]).axes[0].get_lines()
    drawn = [(line.get_color(), line.get_linestyle()) for line in lines]
    assert len(set(drawn)) == 9
    assert {line.get_linewidth() for line in lines} == {LINE_WIDTHS["sample"]}
    assert LINE_WIDTHS["sample"] == 3 * LINE_WIDTHS["line"]
    traces = [t for t in data.interactive().data if t.mode == "lines"]
    assert len({(t.line.color, t.line.dash) for t in traces}) == 9
    # One line per sample draws at twice the files' width on the page (Rytis 2026-10-25).
    assert {t.line.width for t in traces} == {LINE_WIDTHS["sample"] * PAGE_SAMPLE_LINE_SCALE}
    assert {t.opacity for t in traces} == {1.0}


def test_more_runs_than_styles_are_refused():
    assert len(distinct_styles([str(i) for i in range(MAX_DISTINCT)])) == MAX_DISTINCT
    with pytest.raises(ValueError, match=f"at most {MAX_DISTINCT}"):
        _curves(MAX_DISTINCT + 1)


@pytest.mark.parametrize("theme", ["white", "black"])
def test_line_colours_reach_three_to_one(theme):
    background = _luminance(THEMES[theme].background)
    for hexa in LINE_COLOURS:
        high, low = sorted((_luminance(hexa), background), reverse=True)
        assert (high + 0.05) / (low + 0.05) >= 3.0, (hexa, theme)
