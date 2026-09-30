"""Tick labels that collide on a short axis are thinned before the layout check (2026-09-30).

Track A's fastp strip held 990,000 to 1,000,000 reads in steps of 2,000 on a short value axis;
the labels stacked on each other and the page refused the figure.
"""

import io

from ctx.viz.save import layout_problems, thin_tick_labels
from ctx.viz.style import THEMES, apply_matplotlib


def _crowded():
    from matplotlib.figure import Figure  # noqa: PLC0415 - heavy
    from matplotlib.ticker import FixedLocator, StrMethodFormatter  # noqa: PLC0415 - heavy

    apply_matplotlib(THEMES["white"])
    fig = Figure(figsize=(3.35, 1.0))
    ax = fig.subplots()
    ax.plot([0, 1], [990_000, 1_000_000])
    ax.yaxis.set_major_locator(FixedLocator(range(990_000, 1_000_001, 1_000)))
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    fig.tight_layout()
    return fig, ax


def test_crowded_labels_are_thinned_and_stay_thinned():
    fig, ax = _crowded()
    assert any(p.startswith("overlap") for p in layout_problems(fig))
    thin_tick_labels(fig)
    assert layout_problems(fig) == []  # after a second draw, so the thinning survives it
    shown = [t.get_text() for t in ax.get_yticklabels() if t.get_visible()]
    assert 1 < len(shown) < 11 and shown[0] == "990,000"
    svg = io.StringIO()
    fig.savefig(svg, format="svg")
    # Index 1 is hidden at any step of 2 or more; the saved file draws what the check saw.
    assert "990,000" in svg.getvalue() and "991,000" not in svg.getvalue()


def test_labels_that_fit_are_left_alone():
    from matplotlib.figure import Figure  # noqa: PLC0415 - heavy

    apply_matplotlib(THEMES["white"])
    fig = Figure(figsize=(3.35, 3.0))
    ax = fig.subplots()
    ax.plot([0, 1], [0, 10])
    fig.tight_layout()
    thin_tick_labels(fig)
    assert all(t.get_visible() for t in ax.get_yticklabels() if t.get_text())


def test_an_axis_switched_off_has_no_labels_to_check():
    from matplotlib.figure import Figure  # noqa: PLC0415 - heavy

    apply_matplotlib(THEMES["white"])
    fig = Figure(figsize=(1.0, 1.0))
    ax = fig.subplots()
    ax.set_position((0.0, 0.0, 1.0, 1.0))  # the y labels sit left of the page
    assert any(p.startswith("off the page") for p in layout_problems(fig))
    ax.set_axis_off()  # now drawn nowhere, so nothing to report
    assert layout_problems(fig) == []


def test_long_level_names_wrap_inside_their_column():
    """The renamed batch levels ran into the next stage of the parallel sets (2026-09-30)."""
    import dataclasses  # noqa: PLC0415

    from synth import demo_linked  # noqa: PLC0415

    flow = demo_linked().flow()
    long = [dataclasses.replace(n, label=f"Reads ≤ 76 nt (min 35 nt) {n.label}")
            for n in flow.nodes]  # fmt: skip
    fig = dataclasses.replace(flow, nodes=tuple(long)).static(THEMES["white"])
    assert layout_problems(fig) == []
    assert any("\n" in t.get_text() for t in fig.axes[0].texts)


def test_no_detached_multiplier():
    """Step-1's rRNA fractions, near 1e-5 on a linear axis, got one "1e-5" above the axis."""
    from matplotlib.figure import Figure  # noqa: PLC0415 - heavy

    apply_matplotlib(THEMES["white"])
    fig = Figure(figsize=(3.35, 2.5))
    ax = fig.subplots()
    ax.plot([0, 1], [6e-7, 6.8e-5])
    fig.tight_layout()  # the layout sees the full labels, so they stay on the page
    from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415 - heavy

    fig.draw(FigureCanvasAgg(fig).get_renderer())
    assert ax.yaxis.get_offset_text().get_text() == ""
    assert any(t.get_text().startswith("0.0000") for t in ax.get_yticklabels())
    assert layout_problems(fig) == []
