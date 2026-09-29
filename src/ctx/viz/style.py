"""The only place a colour, font, size, line width or figure dimension is set.

Sizes, fonts, line widths, ticks and points are those of the Candida Ribo-seq figures
(Rytis, 2026-09-29, D-024), which are the global figure rules; the palette order follows
CONVENTIONS.md (Figures).
"""

from dataclasses import dataclass
from typing import Any, Literal

MM_PER_INCH = 25.4
# sRGB transfer function breakpoint (WCAG 2 relative luminance).
SRGB_LINEAR_BELOW = 0.04045

# Okabe-Ito in the CONVENTIONS.md order. #000000, the eighth, is ink, never a category.
DATA_COLOURS = ("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
INK_COLOUR = "#000000"
MAX_CATEGORIES = len(DATA_COLOURS)
# 1.32:1 against white (WCAG 2 contrast), the lowest of the seven; unreadable when projected.
PROJECTOR_UNSAFE = frozenset({"#F0E442"})
# Colours by role, never by list position, so a legend cannot drift from its data.
ROLES = {"reference": "#0072B2", "comparison": "#D55E00"}

# Liberation Sans is Arial's metric twin, which a machine without Arial may redistribute.
FONT_CHAIN = ("Arial", "Liberation Sans", "DejaVu Sans")
SIZES = {
    "base": 9.0, "axis_label": 10.0, "tick": 9.0, "legend": 9.0, "annotation": 9.0,
    "caption": 9.0, "panel": 12.0,
}  # fmt: skip
# No text below this; a mathtext superscript on a 9 pt tick would be 6.3 pt, so none is drawn.
MIN_TEXT_PT = 8.0
LINE_WIDTHS = {"axes": 1.0, "line": 1.0, "box": 1.5, "errorbar": 1.5, "marker_edge": 0.0}
TICK_PT = {"major": 3.0, "minor": 1.5}
# Point diameter in pt, and its opacity, so overplotted points still show their density.
MARKER_SIZE = 6.0
POINT_ALPHA = 0.5
PX_PER_PT = 96 / 72
# Plotly marker diameter in pixels: MARKER_SIZE on screen.
PLOTLY_MARKER_PX = MARKER_SIZE * PX_PER_PT
# Panel letters: bold capitals at the top left, nudged up by dy_pt (no figure has panels yet).
PANEL_LABEL = {"size": SIZES["panel"], "weight": "bold", "dx_pt": 0.0, "dy_pt": 2.0}
WIDTHS_MM = {"single": 85.0, "double": 170.0}
ASPECT = {"box": 1.1, "scatter": 1.1, "bar": 1.0}
# Point shapes for a second category, in level order: (matplotlib, plotly).
MARKERS = (
    ("o", "circle"),
    ("^", "triangle-up"),
    ("s", "square"),
    ("D", "diamond"),
    ("v", "triangle-down"),
)
# Average glyph width in em, used to wrap the caption to the width a figure really has.
AVERAGE_CHAR_EM = 0.55
# Above this many features a box plot takes the double column: at 9 pt the p and r labels of
# 5 features overlap in 85 mm and those of 4 do not (measured 2026-09-29).
SINGLE_COLUMN_MAX_FEATURES = 4
# Empty space above the tallest box for the p/r labels, as a fraction of the data range.
HEADROOM = 0.35
LINE_SPACING = 1.3
# Legend geometry in font-size units; marker handles need less than matplotlib's line default.
LEGEND = {"handlelength": 1.0, "columnspacing": 1.0, "handletextpad": 0.4}
LEGEND_MAX_COLUMNS = 4
# Category labels under vertical bars run upward, so long names never collide.
BAR_LABEL_ROTATION = 90.0
# Horizontal bars grow the figure: this much height per category, plus room for the axis,
# legend and caption. 20 categories at 9 pt stop overlapping from 3.0 mm (measured
# 2026-09-29); 4.5 mm holds a label's line pitch, 9 pt x LINE_SPACING = 4.1 mm, with a gap.
BAR_ROW_MM = 4.5
# The share of a figure's height a vertical axis title may use before it wraps.
AXIS_TITLE_SHARE = 0.55
BAR_MARGIN_MM = 45.0


@dataclass(frozen=True)
class Theme:
    """Colours of one theme; the name is the file token (``_white`` / ``_black``)."""

    name: str
    ink: str
    background: str
    grid: str
    reference: str


# One reference grey for both themes: a plotly trace keeps its colour when the page theme
# flips, so it must read on both backgrounds (3.36:1 on white, 5.62:1 on #111111; the 3:1
# graphics threshold). The earlier light grey, #999999, was 2.85:1.
THEMES = {
    "white": Theme(
        "white", ink=INK_COLOUR, background="#FFFFFF", grid="#E6E6E6", reference="#8C8C8C"
    ),
    "black": Theme(
        "black", ink="#F2F2F2", background="#111111", grid="#39424A", reference="#8C8C8C"
    ),
}
# Page furniture per theme, never a data colour, so a link or a pill cannot read as a series.
# Contrast on the background (surface): text 17.40 (16.07), muted 6.69 (6.18), link 8.39 (7.75)
# light; 14.85 (13.40), 7.94 (7.17), 9.20 (8.30) dark; each pill's ink on its fill >= 6.83.
CHROME = {
    "white": {
        "surface": "#F6F6F3", "text": "#1A1A1A", "muted": "#5C5C5C", "link": "#1D4E89",
        "pass-fg": "#1B5E20", "pass-bg": "#E6F2E7", "fail-fg": "#8E1B1B", "fail-bg": "#FBE9E7",
        "recorded-fg": "#4D4D4D", "recorded-bg": "#EDEDED",
    },
    "black": {
        "surface": "#1C1C1C", "text": "#E4E4E4", "muted": "#A8A8A8", "link": "#8CB8F0",
        "pass-fg": "#9BD3A5", "pass-bg": "#16301C", "fail-fg": "#F4ABA3", "fail-bg": "#3D1816",
        "recorded-fg": "#C8C8C8", "recorded-bg": "#262626",
    },
}  # fmt: skip
MONO_CHAIN = ("ui-monospace", "SF Mono", "Menlo", "Consolas", "Liberation Mono", "monospace")
# The plots on a page run 3 pt larger than the figure sizes (Rytis, 2026-09-28). The PDF and
# SVG files keep SIZES, which their 85 mm layout checks are measured at.
PAGE_BOOST_PT = 3.0
# Page text in px, set on its own: it no longer follows the figure sizes.
PAGE_BASE_PX = 18
TYPE_PX = {"xs": 16, "sm": 17, "md": 20, "lg": 24, "xl": 32, "display": 48}


def on_page(size: float) -> float:
    """A figure text size in pt as a plot on the page shows it, in px: 3 pt larger."""
    return (size + PAGE_BOOST_PT) * PX_PER_PT


SPACE_PX = (4, 8, 12, 16, 24, 32, 48, 64, 96)
LAYOUT = {
    "content": "1040px", "rail": "208px", "gap": "56px", "gutter": "16px",
    "gutter-wide": "32px", "key-min": "150px", "card-min": "240px", "reading": "320px",
    "radius": "4px",
    "rule": "1px", "rule-mid": "2px", "rule-strong": "4px",
}  # fmt: skip
# Media queries cannot read var(), so the breakpoints are injected into the stylesheet.
RAIL_FROM_PX, NARROW_PX = 1100, 640
PLOTLY_HEIGHT_PX = {"scatter": 460, "box": 440, "bar": 400}
# Horizontal bars grow per category, as the static figure does with BAR_ROW_MM.
PLOTLY_BAR_BASE_PX, PLOTLY_BAR_ROW_PX = 150, 28
PLOTLY_MARGIN_PX = {"l": 8, "r": 8, "t": 48, "b": 8}  # automargin grows each side as needed
PLOTLY_TICK_PX, PLOTLY_TITLE_STANDOFF_PX = 4, 8
TRANSPARENT = "rgba(0,0,0,0)"
# Dash of a threshold line, in both renderers.
THRESHOLD_DASH: tuple[Literal["--"], Literal["dash"]] = ("--", "dash")


def category_colours(levels: list[str]) -> dict[str, str]:
    """Colour each category level, in the given order, from the data colours.

    Args:
        levels: Category levels in display order.

    Returns:
        Level to hex colour.

    Raises:
        ValueError: Above ``MAX_CATEGORIES``: an eighth colour would have to be ink.
    """
    if len(levels) > MAX_CATEGORIES:
        raise ValueError(f"{len(levels)} categories; at most {MAX_CATEGORIES} can be coloured")
    return dict(zip(levels, DATA_COLOURS, strict=False))


def figure_inches(width: str, kind: str) -> tuple[float, float]:
    """Figure size in inches for a column width and a figure type's aspect ratio."""
    inches = WIDTHS_MM[width] / MM_PER_INCH
    return inches, WIDTHS_MM["single"] / MM_PER_INCH * ASPECT[kind]


def caption_chars(width_inches: float) -> int:
    """How many caption characters fit on one line of a figure this wide."""
    return int(width_inches * 72 / (SIZES["caption"] * AVERAGE_CHAR_EM))


def contrast(first: str, second: str) -> float:
    """WCAG 2 contrast ratio between two hex colours, from 1 to 21."""

    def luminance(colour: str) -> float:
        """Relative luminance of a hex colour (WCAG 2 definition)."""
        channels = [int(colour[i : i + 2], 16) / 255 for i in (1, 3, 5)]
        linear = [
            c / 12.92 if c <= SRGB_LINEAR_BELOW else ((c + 0.055) / 1.055) ** 2.4 for c in channels
        ]
        return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]

    high, low = sorted((luminance(first), luminance(second)), reverse=True)
    return (high + 0.05) / (low + 0.05)


def apply_matplotlib(theme: Theme, fonts: tuple[str, ...] = FONT_CHAIN) -> None:
    """Set matplotlib's rc for a theme from this file's values alone.

    No seaborn context: ``set_context("paper")`` also set 4.8 pt ticks, 0.8 pt minor ticks
    and 4.8 pt markers that no value here named.

    Args:
        theme: One of ``THEMES``.
        fonts: The font chain; tests pass a chain with a missing first font.
    """
    import matplotlib as mpl  # noqa: PLC0415 - heavy; imported after the runtime profile
    from cycler import cycler  # noqa: PLC0415 - ships with matplotlib

    mpl.rcParams.update(
        {
            "xtick.major.size": TICK_PT["major"],
            "ytick.major.size": TICK_PT["major"],
            "xtick.minor.size": TICK_PT["minor"],
            "ytick.minor.size": TICK_PT["minor"],
            "xtick.major.width": LINE_WIDTHS["axes"],
            "ytick.major.width": LINE_WIDTHS["axes"],
            "xtick.minor.width": LINE_WIDTHS["axes"],
            "ytick.minor.width": LINE_WIDTHS["axes"],
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.labelsize": SIZES["tick"],
            "ytick.labelsize": SIZES["tick"],
            "xtick.color": theme.ink,
            "ytick.color": theme.ink,
            "font.family": "sans-serif",
            "font.sans-serif": list(fonts),
            "font.size": SIZES["base"],
            "axes.labelsize": SIZES["axis_label"],
            "axes.titlesize": SIZES["axis_label"],
            "legend.fontsize": SIZES["legend"],
            "legend.title_fontsize": SIZES["legend"],
            "axes.linewidth": LINE_WIDTHS["axes"],
            "axes.formatter.use_mathtext": False,
            "lines.linewidth": LINE_WIDTHS["line"],
            "lines.markersize": MARKER_SIZE,
            "grid.color": theme.grid,
            "grid.linewidth": LINE_WIDTHS["axes"],
            "patch.linewidth": LINE_WIDTHS["box"],
            "lines.markeredgewidth": LINE_WIDTHS["marker_edge"],
            "text.color": theme.ink,
            "axes.labelcolor": theme.ink,
            "axes.edgecolor": theme.ink,
            "axes.prop_cycle": cycler(color=list(DATA_COLOURS)),
            "figure.facecolor": "none",
            "axes.facecolor": "none",
            "savefig.transparent": True,
            "legend.frameon": False,
            "legend.handlelength": LEGEND["handlelength"],
            "legend.columnspacing": LEGEND["columnspacing"],
            "legend.handletextpad": LEGEND["handletextpad"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def _plotly_axis(theme: Theme) -> dict[str, Any]:
    """Axes as the static figures draw them: ink lines and ticks, no grid, room for titles."""
    small = {"size": on_page(SIZES["tick"]), "color": theme.ink}
    return {
        "showgrid": False, "showline": True, "linecolor": theme.ink,
        "linewidth": LINE_WIDTHS["axes"], "ticks": "outside", "ticklen": PLOTLY_TICK_PX,
        "tickcolor": theme.ink, "tickwidth": LINE_WIDTHS["axes"], "tickfont": small,
        "zeroline": False, "automargin": True,
        "title": {"font": {"size": on_page(SIZES["axis_label"]), "color": theme.ink},
                  "standoff": PLOTLY_TITLE_STANDOFF_PX},
    }  # fmt: skip


def plotly_templates() -> dict[str, Any]:
    """Plotly layout templates for the two themes, built on ``plotly_white``.

    Backgrounds are transparent, so a figure takes the page's colour in either theme; the
    page's toggle swaps these templates with ``Plotly.relayout``.
    """
    import plotly.io as pio  # noqa: PLC0415 - heavy; imported after the runtime profile

    templates = {}
    for name, theme in THEMES.items():
        template = pio.templates["plotly_white"].to_plotly_json()
        small = {"size": on_page(SIZES["tick"]), "color": theme.ink}
        template["layout"].update(
            {
                "font": {
                    "family": ", ".join([*FONT_CHAIN, "sans-serif"]),
                    "size": on_page(SIZES["base"]),
                    "color": theme.ink,
                },
                "colorway": list(DATA_COLOURS),
                "paper_bgcolor": TRANSPARENT,
                "plot_bgcolor": TRANSPARENT,
                "margin": dict(PLOTLY_MARGIN_PX),
                "legend": {"orientation": "h", "x": 0, "xanchor": "left", "y": 1.02,
                           "yanchor": "bottom", "font": small, "bgcolor": TRANSPARENT,
                           "title": {"text": ""}},
                "hoverlabel": {"bgcolor": theme.background, "bordercolor": theme.grid,
                               "font": {**small, "family": ", ".join(FONT_CHAIN)},
                               "align": "left"},
                "modebar": {"bgcolor": TRANSPARENT, "color": theme.reference,
                            "activecolor": theme.ink},
                "annotationdefaults": {"font": small, "arrowcolor": theme.ink},
                "shapedefaults": {"line": {"color": theme.reference}},
                "xaxis": _plotly_axis(theme),
                "yaxis": _plotly_axis(theme),
            }
        )  # fmt: skip
        templates[name] = template
    return templates


def css_variables(theme: Theme) -> str:
    """The theme as CSS custom properties, so the HTML and the figures share one source."""
    pairs = {
        "--ink": theme.ink,
        "--background": theme.background,
        "--grid": theme.grid,
        "--reference": theme.reference,
        **{f"--data-{i + 1}": colour for i, colour in enumerate(DATA_COLOURS)},
        **{f"--role-{role}": colour for role, colour in ROLES.items()},
        **{f"--{key}": colour for key, colour in CHROME[theme.name].items()},
        "--font": ", ".join([*(f'"{f}"' for f in FONT_CHAIN), "sans-serif"]),
        "--size-base": f"{PAGE_BASE_PX}px",
    }
    return "; ".join(f"{k}: {v}" for k, v in pairs.items()) + ";"


def css_layout() -> str:
    """Type scale, spacing and layout as CSS custom properties; the same in both themes."""
    pairs = {
        **{f"--fs-{name}": f"{px}px" for name, px in TYPE_PX.items()},
        **{f"--space-{i + 1}": f"{px}px" for i, px in enumerate(SPACE_PX)},
        **{f"--{name}": value for name, value in LAYOUT.items()},
        "--font-mono": ", ".join(f'"{f}"' if " " in f else f for f in MONO_CHAIN),
    }
    return "; ".join(f"{k}: {v}" for k, v in pairs.items()) + ";"
