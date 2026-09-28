"""The only place a colour, font, size, line width or figure dimension is set.

Order and font follow CONVENTIONS.md (Figures); line widths, figure widths, black as ink
and the seven-category cap follow the global figure rules where the project is silent.
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

FONT_CHAIN = ("Helvetica", "Arial", "Liberation Sans", "DejaVu Sans")
SIZES = {"base": 14.0, "axis_label": 14.0, "tick": 12.0, "annotation": 12.0, "caption": 12.0}
LINE_WIDTHS = {"axes": 1.0, "line": 1.0, "box": 1.5, "errorbar": 1.5, "marker_edge": 0.0}
MARKER_SIZE = 4.0
# Plotly marker diameter in pixels, the on-screen size of MARKER_SIZE.
PLOTLY_MARKER_PX = 6.0
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
# Above this many features a box plot takes the double column: at 14 pt the labels of more
# features do not fit in 85 mm (measured 2026-09-26, 8 genes: every label overlapped).
SINGLE_COLUMN_MAX_FEATURES = 3
# Empty space above the tallest box for the p/r labels, as a fraction of the data range.
HEADROOM = 0.35
LINE_SPACING = 1.3
# Legend geometry in font-size units; marker handles need less than matplotlib's line default.
LEGEND = {"handlelength": 1.0, "columnspacing": 1.0, "handletextpad": 0.4}
LEGEND_MAX_COLUMNS = 4
# Category labels under vertical bars run upward, so long names never collide.
BAR_LABEL_ROTATION = 90.0
# Horizontal bars grow the figure: this much height per category, plus room for the axis,
# legend and caption, so 20 categories stay readable at 14 pt.
BAR_ROW_MM = 7.0
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
        "white", ink=INK_COLOUR, background="#FFFFFF", grid="#E5E5E5", reference="#8C8C8C"
    ),
    "black": Theme(
        "black", ink="#F2F2F2", background="#111111", grid="#333333", reference="#8C8C8C"
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
# Everything read on screen, page text and the plots on the page, runs 3 pt larger than the
# figure sizes (Rytis, 2026-09-28): 3 pt is 4 px at CSS's 96 px to 72 pt. The PDF and SVG
# files keep SIZES, which their 85 mm layout checks are measured at.
PAGE_BOOST_PX = 3 * 96 / 72
# Page type scale in px around SIZES["base"], before the boost.
TYPE_PX = {"xs": 12, "sm": 13, "md": 16, "lg": 20, "xl": 28, "display": 44}


def on_page(size: float) -> float:
    """A text size as a page shows it: the figure size plus ``PAGE_BOOST_PX``."""
    return size + PAGE_BOOST_PX


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
    """Set matplotlib's rc for a theme: seaborn's paper context first, then this file's values.

    The order matters: ``set_context("paper")`` rescales fonts and line widths, so applying
    it afterwards would silently replace every size set here.

    Args:
        theme: One of ``THEMES``.
        fonts: The font chain; tests pass a chain with a missing first font.
    """
    import matplotlib as mpl  # noqa: PLC0415 - heavy; imported after the runtime profile
    import seaborn as sns  # noqa: PLC0415 - heavy; imported after the runtime profile
    from cycler import cycler  # noqa: PLC0415 - ships with matplotlib

    sns.set_context("paper")
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": list(fonts),
            "font.size": SIZES["base"],
            "axes.labelsize": SIZES["axis_label"],
            "axes.titlesize": SIZES["axis_label"],
            "xtick.labelsize": SIZES["tick"],
            "ytick.labelsize": SIZES["tick"],
            "legend.fontsize": SIZES["tick"],
            "legend.title_fontsize": SIZES["tick"],
            "axes.linewidth": LINE_WIDTHS["axes"],
            "xtick.major.width": LINE_WIDTHS["axes"],
            "ytick.major.width": LINE_WIDTHS["axes"],
            "lines.linewidth": LINE_WIDTHS["line"],
            "patch.linewidth": LINE_WIDTHS["box"],
            "lines.markeredgewidth": LINE_WIDTHS["marker_edge"],
            "text.color": theme.ink,
            "axes.labelcolor": theme.ink,
            "axes.edgecolor": theme.ink,
            "xtick.color": theme.ink,
            "ytick.color": theme.ink,
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
        "--size-base": f"{on_page(SIZES['base']):g}px",
    }
    return "; ".join(f"{k}: {v}" for k, v in pairs.items()) + ";"


def css_layout() -> str:
    """Type scale, spacing and layout as CSS custom properties; the same in both themes."""
    pairs = {
        **{f"--fs-{name}": f"{on_page(px):g}px" for name, px in TYPE_PX.items()},
        **{f"--space-{i + 1}": f"{px}px" for i, px in enumerate(SPACE_PX)},
        **{f"--{name}": value for name, value in LAYOUT.items()},
        "--font-mono": ", ".join(f'"{f}"' if " " in f else f for f in MONO_CHAIN),
    }
    return "; ".join(f"{k}: {v}" for k, v in pairs.items()) + ";"
