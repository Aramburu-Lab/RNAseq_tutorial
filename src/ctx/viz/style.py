"""The only place a colour, font, size, line width or figure dimension is set.

Sizes, fonts, line widths, ticks and points are those of the Candida Ribo-seq figures
(Rytis, 2026-09-29, D-024), which are the global figure rules; the palette order follows
CONVENTIONS.md (Figures).
"""

from collections.abc import Sequence
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
# A QC band's colour, from the Okabe-Ito data colours (ctx.qc.registry, conf/qc_metrics.toml).
BAND_COLOURS = {"good": "#009E73", "usable": "#E69F00", "bad": "#D55E00"}
# Cividis (Nuñez et al. 2018) for continuous values, in both themes: matplotlib 3.11.2's
# colormap at the midpoints of eight equal bands; its ends are #00224E and #FEE838. The stops
# run from 13.1:1 to 1.4:1 against white and back against black, so an isoline's label sits on
# a box of the theme's grid grey (tests/test_viz.py holds the stops to the colormap).
CIVIDIS = ("#002E6A", "#32436D", "#535A6D", "#6F7073", "#8C8878", "#ABA072", "#CDBB63",
           "#F0D846")  # fmt: skip
# Opacity of a band shaded behind points: a tint, so the points drawn over it keep their colour.
BAND_ALPHA = 0.15

# Liberation Sans is Arial's metric twin, which a machine without Arial may redistribute.
FONT_CHAIN = ("Arial", "Liberation Sans", "DejaVu Sans")
SIZES = {
    "base": 9.0, "axis_label": 10.0, "tick": 9.0, "legend": 9.0, "annotation": 9.0,
    "caption": 9.0, "panel": 12.0,
}  # fmt: skip
# No text below this; a mathtext superscript on a 9 pt tick would be 6.3 pt, so none is drawn.
MIN_TEXT_PT = 8.0
# A shaded band has no outline, as a point has no edge.
# "sample": one line per sample, three times a plot line so each stays traceable (Rytis
# 2026-09-30, the base-quality curves).
LINE_WIDTHS = {"axes": 1.0, "line": 1.0, "box": 1.5, "errorbar": 1.5, "marker_edge": 0.0,
               "band_edge": 0.0, "sample": 3.0}  # fmt: skip
TICK_PT = {"major": 3.0, "minor": 1.5}
# Point diameter in pt, and its opacity, so overplotted points still show their density.
MARKER_SIZE = 6.0
# A bubble's diameter in pt at a size value of 0 and at the largest value; area grows linearly.
BUBBLE_PT = {"min": 3.0, "max": 14.0}
# Size values the legend shows, as shares of the largest.
BUBBLE_KEYS = (0.25, 0.5, 1.0)
POINT_ALPHA = 0.5
PX_PER_PT = 96 / 72
# A plot on the page draws its text, points and marks 25 % larger than the figure file does
# (Rytis 2026-09-30); the PDF and SVG keep the journal sizes above.
PAGE_SCALE = 1.25
# Plotly marker diameter in pixels: MARKER_SIZE on screen, at the page's scale.
PLOTLY_MARKER_PX = MARKER_SIZE * PX_PER_PT * PAGE_SCALE
# Panel letters: bold capitals at the top left, nudged up by dy_pt (no figure has panels yet).
PANEL_LABEL = {"size": SIZES["panel"], "weight": "bold", "dx_pt": 0.0, "dy_pt": 2.0}
WIDTHS_MM = {"single": 85.0, "double": 170.0}
# Every figure takes the double column (Rytis, 2026-10-01): a height is the single column's
# 85 mm times the type's aspect, so at 85 mm wide every figure stood taller than wide.
FIGURE_WIDTH = "double"
ASPECT = {"box": 1.1, "scatter": 1.1, "bar": 1.0, "strip": 1.1, "flow": 1.1, "contour": 1.1}
# Point shapes for a second category, in level order: (matplotlib, plotly).
MARKERS = (
    ("o", "circle"),
    ("^", "triangle-up"),
    ("s", "square"),
    ("D", "diamond"),
    ("v", "triangle-down"),
)
# A heatmap's cells: the three opacities of a stepped diverging scale, weakest first; the row
# pitch that holds a 9 pt label (as BAR_ROW_MM); room for axes, tracks, legend and caption;
# and the most rows or columns that are labelled, above which hover and the CSV name them.
HEAT_ALPHAS = (0.33, 0.67, 1.0)
HEAT = {"row_mm": 4.5, "margin_mm": 70.0, "label_rows": 60, "label_columns": 40,
        "max_rows_drawn": 60}  # fmt: skip
# A flow diagram: node width and the share of the width kept for the last stage's labels,
# both as fractions of the axes; the gap between nodes is one label's pitch.
# label_share: the most of a flow's height its labels and gaps may take; the rest is bars.
FLOW = {"node": 0.02, "label_room": 0.3, "pad": 0.01, "label_share": 0.6}
PLOTLY_FLOW = {"node_px": 16, "pad_px": 18}
# A density's outer highest-density region is dashed, its inner one solid.
DENSITY_DASH: dict[str, tuple[Literal["-", "--"], Literal["solid", "dash"]]] = {
    "inner": ("-", "solid"), "outer": ("--", "dash"),
}  # fmt: skip
# Space a label box leaves around its text, in font-size units.
LABEL_BOX_PAD = 0.15
# A strip plot spreads a category's points over this share of its slot either side of the
# centre, and draws its median ± k MAD bar over a little more.
STRIP = {"jitter": 0.3, "range": 0.4}
# Average glyph width in em, used to wrap axis titles and labels to the room they have.
AVERAGE_CHAR_EM = 0.55
# The caption's inset from each side, as a share of the figure's width; annotate.place_caption
# draws at it and wraps inside it.
CAPTION_INSET = 0.02
# Empty space above the tallest box for the p/r labels, as a fraction of the data range.
HEADROOM = 0.35
LINE_SPACING = 1.3
# Legend geometry in font-size units; marker handles need less than matplotlib's line default.
# handleheight is matplotlib's own 0.7, set here so the page's swatch reads the same value.
LEGEND = {"handlelength": 1.0, "handleheight": 0.7, "columnspacing": 1.0, "handletextpad": 0.4}
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
# Page furniture per theme, from the ALI report's design (Rytis, 2026-09-29, D-024): a
# near-white page, white cards on it, an inset grey, hairlines. Never a data colour for a link
# or a pill, so neither reads as a series; ALI's blue accent is #0072B2, a data colour, so the
# link keeps its own blue, and ALI's faint grey (3.00:1 on the inset) is not used for text.
# tests/test_report.py::test_text_contrast measures every text colour on every surface.
CHROME = {
    "white": {
        "page": "#FCFCFD", "surface": "#FFFFFF", "surface-2": "#F6F7F9",
        "text": "#181B1F", "muted": "#5B6168", "line": "#E4E7EB", "line-strong": "#D2D7DD",
        "link": "#1D4E89", "warn": "#B45309", "warn-fill": "#E69F00", "good-fill": "#009E73",
        "pass-fg": "#1B5E20", "pass-bg": "#E6F2E7", "fail-fg": "#8E1B1B", "fail-bg": "#FBE9E7",
        "recorded-fg": "#4D4D4D", "recorded-bg": "#EDEDED",
        "made-fg": "#7A4B00", "made-bg": "#FBEFD9",
        "shadow": "0 1px 2px rgba(16,24,40,.04), 0 8px 24px -12px rgba(16,24,40,.18)",
    },
    "black": {
        "page": "#101316", "surface": "#1A1E22", "surface-2": "#23282E",
        "text": "#E7E9EC", "muted": "#9BA2AA", "line": "#2A2F36", "line-strong": "#3A4047",
        "link": "#8CB8F0", "warn": "#F0B45E", "warn-fill": "#E69F00", "good-fill": "#34D399",
        "pass-fg": "#9BD3A5", "pass-bg": "#16301C", "fail-fg": "#F4ABA3", "fail-bg": "#3D1816",
        "recorded-fg": "#C8C8C8", "recorded-bg": "#262626",
        "made-fg": "#F0C27A", "made-bg": "#3A2A10",
        "shadow": "0 1px 2px rgba(0,0,0,.4), 0 14px 34px -18px rgba(0,0,0,.7)",
    },
}  # fmt: skip
# Text colour, and the surfaces it is read on, for the contrast test; tints are mixed in.
TEXT_ON = {
    "text": ("page", "surface", "surface-2"), "muted": ("page", "surface", "surface-2"),
    "link": ("page", "surface", "surface-2"), "warn": ("page", "surface-2"),
    "pass-fg": ("pass-bg",), "fail-fg": ("fail-bg",), "recorded-fg": ("recorded-bg",),
    "made-fg": ("made-bg",),
}  # fmt: skip
# How much of a fill colour tints a callout, and its border, over the page (ALI's color-mix).
TINT = {"fill": "9%", "border": "45%"}
# The page text's own chain, ALI's; Arial, first in the figures' chain, has Helvetica's widths.
PAGE_FONT_CHAIN = ("Helvetica Neue", "Helvetica", "Arial", "Liberation Sans")
MONO_CHAIN = ("ui-monospace", "SF Mono", "Menlo", "Consolas", "Liberation Mono", "monospace")
# The plots on a page run 3 pt larger than the figure sizes (Rytis, 2026-09-28). The PDF and
# SVG files keep SIZES, which their layout checks are measured at.
PAGE_BOOST_PT = 3.0
# Page text in px, set on its own: it no longer follows the figure sizes.
PAGE_BASE_PX = 18
TYPE_PX = {"xs": 16, "sm": 17, "md": 20, "lg": 24, "xl": 32, "display": 48}


def on_page(size: float) -> float:
    """A figure text size in pt as a plot on the page shows it, in px: 3 pt larger, then scaled."""
    return (size + PAGE_BOOST_PT) * PX_PER_PT * PAGE_SCALE


SPACE_PX = (4, 8, 12, 16, 24, 32, 48, 64, 96)
# ALI's page geometry: a 1120 px column with 24 px gutters, 14 px card and 9 px inset radii,
# a sticky top bar blurred behind; the top bar's height offsets every sticky or anchored thing.
LAYOUT = {
    "content": "1072px", "rail": "208px", "gap": "56px", "gutter": "16px",
    "gutter-wide": "24px", "key-min": "160px", "card-min": "240px", "reading": "320px",
    "radius": "14px", "radius-sm": "9px", "radius-xs": "6px",
    "rule": "1px", "rule-mid": "2px", "rule-accent": "3px", "rule-strong": "4px",
    "topbar": "60px", "topbar-blur": "8px", "dot": "11px",
    # A glossary icon's target, WCAG 2.2's minimum target size (2.5.8), and its tooltip's width.
    "hit": "24px", "tip": "352px",
}  # fmt: skip
# Media queries cannot read var(), so the breakpoints are injected into the stylesheet.
RAIL_FROM_PX, NARROW_PX = 1100, 640
PLOTLY_HEIGHT_PX = {"scatter": 460, "box": 440, "bar": 400, "strip": 440, "flow": 460,
                    "contour": 460, "linked": 820}  # fmt: skip
# Horizontal bars grow per category, as the static figure does with BAR_ROW_MM.
PLOTLY_BAR_BASE_PX, PLOTLY_BAR_ROW_PX = 150, 28
PLOTLY_MARGIN_PX = {"l": 8, "r": 8, "t": 48, "b": 8}  # automargin grows each side as needed
PLOTLY_TICK_PX, PLOTLY_TITLE_STANDOFF_PX = 4, 8
# A legend entry's line on the page: plotly.js adds 3 px to each entry's text line
# (legend/draw.js, Math.max(h, 16) + 3; 29 px at a 20 px font, measured 2026-10-01), draws
# its symbol in 40 px before the text (itemwidth 30 and an item gap of 5 each side), and pads
# the legend by that gap above and below.
PLOTLY_LEGEND_LINE_PX = on_page(SIZES["legend"]) * LINE_SPACING + 3
PLOTLY_LEGEND_SYMBOL_PX, PLOTLY_LEGEND_PAD_PX = 40, 10
# A legend swatch on the page is the static's patch height at the legend's page size:
# plotly's default 6 px square beside a 20 px label hid the palest levels (2026-10-01).
PLOTLY_LEGEND_SWATCH_PX = LEGEND["handleheight"] * on_page(SIZES["legend"])
# The narrowest page a figure is checked at is a 375 px phone. Its plot is that less the
# gutters, the card's padding (--space-3 below NARROW_PX) and its border on each side: 317 px,
# as measured in Chromium (2026-10-01). Text a plot draws beside its cells wraps to fit it.
# The stylesheet sets those three again; tests/test_viz_heatmap.py reads them against this.
PHONE_PX = 375
PLOTLY_NARROW_PX = PHONE_PX - 2 * (int(LAYOUT["gutter"].removesuffix("px")) + SPACE_PX[2]
                                   + int(LAYOUT["rule"].removesuffix("px")))  # fmt: skip
# A hover line on that phone: plotly pads the text 3 px inside its box and draws a 6 px arrow
# (fx/constants.js HOVERTEXTPAD, HOVERARROWSIZE); a wider box ran off the plot (2026-10-01).
PLOTLY_HOVER_PX = PLOTLY_NARROW_PX - 2 * (3 + 6)
# A heatmap row on the page holds its label as HEAT's row_mm holds a 9 pt one in the file:
# the same pitch per point of label, at the label's page size. 18 px held a 20 px label
# only on every other row, and plotly dropped the rest (2026-10-01).
PLOTLY_HEAT_ROW_PX = HEAT["row_mm"] / MM_PER_INCH * 72 / SIZES["tick"] * on_page(SIZES["tick"])
# Below a heatmap's column labels and title: the bottom margin, the title's standoff and the
# tick labels' pad. Above it the legend takes the top margin, or the margin the legend.
PLOTLY_HEAT_BASE_PX = PLOTLY_MARGIN_PX["b"] + PLOTLY_TITLE_STANDOFF_PX + PLOTLY_TICK_PX
TRANSPARENT = "rgba(0,0,0,0)"
# Dash of a threshold line, in both renderers.
THRESHOLD_DASH: tuple[Literal["--"], Literal["dash"]] = ("--", "dash")
# One line per sample, each a unique colour and dash pair. The colours are the data colours with
# at least 3:1 contrast on both theme backgrounds, measured 2026-09-30 (white, black): #009E73
# 3.42, 5.52; #0072B2 5.19, 3.64; #D55E00 3.87, 4.88; #CC79A7 3.06, 6.17. The yellow measured
# 1.32 on white. Dashes as (matplotlib, plotly) names; all colours solid first, then dashed.
LINE_COLOURS = ("#009E73", "#0072B2", "#D55E00", "#CC79A7")
MplDash = Literal["solid", "dashed", "dotted", "dashdot"]
DASHES: tuple[tuple[MplDash, str], ...] = (
    ("solid", "solid"), ("dashed", "dash"), ("dotted", "dot"), ("dashdot", "dashdot"),
)  # fmt: skip
MAX_DISTINCT = len(LINE_COLOURS) * len(DASHES)


def distinct_styles(levels: Sequence[str]) -> dict[str, tuple[str, tuple[MplDash, str]]]:
    """A unique (colour, (matplotlib dash, plotly dash)) per level, in the order given.

    Raises:
        ValueError: If there are more levels than colour and dash pairs.
    """
    if len(levels) > MAX_DISTINCT:
        raise ValueError(f"{len(levels)} lines; at most {MAX_DISTINCT} distinct styles")
    n = len(LINE_COLOURS)
    return {level: (LINE_COLOURS[i % n], DASHES[i // n]) for i, level in enumerate(levels)}


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


def tint(colour: str, alpha: float) -> str:
    """A hex colour as CSS rgba at an opacity: a filled area's tint, its outline left whole."""
    red, green, blue = (int(colour[i : i + 2], 16) for i in (1, 3, 5))
    return f"rgba({red},{green},{blue},{alpha})"


def _translucent(targets: tuple[str, str], grounds: tuple[str, str]) -> tuple[str, float]:
    """The one colour and opacity that, laid over each ground, land nearest its target.

    Over ground g at opacity a, colour c shows as a·c + (1 − a)·g. The two grounds share a
    and c, so 1 − a is the least-squares slope of the targets' difference on the grounds'
    difference; c then splits each channel's error evenly between them.

    Returns:
        The colour as hex, and its opacity.
    """
    (t1, t2), (g1, g2) = ([[int(c[i : i + 2], 16) for i in (1, 3, 5)] for c in pair]
                          for pair in (targets, grounds))  # fmt: skip
    channels = list(zip(t1, t2, g1, g2, strict=True))
    keep = sum((a - b) * (c - d) for a, b, c, d in channels) / sum(
        (c - d) ** 2 for _, _, c, d in channels
    )
    colour = [round((a + b - keep * (c + d)) / 2 / (1 - keep)) for a, b, c, d in channels]
    return "#" + "".join(f"{min(max(v, 0), 255):02X}" for v in colour), round(1 - keep, 3)


# A heatmap's "nothing to show" level on the page. A plotly trace keeps its colour when the
# theme flips (as the reference grey does), so it is the one translucent colour that, over
# each theme's figure card, lands nearest that theme's grid: within 3 of 255 per channel on
# both. The white grid had been the brightest cell on a dark page (2026-10-01).
PLOTLY_NO_LEVEL = _translucent(
    (THEMES["white"].grid, THEMES["black"].grid),
    (CHROME["white"]["surface"], CHROME["black"]["surface"]),
)


def plotly_steps(colours: list[str]) -> list[list[Any]]:
    """A plotly colorscale over 0 to n - 1 in which integer k, and only k, takes colour k."""
    n = len(colours)
    scale: list[list[Any]] = []
    for k, colour in enumerate(colours):
        scale += [[k / n, colour], [(k + 1) / n, colour]]
    return scale


def figure_inches(width: str, kind: str) -> tuple[float, float]:
    """Figure size in inches for a column width and a figure type's aspect ratio."""
    inches = WIDTHS_MM[width] / MM_PER_INCH
    return inches, WIDTHS_MM["single"] / MM_PER_INCH * ASPECT[kind]


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
            # Every tick carries its full value: no detached multiplier or offset above the
            # axis, which on Step-1's rRNA strip (values near 1e-5) landed on the legend.
            "axes.formatter.useoffset": False,
            "axes.formatter.limits": (-9, 9),
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
            # Layout measured as the PDF and SVG set text: hinted glyphs measured a caption line
            # 238.3 pt that both files set at 227.8 pt (2026-10-01), and only vector files ship.
            "text.hinting": "no_hinting",
            "legend.frameon": False,
            "legend.handlelength": LEGEND["handlelength"],
            "legend.handleheight": LEGEND["handleheight"],
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
        "--font-page": ", ".join([*(f'"{f}"' for f in PAGE_FONT_CHAIN), "sans-serif"]),
        "--size-base": f"{PAGE_BASE_PX}px",
    }
    return "; ".join(f"{k}: {v}" for k, v in pairs.items()) + ";"


def css_layout() -> str:
    """Type scale, spacing and layout as CSS custom properties; the same in both themes."""
    pairs = {
        **{f"--fs-{name}": f"{px}px" for name, px in TYPE_PX.items()},
        **{f"--space-{i + 1}": f"{px}px" for i, px in enumerate(SPACE_PX)},
        **{f"--{name}": value for name, value in LAYOUT.items()},
        **{f"--tint-{name}": value for name, value in TINT.items()},
        "--font-mono": ", ".join(f'"{f}"' if " " in f else f for f in MONO_CHAIN),
    }
    return "; ".join(f"{k}: {v}" for k, v in pairs.items()) + ";"
