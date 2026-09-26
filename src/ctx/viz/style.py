"""The only place a colour, font, size, line width or figure dimension is set.

Order and font follow CONVENTIONS.md (Figures); line widths, figure widths, black as ink
and the seven-category cap follow the global figure rules where the project is silent.
"""

from dataclasses import dataclass
from typing import Any

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
WIDTHS_MM = {"single": 85.0, "double": 170.0}
ASPECT = {"box": 1.1}
# Average glyph width in em, used to wrap the caption to the width a figure really has.
AVERAGE_CHAR_EM = 0.55
# Above this many features a box plot takes the double column: at 14 pt the labels of more
# features do not fit in 85 mm (measured 2026-09-26, 8 genes: every label overlapped).
SINGLE_COLUMN_MAX_FEATURES = 3
# Empty space above the tallest box for the p/r labels, as a fraction of the data range.
HEADROOM = 0.35
LINE_SPACING = 1.3


@dataclass(frozen=True)
class Theme:
    """Colours of one theme; the name is the file token (``_white`` / ``_black``)."""

    name: str
    ink: str
    background: str
    grid: str
    reference: str


THEMES = {
    "white": Theme(
        "white", ink=INK_COLOUR, background="#FFFFFF", grid="#E5E5E5", reference="#999999"
    ),
    "black": Theme(
        "black", ink="#F2F2F2", background="#111111", grid="#333333", reference="#8C8C8C"
    ),
}


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
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def plotly_templates() -> dict[str, Any]:
    """Plotly layout templates for the two themes, built on ``plotly_white``."""
    import plotly.io as pio  # noqa: PLC0415 - heavy; imported after the runtime profile

    templates = {}
    for name, theme in THEMES.items():
        template = pio.templates["plotly_white"].to_plotly_json()
        template["layout"].update(
            {
                "font": {
                    "family": ", ".join([*FONT_CHAIN, "sans-serif"]),
                    "size": SIZES["base"],
                    "color": theme.ink,
                },
                "colorway": list(DATA_COLOURS),
                "paper_bgcolor": theme.background,
                "plot_bgcolor": theme.background,
            }
        )
        for axis in ("xaxis", "yaxis"):
            template["layout"][axis].update(
                {"gridcolor": theme.grid, "linecolor": theme.ink, "zerolinecolor": theme.grid}
            )
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
        "--font": ", ".join([*(f'"{f}"' for f in FONT_CHAIN), "sans-serif"]),
        "--size-base": f"{SIZES['base']}px",
    }
    return "; ".join(f"{k}: {v}" for k, v in pairs.items()) + ";"
