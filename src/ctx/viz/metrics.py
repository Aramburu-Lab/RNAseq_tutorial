"""How wide a text is set, from the advance widths of Liberation Sans, Arial's metric twin.

The page sets its plots in the Arial chain, and a page's own PDF embeds these widths as Ctx Sans,
its subset of Liberation Sans renamed (D-033), so they are the drawn widths, where
``style.AVERAGE_CHAR_EM`` is a guess: the band map's 19-character run names measured 246.75 px
at the page's 20 px, against 209 estimated. Kerning is left out: Chromium drew those names and
the parallel sets' titles within 0.5 px of these sums (2026-10-01).
"""

import io
from functools import cache

from ctx.report.export import FONTS, vendored


@cache
def _advances(style: str) -> tuple[dict[int, int], int, int]:
    """(advance per code point, the advance of a glyph the font lacks, units per em)."""
    from fontTools.ttLib import TTFont  # type: ignore[import-untyped]  # noqa: PLC0415 - heavy

    font = TTFont(io.BytesIO(vendored(FONTS[style])))
    metrics = font["hmtx"].metrics
    widths = {code: metrics[glyph][0] for code, glyph in font.getBestCmap().items()}
    return widths, metrics[".notdef"][0], font["head"].unitsPerEm


def text_width(text: str, size: float, bold: bool = False) -> float:
    """The advance width of one line of ``text`` set at ``size``, in the unit of ``size``."""
    widths, missing, per_em = _advances("bold" if bold else "normal")
    return size * sum(widths.get(ord(c), missing) for c in text) / per_em


def wrap_to(text: str, room: float, size: float) -> list[str]:
    """Words greedily onto lines no wider than ``room`` at ``size``; a word is never split."""
    lines: list[str] = []
    for word in text.split():
        trial = f"{lines[-1]} {word}" if lines else word
        if lines and text_width(trial, size) <= room:
            lines[-1] = trial
        else:
            lines.append(word)
    return lines
