"""How wide a text is set: the one measure of text width, from Liberation Sans, Arial's twin.

The figure files' captions, the page plots' wrapped labels and a page's own PDF and SVG of a
changed view (D-033) all wrap with ``wrap_to``, so the view's caption breaks where the files'
caption breaks: the page reads the same numbers (``table``). A width is the sum of the glyphs'
advances and of the font's kerning pairs, its legacy ``kern`` table, which FreeType applies
when matplotlib sets the files' text. Arial's table differs from it only in pairs with a
no-break space, a soft hyphen or a semicolon after T, V, W or Y. On the 99 captions of the
course's figures (2026-10-07) matplotlib set every line 0.0004 to 0.038 pt narrower than this
measure at 170 mm, and the breaks matched at 170 and 85 mm; without the kerning the measure
ran up to 4.33 pt wide and broke 5 captions a word early at 85 mm. A page's PDF embeds these
fonts as Ctx Sans, its subset renamed; Chromium drew the band map's run names and the parallel
sets' titles within 0.5 px of the advances (2026-10-01).

A vertical axis title and a flow's labels wrap with it too, from ``words``, which holds a
parenthesised span such as a unit whole, as the page's ``ctxExport.words`` does.
"""

import io
from functools import cache
from typing import Any

from ctx.report.export import FONTS, VENDOR, vendored


@cache
def register_fonts() -> None:
    """Give matplotlib the vendored faces, so the files' text is set in the face measured here.

    The files name Liberation Sans first (``style.FONT_CHAIN``), and only the vendored files
    answer to that name: an installed copy of another release would set other bytes, and
    matplotlib takes the first of two equal matches, which is the installed one. Before the
    faces were registered, a machine with neither Arial nor Liberation Sans fell back to DejaVu
    Sans, wider than this measure, and a caption wrapped by it would cross its inset.
    """
    from matplotlib import font_manager  # noqa: PLC0415 - heavy

    manager = font_manager.fontManager
    manager.ttflist = [f for f in manager.ttflist if f.name != "Liberation Sans"]
    for name in FONTS.values():
        vendored(name)  # refuses a file that is not the bytes its manifest names
        font_manager.fontManager.addfont(str(VENDOR / name))


@cache
def _metrics(style: str) -> tuple[dict[str, int], int, dict[str, int], int]:
    """The vendored face ``style``'s numbers a width is measured from.

    Its advance per character, the advance of a character it lacks, its kerning per pair of
    characters and its units per em.
    """
    from fontTools.ttLib import TTFont  # type: ignore[import-untyped]  # noqa: PLC0415 - heavy

    font = TTFont(io.BytesIO(vendored(FONTS[style])))
    metrics = font["hmtx"].metrics
    glyphs = font.getBestCmap()
    advances = {chr(code): metrics[glyph][0] for code, glyph in glyphs.items()}
    by_glyph: dict[str, list[str]] = {}
    for code, glyph in glyphs.items():
        by_glyph.setdefault(glyph, []).append(chr(code))
    kern = {
        left + right: value
        for table in font["kern"].kernTables
        for (first, second), value in table.kernTable.items()
        for left in by_glyph.get(first, ())
        for right in by_glyph.get(second, ())
    }
    return advances, metrics[".notdef"][0], kern, font["head"].unitsPerEm


def text_width(text: str, size: float, bold: bool = False) -> float:
    """The width of one line of ``text`` set at ``size``, kerned, in the unit of ``size``."""
    advances, missing, kern, per_em = _metrics("bold" if bold else "normal")
    total = sum(advances.get(c, missing) for c in text)
    total += sum(kern.get(text[i : i + 2], 0) for i in range(len(text) - 1))
    return size * total / per_em


def words(text: str) -> list[str]:
    """``text``'s words, a parenthesised span held as one, as ``ctxExport.words`` splits a title.

    A unit split at its spaces, ``(log2(count + 1))``, left ``1))`` alone on a line (M3.6), and
    a flow's level ``Reads > 76 nt (max 500 nt)`` ended a line at ``(max`` (2026-10-07).
    """
    held: list[str] = []
    word, depth = "", 0
    for char in text:
        if char.isspace() and depth <= 0:
            held.append(word)
            word = ""
        else:
            depth += (char == "(") - (char == ")")
            word += char
    return [w for w in (*held, word) if w]


def wrap_to(text: str | list[str], room: float, size: float) -> list[str]:
    """Words greedily onto lines no wider than ``room`` at ``size``; a word is never split.

    A line ends only where its next word would cross ``room``; a word wider than ``room``
    takes a line of its own. A text is split at its white space, any run of it one space;
    a list is its words already, as ``words`` holds a unit, and as ``ctxExport.wrap`` takes
    either.
    """
    lines: list[str] = []
    for word in text.split() if isinstance(text, str) else text:
        trial = f"{lines[-1]} {word}" if lines else word
        if lines and text_width(trial, size) <= room:
            lines[-1] = trial
        else:
            lines.append(word)
    return lines


def table(characters: frozenset[str]) -> dict[str, Any]:
    """The regular face's numbers for these characters, as a page's export measures a line.

    ``ctxExport.width`` adds them as ``text_width`` does, in the same order of operations.
    """
    advances, missing, kern, per_em = _metrics("normal")
    return {
        "per_em": per_em,
        "missing": missing,
        "advances": {c: advances[c] for c in sorted(characters) if c in advances},
        "kern": {pair: v for pair, v in sorted(kern.items()) if set(pair) <= characters},
    }
