"""Inline markup in page text: backticked code, and glossary marks ``[[term]]`` (D-048).

One reading of both, ``pieces()``, used by the prose lint, the glossary lint, the page renderer
and every door out of HTML: a PDF or SVG and ``--help`` print a mark as its display text.
"""

import re
from collections.abc import Callable
from typing import NamedTuple

from markupsafe import Markup, escape

# [[key]] shows the key as written, [[key|display]] shows the display text; a "!" after the key
# gives this mention the icon even when its section has shown the term already. The key is an
# id, a term or an alias of a glossary entry, in any case.
MARK = re.compile(r"\[\[(?P<key>[^\[\]|!]+?)(?P<force>!?)(?:\|(?P<display>[^\[\]|]+))?\]\]")
_CODE = re.compile(r"`([^`]+)`")


class Mark(NamedTuple):
    """One glossary mark: the name it resolves by, what the page shows, and a forced icon."""

    key: str
    display: str
    force: bool


class Piece(NamedTuple):
    """A stretch of page text as written: words, a backticked code span, or a glossary mark."""

    raw: str
    code: bool = False
    mark: Mark | None = None


def _mark(match: re.Match[str]) -> Mark:
    """The mark a match of ``MARK`` holds."""
    key = match["key"].strip()
    return Mark(key, (match["display"] or key).strip(), bool(match["force"]))


def _words(stretch: str) -> list[Piece]:
    """A stretch outside code as words and marks, in order."""
    out, at = [], 0
    for match in MARK.finditer(stretch):
        out += [Piece(stretch[at : match.start()]), Piece(match.group(0), mark=_mark(match))]
        at = match.end()
    return [p for p in [*out, Piece(stretch[at:])] if p.raw]


def pieces(text: str) -> list[Piece]:
    """The text in order as words, code and marks.

    A mark inside backticks is code, not a mark: the page prints it literally, so nothing may
    read it as one, neither the lint nor the tooltips a page embeds.
    """
    out, at = [], 0
    for code in _CODE.finditer(text):
        out += [*_words(text[at : code.start()]), Piece(code.group(0), code=True)]
        at = code.end()
    return out + _words(text[at:])


def marks(text: str) -> list[Mark]:
    """Every glossary mark in a text, in order."""
    return [p.mark for p in pieces(text) if p.mark is not None]


def plain(text: str) -> str:
    """The text as a reader sees it outside a page: each mark replaced by its display text."""
    return "".join(p.raw if p.mark is None else p.mark.display for p in pieces(text))


def html(text: str, term: Callable[[Mark], str | None] | None = None) -> Markup:
    """Page text as HTML: escaped, each backticked span as code, each mark rendered by ``term``.

    Args:
        text: The page text, with its marks.
        term: Returns a mark's HTML, or None to show its display text alone; without it
            every mark shows its display text.

    Returns:
        Safe HTML.
    """
    out = []
    for piece in pieces(text):
        if piece.code:
            out.append(f"<code>{escape(piece.raw[1:-1])}</code>")
        elif piece.mark is None:
            out.append(str(escape(piece.raw)))
        else:
            shown = term(piece.mark) if term else None
            out.append(shown if shown is not None else str(escape(piece.mark.display)))
    return Markup("".join(out))
