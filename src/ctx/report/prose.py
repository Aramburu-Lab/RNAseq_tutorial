"""The prose rules of CONVENTIONS.md, applied to the page text in ``content/*.toml``.

A page file carries the prose sections of the nine-label skeleton. Each section has a word
budget; every sentence is at most 25 words and every paragraph at most three sentences; no
banned word appears, a figurative verb in none of its forms; a number in the prose comes from a
``{placeholder}`` or sits in ``code``, never typed as a literal; and no text but the Questions
asks one: a question mark in a title, a heading or a sentence is a violation. A glossary mark
counts as its display text (D-048). Run with ``python -m ctx.report.prose``; exit 1 lists
every violation.
"""

import re
import sys
import tomllib
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

from ctx.report.marks import plain

CONTENT = Path(__file__).parent / "content"
# Words per section, per item for lists (CONVENTIONS.md, Report prose).
BUDGETS = {"purpose": 30, "input": 40, "method": 90, "interpretation": 120}
ITEMS = {"questions": (3, 20), "pitfalls": (3, 20)}
CAPTION_WORDS = 40
CARD_WORDS = 12
SENTENCE_WORDS = 25
PARAGRAPH_SENTENCES = 3
BANNED = (
    "we", "let's", "you will", "in this module", "crucial", "powerful", "robust", "seamless",
    "simply", "actually", "very", "unlock", "dive", "journey", "landscape", "leverage",
    "delve", "awesome", "exciting",
    # Software and data do not know, bring, want, tell or arrive, in any form (D-048).
    "know", "knows", "knew", "known", "knowing", "bring", "brings", "brought", "bringing",
    "want", "wants", "wanted", "wanting", "tell", "tells", "told", "telling",
    "arrive", "arrives", "arrived", "arriving",
)  # fmt: skip
_BANNED = re.compile(r"\b(" + "|".join(re.escape(b) for b in BANNED) + r")\b", re.I)
# A placeholder, or text in backticks, is not prose: its digits come from code.
_NOT_PROSE = re.compile(r"\{[^{}]*\}|`[^`]*`")
# A sentence may start in lower case (fastp, lncRNAs, nf-core), so a stop before a lower-case
# word ends one too, unless it closes an abbreviation that runs on (e.g., i.e., vs., cf., et al.).
_SENTENCE = re.compile(
    r"(?<=[.!?])\s+(?=[A-Z`{])"
    r"|(?<![Ee]\.g\.)(?<!\bi\.e\.)(?<!\bvs\.)(?<!\bal\.)(?<!\bcf\.)(?<=[.!?])\s+(?=[a-z])"
)
# A quantity standing alone; digits inside a name (Step-1, M5, SHA-256, chr22) are not one,
# nor a strand end (5', 3'). A trailing full stop or comma ends the sentence, not the number.
_NUMBER = re.compile(r"(?<![\w-])\d+(?:[,.]\d+)*(?![\w'-])")


def words(text: str) -> int:
    """Words in text, a placeholder counted as one."""
    return len(_NOT_PROSE.sub("x", text).split())


def _sentences(paragraph: str) -> list[str]:
    """The sentences of one paragraph, split after a full stop, question or exclamation mark."""
    return [s for s in _SENTENCE.split(paragraph.strip()) if s]


def outside_code(text: str) -> str:
    """The text with each placeholder and backticked span blanked: what is left is prose."""
    return _NOT_PROSE.sub(" ", text)


def quantities(text: str) -> list[str]:
    """The numbers typed in a text, outside placeholders and code, a mark read as displayed."""
    return _NUMBER.findall(outside_code(plain(text)))


def check_text(where: str, text: str, budget: int | None, *, numbers: bool = True) -> Iterator[str]:
    """Every violation in one piece of prose, each naming where it is.

    Args:
        where: The file and key, printed with each violation.
        text: The prose; a glossary mark counts as its display text.
        budget: Most words allowed, or None.
        numbers: Whether a typed number is a violation; a glossary definition states numbers
            and gives their source instead (D-048).
    """
    text = plain(text)
    if budget is not None and (n := words(text)) > budget:
        yield f"{where}: {n} words; at most {budget}"
    for paragraph in text.split("\n\n"):
        sentences = _sentences(paragraph)
        if len(sentences) > PARAGRAPH_SENTENCES:
            n = len(sentences)
            yield f"{where}: a paragraph of {n} sentences; at most {PARAGRAPH_SENTENCES}"
        for sentence in sentences:
            if (n := words(sentence)) > SENTENCE_WORDS:
                yield f"{where}: a sentence of {n} words; at most {SENTENCE_WORDS}: {sentence[:60]}"
    prose = _NOT_PROSE.sub(" ", text)
    for match in _BANNED.finditer(prose):
        yield f"{where}: banned word {match.group(0)!r}"
    if numbers and quantities(text):
        yield f"{where}: a typed number; take it from a {{placeholder}} or put it in `code`"


def asks(where: str, text: str) -> Iterator[str]:
    """A question mark outside placeholders and code: only the Questions section asks."""
    if "?" in outside_code(plain(text)):
        yield f"{where}: a question; state it, or ask it under Questions"


def check_asks(name: str, node: Any) -> Iterator[str]:
    """Every string of a content file but the Questions: titles, headings, leads and prose."""
    if isinstance(node, str):
        yield from asks(name, node)
    elif isinstance(node, list):
        for i, item in enumerate(node):
            yield from check_asks(f"{name}[{i}]", item)
    elif isinstance(node, dict):
        for key, value in node.items():
            if key != "questions":
                yield from check_asks(f"{name}.{key}", value)


def check_page(name: str, raw: Mapping[str, Any]) -> Iterator[str]:
    """Violations in a page file: the sections, their lists, and the figure captions.

    A table under ``variants`` replaces some of the page's sections for one design
    (``text.skeleton``), so each is held to the same budgets as the sections it replaces.
    """
    for variant, sections in raw.get("variants", {}).items():
        yield from check_page(f"{name}: variants.{variant}", sections)
    for key, budget in BUDGETS.items():
        if key in raw:
            yield from check_text(f"{name}: {key}", raw[key], budget)
    for key, (most, each) in ITEMS.items():
        found = raw.get(key, [])
        # A table of variants (M5's, M6's questions) is asked by what the data hold: each
        # variant is held to the item budget, and the count to the page that asks them.
        items = list(found.values()) if isinstance(found, dict) else found
        if isinstance(found, list) and len(items) > most:
            yield f"{name}: {key}: {len(items)} items; at most {most}"
        for i, item in enumerate(items):
            text = item if isinstance(item, str) else item.get("prompt", "")
            yield from check_text(f"{name}: {key}[{i}]", text, each)
    for fig, entry in raw.get("figures", {}).items():
        if isinstance(entry, dict) and "caption" in entry:
            where = f"{name}: figures.{fig}.caption"
            yield from check_text(where, entry["caption"], CAPTION_WORDS)
            # The PDF prints a caption as it is, so a backtick would show there.
            if "`" in entry["caption"]:
                yield f"{where}: a backtick; the PDF prints it, so name the literal in words"


def check_cards(raw: Mapping[str, Any]) -> Iterator[str]:
    """Parameter cards: meaning and effect, each at most twelve words, per-command effects too."""
    for flag, card in raw.get("options", {}).items():
        for part in ("meaning", "effect"):
            if part in card:
                yield from check_text(f"cards: {flag}.{part}", card[part], CARD_WORDS)
    for command, effects in raw.get("effects", {}).items():
        for flag, effect in effects.items():
            yield from check_text(f"cards: effects.{command}.{flag}", effect, CARD_WORDS)


def check_all(name: str, node: Any) -> Iterator[str]:
    """Every string in a file with no skeleton, such as the course page and the kit README."""
    if isinstance(node, str):
        yield from check_text(name, node, None)
    elif isinstance(node, dict):
        for key, value in node.items():
            yield from check_all(f"{name}.{key}", value)


def lint(content: Path = CONTENT) -> list[str]:
    """Every violation in the page files that follow the skeleton, the cards and the course."""
    problems: list[str] = []
    for path in sorted(content.glob("*.toml")):
        raw = tomllib.loads(path.read_text())
        if "purpose" in raw:
            problems += check_page(path.stem, raw)
        if path.stem == "cards":
            problems += check_cards(raw)
        if path.stem == "course":
            problems += check_all("course", raw)
        problems += check_asks(path.stem, raw)
    return problems


if __name__ == "__main__":
    found = lint()
    print("\n".join(found) or "prose lint: no violations")
    sys.exit(1 if found else 0)
