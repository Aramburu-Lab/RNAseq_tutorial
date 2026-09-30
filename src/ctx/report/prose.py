"""The prose rules of CONVENTIONS.md, applied to the page text in ``content/*.toml``.

A page file carries the prose sections of the nine-label skeleton. Each section has a word
budget; every sentence is at most 25 words and every paragraph at most three sentences; no
banned word appears; and a number in the prose comes from a ``{placeholder}`` or sits in
``code``, never typed as a literal. Run with ``python -m ctx.report.prose``; exit 1 lists
every violation.
"""

import re
import sys
import tomllib
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

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
    "delve", "awesome", "exciting", "arrive", "arrives", "bring", "brings", "knows", "wants",
    "tells",
)  # fmt: skip
_BANNED = re.compile(r"\b(" + "|".join(re.escape(b) for b in BANNED) + r")\b", re.I)
# A placeholder, or text in backticks, is not prose: its digits come from code.
_NOT_PROSE = re.compile(r"\{[^{}]*\}|`[^`]*`")
_SENTENCE = re.compile(r"(?<=[.!?])\s+(?=[A-Z`{])")
# A quantity standing alone; digits inside a name (Step-1, M5, SHA-256, chr22) are not one.
_NUMBER = re.compile(r"(?<![\w-])\d[\d,.]*(?![\w-])")


def words(text: str) -> int:
    """Words in text, a placeholder counted as one."""
    return len(_NOT_PROSE.sub("x", text).split())


def _sentences(paragraph: str) -> list[str]:
    """The sentences of one paragraph, split after a full stop, question or exclamation mark."""
    return [s for s in _SENTENCE.split(paragraph.strip()) if s]


def check_text(where: str, text: str, budget: int | None) -> Iterator[str]:
    """Every violation in one piece of prose, each naming where it is."""
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
    if _NUMBER.search(prose):
        yield f"{where}: a typed number; take it from a {{placeholder}} or put it in `code`"


def check_page(name: str, raw: Mapping[str, Any]) -> Iterator[str]:
    """Violations in a page file: the sections, their lists, and the figure captions."""
    for key, budget in BUDGETS.items():
        if key in raw:
            yield from check_text(f"{name}: {key}", raw[key], budget)
    for key, (most, each) in ITEMS.items():
        items = raw.get(key, [])
        if len(items) > most:
            yield f"{name}: {key}: {len(items)} items; at most {most}"
        for i, item in enumerate(items):
            text = item if isinstance(item, str) else item.get("prompt", "")
            yield from check_text(f"{name}: {key}[{i}]", text, each)
    for fig, entry in raw.get("figures", {}).items():
        if isinstance(entry, dict) and "caption" in entry:
            yield from check_text(f"{name}: figures.{fig}.caption", entry["caption"], CAPTION_WORDS)


def check_cards(raw: Mapping[str, Any]) -> Iterator[str]:
    """Parameter cards: meaning and effect, each at most twelve words, per-command effects too."""
    for flag, card in raw.get("options", {}).items():
        for part in ("meaning", "effect"):
            if part in card:
                yield from check_text(f"cards: {flag}.{part}", card[part], CARD_WORDS)
    for command, effects in raw.get("effects", {}).items():
        for flag, effect in effects.items():
            yield from check_text(f"cards: effects.{command}.{flag}", effect, CARD_WORDS)


def lint(content: Path = CONTENT) -> list[str]:
    """Every violation in the page files that follow the skeleton, and in the cards."""
    problems: list[str] = []
    for path in sorted(content.glob("*.toml")):
        raw = tomllib.loads(path.read_text())
        if "purpose" in raw:
            problems += check_page(path.stem, raw)
        if path.stem == "cards":
            problems += check_cards(raw)
    return problems


if __name__ == "__main__":
    found = lint()
    print("\n".join(found) or "prose lint: no violations")
    sys.exit(1 if found else 0)
