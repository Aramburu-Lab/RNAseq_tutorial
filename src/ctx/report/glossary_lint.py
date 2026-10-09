"""The glossary lint (D-048): ``python -m ctx.report.glossary_lint``; exit 1 lists every error.

Errors: the schema (``glossary.structure``); a ``short`` over 20 words or a ``long`` over 60; a
definition that breaks the prose rules, all but the typed number, since a definition states
numbers and gives each its ``source``; a mark in a page file that names no entry, or sits
where the page prints it raw. At a release tag, also a ``confirm`` flag and a number with no
source. Warnings, printed without failing: those two before a release tag, and each glossary
term whose first mention in a section of a page file is bare, marked later or not at all.
"""

import os
import re
import subprocess
import sys
import tomllib
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from ctx.report import glossary, marks, prose

REPO = Path(__file__).resolve().parents[3]
# Words per field; typical and rna_vs_dna keep the sentence rules without a total.
DEFINITIONS = {"short": 20, "long": 60, "typical": None, "rna_vs_dna": None}
# Where a page renders a mark, by the section that carries its icon; anywhere else the page
# would print the brackets.
SECTIONS = {"purpose", "input", "method", "interpretation", "pitfalls", "questions"}
_ANSWER = {"prompt", "answer"}
# options."--x".meaning and effects.<command>."--x": a card's text is three keys deep.
_CARD_DEPTH = 3
# A name of three characters or fewer, or with two capitals (GO, OS, TPM, MultiQC), is matched
# as written: "go" and "os" are words, not the Gene Ontology or overall survival.
_EXACT = re.compile(r"^.{1,3}$|[A-Z].*[A-Z]")


def definitions(raw: list[dict[str, Any]]) -> Iterator[str]:
    """Word budgets and the prose rules on every definition."""
    for entry in raw:
        for field, budget in DEFINITIONS.items():
            if value := entry.get(field):
                where = f"{entry['id']}.{field}"
                yield from prose.check_text(where, str(value), budget, numbers=False)
                yield from prose.asks(where, str(value))


def unverified(raw: list[dict[str, Any]]) -> Iterator[str]:
    """Entries marked ``confirm``, and numbers in short, long or typical with no ``source``."""
    for entry in raw:
        if entry.get("confirm"):
            yield f"{entry['id']}: confirm: true; verify it and give its source"
        if entry.get("source"):
            continue
        for field in ("short", "long", "typical"):
            if found := prose.quantities(str(entry.get(field, ""))):
                yield f"{entry['id']}.{field}: the number {found[0]} has no source"
                break


def section(path: tuple[str, ...], page: bool) -> str | None:
    """The section a content string is rendered in with marks, or None where it is not.

    Args:
        path: The string's keys and indices in its file.
        page: Whether the file is a page file (it has a ``purpose``); else the cards.
    """
    head = path[0]
    if page and head in SECTIONS and (head != "questions" or path[-1] in _ANSWER):
        return head
    if page and head == "figures" and path[-1] == "caption":
        return "results"
    if not page and head in {"options", "effects"} and len(path) == _CARD_DEPTH:
        return "parameters"
    return None


def strings(node: Any, path: tuple[str, ...] = ()) -> Iterator[tuple[tuple[str, ...], str]]:
    """Every string of a parsed file, with its path."""
    if isinstance(node, str):
        yield path, node
    elif isinstance(node, list):
        for i, item in enumerate(node):
            yield from strings(item, (*path, str(i)))
    elif isinstance(node, dict):
        for key, value in node.items():
            yield from strings(value, (*path, key))


def _names() -> tuple[re.Pattern[str], re.Pattern[str]]:
    """Every term and alias, as matched as written and as matched in any case."""
    names = [n for e in glossary.load().values() for n in (e.term, *e.aliases)]
    exact = sorted({n for n in names if _EXACT.search(n)}, key=len, reverse=True)
    loose = sorted({n for n in names if not _EXACT.search(n)}, key=len, reverse=True)

    def pattern(found: list[str], flags: int) -> re.Pattern[str]:
        """One pattern for many names, each a whole word; the longest tried first."""
        alternatives = "|".join(re.escape(n) for n in found)
        return re.compile(rf"(?<![\w-])(?:{alternatives})(?![\w-])", flags)

    return pattern(exact, 0), pattern(loose, re.I)


def page_marks(content: Path = prose.CONTENT) -> tuple[list[str], list[str]]:
    """Errors (unknown or unrendered marks) and warnings (terms left unmarked) of the pages."""
    errors: list[str] = []
    groups: dict[tuple[str, str], list[str]] = {}
    for file in sorted(content.glob("*.toml")):
        raw = tomllib.loads(file.read_text(encoding="utf-8"))
        page = "purpose" in raw
        for path, value in strings(raw):
            where = f"{file.stem}: {'.'.join(path)}"
            place = section(path, page) if page or file.stem == "cards" else None
            found = marks.marks(value)
            if found and place is None:
                errors.append(f"{where}: a mark where the page prints it raw")
            for mark in found:
                try:
                    glossary.lookup(mark.key)
                except glossary.GlossaryError as error:
                    errors.append(f"{where}: {error}")
            if place is not None:
                groups.setdefault((file.stem, place), []).append(value)
    return errors, list(_unmarked(groups)) if not errors else []


# Per page file, a glossary name that means something else there (L.6). Matching by case is
# not enough for these: setup's "OS" is the operating system, written as survival's alias is.
OTHER_MEANINGS: dict[str, frozenset[str]] = {"setup": frozenset({"OS"})}


def _unmarked(groups: dict[tuple[str, str], list[str]]) -> Iterator[str]:
    """Per section, each entry whose first mention in reading order carries no mark."""
    patterns = _names()
    for (stem, place), texts in groups.items():
        found = _first_mentions(texts, patterns, OTHER_MEANINGS.get(stem, frozenset()))
        yield from (f"{stem}: {place}: {w}" for w in found)


def _first_mentions(
    texts: list[str], patterns: tuple[re.Pattern[str], ...], other: frozenset[str] = frozenset()
) -> Iterator[str]:
    """A section's entries whose first mention, in reading order, is bare.

    The icon sits on the first mark, so a bare mention before it leaves the reader's first
    meeting with the term unexplained; an entry the section never marks is unmarked. A name
    in ``other`` means something else on this page and is not a mention.
    """
    marked = glossary.used(texts)
    seen: set[str] = set()
    for piece in (p for value in texts for p in marks.pieces(value) if not p.code):
        if piece.mark is not None:
            seen.add(glossary.lookup(piece.mark.key).id)
            continue
        bare = prose.outside_code(piece.raw)
        for match in (m for pattern in patterns for m in pattern.finditer(bare)):
            if match.group(0) in other:
                continue
            entry = glossary.lookup(match.group(0)).id
            if entry not in seen:
                seen.add(entry)
                how = "marked only at a later mention" if entry in marked else "unmarked"
                yield f"{match.group(0)!r} is {entry}, {how}"


def release_tag() -> str | None:
    """The ``v*`` tag this lint runs at, or None: GitHub's ref, else ``git describe``."""
    name = os.environ.get("GITHUB_REF_NAME", "")
    # Only a v* tag is a release, as in git describe's --match below: not a "nightly" tag.
    if os.environ.get("GITHUB_REF_TYPE") == "tag" and name.startswith("v"):
        return name
    done = subprocess.run(
        ["git", "-C", str(REPO), "describe", "--exact-match", "--tags", "--match", "v*", "HEAD"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    return done.stdout.strip() or None


def lint(
    raw: Any = None, content: Path = prose.CONTENT, *, release: bool = False
) -> tuple[list[str], list[str]]:
    """Every error and every warning of the glossary and the marks in the page files.

    Args:
        raw: A parsed glossary file; the packaged one when None.
        content: The page files.
        release: At a release tag, unverified entries and unsourced numbers are errors.

    Returns:
        (errors, warnings).
    """
    raw = glossary.read() if raw is None else raw
    if errors := glossary.structure(raw):
        return errors, []
    errors = list(definitions(raw))
    pending = list(unverified(raw))
    found, warnings = page_marks(content)
    errors += found + (pending if release else [])
    return errors, ([] if release else pending) + warnings


if __name__ == "__main__":
    tag = release_tag()
    errors, warnings = lint(release=tag is not None)
    print("\n".join(f"warning: {w}" for w in warnings))
    print("\n".join(errors) or f"glossary lint: no errors{f' at {tag}' if tag else ''}")
    sys.exit(1 if errors else 0)
