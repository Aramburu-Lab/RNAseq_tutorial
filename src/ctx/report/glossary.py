"""The glossary: ``glossary.yaml`` read once, looked up by id, term or alias, and explained.

Rytis's draft v3.12 (D-048). Page text marks jargon as ``[[term]]`` (``ctx.report.marks``); a
page shows the first mention per section with an info icon whose tooltip holds ``short`` and
a link to the glossary page, which shows every field. In a notebook, ``explain("phred")``
prints one entry. ``ctx.report.glossary_lint`` checks the file and every mark in the pages.
"""

import os
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from datetime import date
from functools import cache
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]
from markupsafe import Markup, escape

from ctx.io.naming import MODULES, stem
from ctx.report import marks, text

CATEGORIES = (
    "sequencing", "quality", "library", "genome", "alignment", "quantification", "qc_metrics",
    "statistics", "survival", "function", "proteomics", "data", "clinical", "tools",
)  # fmt: skip
# The draft's header marks no field of these optional; aliases may be an empty list.
REQUIRED = {"id": str, "term": str, "aliases": list, "category": str, "short": str, "long": str}
OPTIONAL = {"rna_vs_dna": str, "typical": str, "source": str, "see_also": list, "confirm": bool}
_ID = re.compile(r"[a-z0-9]+(?:_[a-z0-9]+)*")
# A name of three characters or fewer, or with two capitals (GO, OS, TPM, MultiQC), is matched
# as written, in a mark and in the lint's search for bare mentions: "go" and "os" are words,
# not the Gene Ontology or overall survival (L.6, D-104).
_EXACT = re.compile(r"^.{1,3}$|[A-Z].*[A-Z]")
# The glossary page sits beside the hub, in the report module's folder of each ctx root.
REPORT = "report"
# The icon: a circled "i" drawn in the text colour the button sets (the link colour).
ICON = (
    '<svg class="info-icon" viewBox="0 0 16 16" aria-hidden="true" focusable="false">'
    '<circle cx="8" cy="8" r="6.5"/><path d="M8 7.2v4.4M8 4.4v.2"/></svg>'
)


class GlossaryError(ValueError):
    """The glossary file breaks its schema, or a mark names no entry."""


@dataclass(frozen=True)
class Entry:
    """One glossary entry, as ``glossary.yaml`` holds it."""

    id: str
    term: str
    category: str
    short: str
    long: str
    aliases: tuple[str, ...] = ()
    rna_vs_dna: str = ""
    typical: str = ""
    source: str = ""
    see_also: tuple[str, ...] = ()
    confirm: bool = False

    def names(self) -> tuple[str, ...]:
        """Every name a mark may use for this entry: its id, its term and its aliases."""
        return (self.id, self.term, *self.aliases)

    def search_text(self) -> str:
        """What the glossary page's search reads: the entry's own words, not the page's labels."""
        fields = (self.term, *self.aliases, self.short, self.long, self.rna_vs_dna, self.typical,
                  self.source)  # fmt: skip
        return marks.plain(" ".join(f for f in fields if f)).replace("`", "")


def structure(raw: Any) -> list[str]:
    """Every schema problem: fields, types, ids, categories, names shared, see_also targets.

    A name (an id, term or alias, in any case) that two entries share makes a mark ambiguous,
    so it is a problem even when it is one entry's alias and another's term.
    """
    if not isinstance(raw, list) or not all(isinstance(e, dict) for e in raw):
        return ["the glossary is not a list of entries"]
    problems: list[str] = []
    owner: dict[str, Any] = {}
    for n, entry in enumerate(raw):
        where = f"entry {n + 1} ({entry.get('id', '?')})"
        problems += _fields(where, entry)
        for name in {str(v).lower() for v in _raw_names(entry)}:
            if owner.setdefault(name, entry.get("id")) != entry.get("id"):
                problems.append(f"{where}: {name!r} also names {owner[name]}")
    ids = [e.get("id") for e in raw]
    problems += [f"id {i!r} is used {ids.count(i)} times" for i in dict.fromkeys(ids)
                 if ids.count(i) > 1]  # fmt: skip
    for entry in raw:
        for target in entry.get("see_also") or []:
            if target not in ids or target == entry.get("id"):
                problems.append(f"{entry.get('id')}: see_also {target!r} is no other entry")
    return problems


def _fields(where: str, entry: Mapping[str, Any]) -> list[str]:
    """One entry's schema problems: unknown or missing fields, types, id and category."""
    fields = {**REQUIRED, **OPTIONAL}
    problems = [f"{where}: unknown field {k!r}" for k in entry if k not in fields]
    problems += [f"{where}: no {k!r}" for k in REQUIRED if k not in entry]
    for key, kind in fields.items():
        value = entry.get(key)
        # An optional field may be left empty (null); a required one may not.
        if not isinstance(value, kind) and (value is not None or key in REQUIRED):
            problems += [f"{where}: {key} is not {kind.__name__}"] if key in entry else []
        elif isinstance(value, list) and not all(isinstance(v, str) for v in value):
            problems.append(f"{where}: {key} holds a value that is not text")
        elif key in REQUIRED and isinstance(value, str) and not value.strip():
            problems.append(f"{where}: {key} is blank")
    if not _ID.fullmatch(str(entry.get("id", ""))):
        problems.append(f"{where}: id is not snake_case")
    if entry.get("category") not in CATEGORIES:
        problems.append(f"{where}: category {entry.get('category')!r} is not declared")
    return problems


def _raw_names(entry: Mapping[str, Any]) -> list[Any]:
    """An entry's id, term and aliases, as the file holds them."""
    return [entry.get("id", ""), entry.get("term", ""), *(entry.get("aliases") or [])]


def read(path: Path | None = None) -> Any:
    """The glossary file parsed, unchecked: the packaged one, or another for a test."""
    source = files("ctx.report").joinpath("glossary.yaml") if path is None else path
    return yaml.safe_load(source.read_text(encoding="utf-8"))


@cache
def load() -> dict[str, Entry]:
    """Every entry by id, in file order.

    Raises:
        GlossaryError: If the file breaks its schema; the message lists every problem.
    """
    raw = read()
    if problems := structure(raw):
        raise GlossaryError("glossary.yaml:\n" + "\n".join(problems))
    entries: dict[str, Entry] = {}
    for e in raw:
        values: dict[str, Any] = {k: tuple(v) if isinstance(v, list) else v for k, v in e.items()}
        entries[e["id"]] = Entry(**values)
    return entries


def exact(name: str) -> bool:
    """Whether a name matches only as written: three characters or fewer, or two capitals."""
    return bool(_EXACT.search(name))


@cache
def _index() -> tuple[dict[str, str], dict[str, str]]:
    """Every name to its entry's id: those matched as written, and the rest lower-cased."""
    names = [(name, e.id) for e in load().values() for name in e.names()]
    as_written = {name: i for name, i in names if exact(name)}
    return as_written, {name.lower(): i for name, i in names if not exact(name)}


def lookup(name: str) -> Entry:
    """The entry an id, term or alias names: as written if ``exact``, else in any case.

    Raises:
        GlossaryError: If no entry has that name.
    """
    as_written, any_case = _index()
    found = as_written.get(name.strip()) or any_case.get(name.strip().lower())
    if found is None:
        raise GlossaryError(f"{name!r} is not in the glossary")
    return load()[found]


def explain(term: str) -> None:
    """Print one entry for a notebook: term, category, short, long and every optional field.

    Args:
        term: An id, term or alias, such as ``"phred"`` or ``"Q score"``.
    """
    entry = lookup(term)
    labels = text.raw("ui", "glossary")
    lines = [f"{entry.term} ({labels['categories'][entry.category]})", entry.short, "", entry.long]
    for key in ("rna_vs_dna", "typical", "source"):
        if value := getattr(entry, key):
            lines += ["", f"{labels[key]} {value}"]
    if entry.aliases:
        lines += ["", f"{labels['aliases']} {', '.join(entry.aliases)}"]
    if entry.see_also:
        lines += [f"{labels['see_also']} {', '.join(load()[s].term for s in entry.see_also)}"]
    if entry.confirm:
        lines += ["", labels["unverified"]]
    print("\n".join(lines))


def grouped() -> tuple[tuple[str, tuple[Entry, ...]], ...]:
    """The entries by category, in the declared order of categories and the file's order."""
    entries = load().values()
    return tuple((c, tuple(e for e in entries if e.category == c)) for c in CATEGORIES)


def page_path(page: Path, day: date | None = None) -> Path | None:
    """The glossary page that a page at ``page`` links to (D-048), or None outside a ctx root.

    Pages sit in ``outputs/NN_<module>/`` of a ctx root, and the glossary beside the hub in
    ``outputs/09_report/``, dated and versioned like every output. A module page links
    ``../09_report/<name>``; a student kit copies the file with the pages that link it. A page
    written anywhere else, such as a figure ``ctx figure`` writes into the reader's own folder
    (Phase 2D), keeps its tooltips and links no glossary: none sits beside it.

    Args:
        page: Where the linking page is written.
        day: The date in the glossary's name; today when omitted.

    Returns:
        The glossary page's path, or None when ``page`` is not in a module folder.
    """
    folders = {f"{number:02d}_{module}" for module, number in MODULES}
    if page.parent.name not in folders:
        return None
    report = f"{dict(MODULES)[REPORT]:02d}_{REPORT}"
    name = stem(REPORT, "glossary", 1, day=day or date.today())
    return page.parent.parent / report / f"{name}.html"


def href(page: Path, glossary: Path) -> str:
    """The glossary's address as the page at ``page`` links it: relative, with slashes."""
    return Path(os.path.relpath(glossary, page.parent)).as_posix()


def used(texts: Iterable[str]) -> set[str]:
    """The ids of every entry the texts mark.

    Raises:
        GlossaryError: If a mark names no entry.
    """
    return {lookup(m.key).id for t in texts for m in marks.marks(t)}


def tooltips(ids: Iterable[str], link: str) -> dict[str, Any]:
    """What a page embeds for its tooltips: only the entries it marks, and the glossary link."""
    return {
        "href": link,
        "link": text.word("ui", "glossary.link"),
        "entries": {i: {"term": load()[i].term, "short": load()[i].short} for i in sorted(ids)},
    }


class Terms:
    """The marks of one page as HTML: the first mention per section carries the info icon."""

    def __init__(self) -> None:
        """Start a page with no term shown."""
        self._shown: set[tuple[str, str]] = set()

    def html(self, value: str, section: str) -> Markup:
        """Page text as HTML, each mark with its icon at its first mention in ``section``.

        Raises:
            GlossaryError: If a mark names no entry.
        """

        def term(mark: marks.Mark) -> str | None:
            """The term with its icon, or None when the section has shown it already."""
            entry = lookup(mark.key)
            if (section, entry.id) in self._shown and not mark.force:
                return None
            self._shown.add((section, entry.id))
            label = escape(text.fill("ui", "glossary.explain", term=entry.term))
            return (
                f'<span class="term">{escape(mark.display)}<button class="info" type="button" '
                f'data-term="{entry.id}" aria-label="{label}" aria-expanded="false">{ICON}'
                "</button></span>"
            )

        return marks.html(value, term)
