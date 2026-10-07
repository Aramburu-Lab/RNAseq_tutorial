"""A figure's lead whose numbers follow its filter: named aggregates over the rows shown (D-051).

A bar filter changes which rows a figure draws. Its lead keeps its words, from content, and
recomputes its numbers over the rows the reader keeps: in Python here, over every row, and in
the page by ``ctxFilter.aggregate`` (``templates/_filter.js.j2``), the twin of ``fill``, tested
under node against it at every filter level (``tests/test_page_script.py``). A number prints
through ``ctx.viz.bar.significant``, the rounding of the filter's live line; a count prints with
thousands separators; a word from a row prints as the row holds it. The caption keeps the view
every row gives (D-051, answer 4); ``Lead.note`` says so in the live line while a filter is set.
"""

import statistics
import string
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from typing import Any

from ctx.viz.bar import significant

KINDS = ("count", "median", "below", "min", "max", "mode")
Row = Mapping[str, Any]


@dataclass(frozen=True)
class Aggregate:
    """One named number of a lead, over the rows shown that match ``where``.

    Attributes:
        name: The placeholder it fills.
        kind: ``count`` the rows; ``below``, the rows whose ``column`` is below ``value``;
            ``median``, ``min`` or ``max`` of ``column``; ``mode``, its most common word, the
            first seen of a tie.
        column: The column it reads; none for a count.
        where: Column to value, compared as text as the filter compares: only rows matching
            every pair count, such as one outcome of a stacked bar.
        value: The bound of ``below``.
        row: For ``min`` and ``max``, (placeholder, column): the first such row's word in that
            column fills the placeholder.
        count: For ``mode``, the placeholder its count fills.
    """

    name: str
    kind: str
    column: str = ""
    where: Mapping[str, str] = field(default_factory=dict)
    value: float = 0.0
    row: tuple[str, str] | None = None
    count: str = ""

    def __post_init__(self) -> None:
        """Refuse an aggregate the page script could not compute the same way."""
        if self.kind not in KINDS:
            raise ValueError(f"{self.name}: kind {self.kind!r} is not one of {KINDS}")
        if (self.kind == "count") == bool(self.column):
            raise ValueError(f"{self.name}: a count reads no column; every other kind reads one")
        if self.row is not None and self.kind not in ("min", "max"):
            raise ValueError(f"{self.name}: only min and max name a row")
        if bool(self.count) != (self.kind == "mode"):
            raise ValueError(f"{self.name}: a mode, and only a mode, names its count")

    def names(self) -> tuple[str, ...]:
        """Every placeholder it fills."""
        row = (self.row[0],) if self.row else ()
        return (self.name, *row, *((self.count,) if self.count else ()))


def _matching(rows: Sequence[Row], where: Mapping[str, str]) -> list[Row]:
    """The rows whose every ``where`` column holds its value, as text."""
    return [r for r in rows if all(str(r.get(c)) == v for c, v in where.items())]


def _one(rows: Sequence[Row], spec: Aggregate) -> dict[str, str] | None:
    """One aggregate's placeholders as text, or None where it has no value to be taken over."""
    picked = _matching(rows, spec.where)
    if spec.kind == "count":
        return {spec.name: f"{len(picked):,}"}
    held = [r for r in picked if r.get(spec.column) is not None]
    if spec.kind == "below":
        return {spec.name: f"{sum(r[spec.column] < spec.value for r in held):,}"}
    if not held:
        return None
    if spec.kind == "median":
        return {spec.name: significant(statistics.median(r[spec.column] for r in held))}
    if spec.kind == "mode":
        counts: dict[str, int] = {}
        for r in held:
            counts[str(r[spec.column])] = counts.get(str(r[spec.column]), 0) + 1
        top = max(counts.values())
        word = next(w for w, n in counts.items() if n == top)
        return {spec.name: word, spec.count: f"{top:,}"}
    # Strictly smaller or larger only, so a tie keeps the first row, as the page does.
    best = held[0]
    for r in held[1:]:
        value, kept = r[spec.column], best[spec.column]
        if value < kept if spec.kind == "min" else value > kept:
            best = r
    out = {spec.name: significant(best[spec.column])}
    if spec.row:
        out[spec.row[0]] = str(best[spec.row[1]])
    return out


def fill(rows: Sequence[Row], aggregates: Sequence[Aggregate]) -> dict[str, str] | None:
    """Every aggregate's placeholders as text over ``rows``, or None if one has no value.

    Args:
        rows: The rows shown, in the page's order, each as the page embeds it.
        aggregates: What to compute.

    Returns:
        Placeholder to text; None when a median, minimum, maximum or mode has no row.
    """
    values: dict[str, str] = {}
    for spec in aggregates:
        found = _one(rows, spec)
        if found is None:
            return None
        values |= found
    return values


@dataclass(frozen=True)
class Lead:
    """A lead whose numbers the page recomputes over the rows a filter keeps.

    Attributes:
        template: The lead from content, its ``{placeholders}`` plain names.
        aggregates: The numbers recomputed over the rows shown.
        fixed: Placeholders that hold one text whatever is shown.
        note: What the live line adds while a filter is set: what the caption counts.
    """

    template: str
    aggregates: tuple[Aggregate, ...]
    fixed: Mapping[str, str] = field(default_factory=dict)
    note: str = ""

    def __post_init__(self) -> None:
        """Refuse a template the page could refill wrongly: a spec, a stray or missing name."""
        parsed = list(string.Formatter().parse(self.template))
        if any(spec or conversion for _, name, spec, conversion in parsed if name):
            raise ValueError("a lead's placeholders take text: no format spec or conversion")
        named = {n for spec in self.aggregates for n in spec.names()}
        if named & set(self.fixed):
            raise ValueError(f"filled twice: {sorted(named & set(self.fixed))}")
        used, given = {name for _, name, _, _ in parsed if name}, named | set(self.fixed)
        if used != given:
            raise ValueError(f"placeholders {sorted(used)}, values {sorted(given)}")
        # A span holds one placeholder: a mark or a code span around one would be cut by it.
        for literal, *_ in parsed:
            if literal.count("`") % 2 or literal.count("[[") != literal.count("]]"):
                raise ValueError("a placeholder inside a code span or a glossary mark")

    def text(self, rows: Sequence[Row]) -> str:
        """The lead filled over ``rows``.

        Raises:
            ValueError: If a number has no row to be taken over.
        """
        values = fill(rows, self.aggregates)
        if values is None:
            raise ValueError(f"the lead {self.template!r} has a number with no row")
        return self.template.format_map({**self.fixed, **values})

    def controls(self, rows: Sequence[Row]) -> dict[str, Any]:
        """What the page script refills the lead from, and the numbers over every row."""
        return {"template": self.template, "aggregates": [asdict(a) for a in self.aggregates],
                "fixed": dict(self.fixed), "note": self.note,
                "start": fill(rows, self.aggregates)}  # fmt: skip
