"""A gene selection: who a gene is, its sign, the link that carries it, and what fills it (D-049).

A gene is its unversioned Ensembl number, signed by its raw log2 FC in the source contrast:
STAT1, ENSG00000115415, down in Step-1's ndn_gcsf, is -115415 there. Plotly's point numbers
change whenever a volcano redraws, and a version suffix changes between annotation releases,
so neither can name a gene from one page to the next. The page twin is ``ctxSelect`` in
``templates/_select.js.j2``; ``tests/test_select.py`` drives both with the same cases.
"""

import math
import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from ctx.io.naming import DESCRIPTION
from ctx.report import text

VERSION = 1
# Longest fragment read, in UTF-16 units as a page's location.hash counts them, '#' included.
# Step-1's longest page list, the 17,830 genes with a padj in any of its 10 contrasts, encodes
# in 123,068 characters unsigned and the end mark's 4 (its 05_de results tables, measured for
# GS.1). What browsers carry on file:// was measured for the GS plan, outside the repo.
LIMIT = 500_000
# A fragment within LIMIT holds fewer than LIMIT / 2 numbers, each a digit and a comma at
# least, so a longer count is refused before it is formatted: grouping thousands by regex
# takes time quadratic in the number of digits.
_DIGITS = len(str(LIMIT // 2))
# ASCII case only, as the page's /i folds: Python alone would read a long s, U+017F, as S.
_ID = re.compile(r"ENSG([0-9]{11})(?:\.[0-9]+)?", re.IGNORECASE | re.ASCII)
# A dataset slug and a contrast name, each an output name's description: both then sit in a
# link unescaped, so encoding needs no percent codes to stay an exact inverse.
_SOURCE = re.compile(rf"{DESCRIPTION}:{DESCRIPTION}")
_HEAD = re.compile(r"#v=([^&]+)(?:&|\Z)")
# Written last by encode: no field holds an '&', so no proper prefix of a fragment ends with
# it, and a link cut anywhere, inside the last number too, lacks it.
END = "&end"
# The version is _HEAD's to check, against VERSION alone.
_FRAGMENT = re.compile(r"#v=[^&]+&src=([^&]*)&n=([^&]*)&g=([^&]*)" + re.escape(END))
# [0-9], never \d: Python's \d takes other scripts' digits, and the page's does not.
_COUNT = re.compile(r"0|[1-9][0-9]*")
_GENE = re.compile(r"-?[1-9][0-9]{0,10}")
_GROUPS = re.compile(r"\B(?=(?:[0-9]{3})+(?![0-9]))")
# A pasted word starts and ends with a letter or digit, so a full stop, dash or slash beside
# a name separates it. Inside, it holds . _ / -, as symbols do (THRA1/BTR in GENCODE v36,
# Y_RNA, LL22NC01-81G9.3); no symbol of Step-1's gene table or of GENCODE v36 starts or ends
# with one. A word joined by slashes that is no symbol is read as its _PARTs (ISG15/IFI6).
_WORD = re.compile(r"[A-Za-z0-9](?:[A-Za-z0-9._/-]*[A-Za-z0-9])?")
_PART = re.compile(r"[A-Za-z0-9](?:[A-Za-z0-9._-]*[A-Za-z0-9])?")
_LINE = re.compile(r"\r\n|\r|\n")
_HEADER = ("source", "time", "count")


class SelectionError(ValueError):
    """A fragment that cannot be read; the message is page text from ``ui.toml [select]``."""

    def __init__(self, key: str, **values: str) -> None:
        """Fill the message for ``key`` from ``ui.toml``; ``values`` fill its placeholders."""
        super().__init__(text.fill("ui", key, **values))
        self.key, self.values = key, values


class GeneListError(ValueError):
    """A ``--genes-from`` list file that is not UTF-8 text, or whose header contradicts it."""


def number(gene_id: str) -> int:
    """The unversioned Ensembl number of a gene ID in any case: ENSG00000123119.14 is 123119.

    Raises:
        ValueError: If the ID is not ENSG, eleven digits and an optional version. A
            ``_PAR_Y`` copy is refused, since it shares its number with the X copy.
    """
    if match := _ID.fullmatch(gene_id):
        return int(match[1])
    raise ValueError(f"{gene_id!r} is not an Ensembl gene ID (ENSG, 11 digits, a version)")


def ensembl(gene: int) -> str:
    """The unversioned Ensembl ID of a gene, signed or not: -123119 is ENSG00000123119."""
    return f"ENSG{abs(gene):011d}"


def group(value: int | str) -> str:
    """A whole number with commas between thousands, from its digits: 12796 is 12,796."""
    return _GROUPS.sub(",", str(value))


@dataclass(frozen=True)
class GeneTable:
    """The genes a page can select, made by ``table``.

    Attributes:
        numbers: Unversioned Ensembl numbers, one per row.
        symbols: Gene symbols; an empty one names no gene.
        lfc: Raw log2 FC in the source contrast, or None for a table that cannot sign.
        row: Each number's row.
        by_symbol: Each upper-case symbol's numbers, ascending.
    """

    numbers: tuple[int, ...]
    symbols: tuple[str, ...]
    lfc: tuple[float, ...] | None
    row: Mapping[int, int]
    by_symbol: Mapping[str, tuple[int, ...]]


def table(
    ids: Sequence[str], symbols: Sequence[str | None], lfc: Sequence[float | None] | None = None
) -> GeneTable:
    """The lookup a selection, the paste box and a list file resolve genes through.

    Args:
        ids: Ensembl gene IDs, versioned or not.
        symbols: The gene symbol of each, or None.
        lfc: The raw log2 FC of each in the source contrast, or None.

    Returns:
        The table.

    Raises:
        ValueError: If the columns differ in length, two IDs share a number, or a gene has
            no finite log2 FC to sign it.
    """
    numbers = tuple(number(i) for i in ids)
    if len(symbols) != len(numbers) or (lfc is not None and len(lfc) != len(numbers)):
        raise ValueError("gene IDs, symbols and log2 FC differ in length")
    if twice := [ensembl(n) for n, k in Counter(numbers).items() if k > 1]:
        raise ValueError(f"numbers naming more than one row: {len(twice)}, such as {twice[:3]}")
    finite = None if lfc is None else [v for v in lfc if v is not None and math.isfinite(v)]
    if finite is not None and len(finite) < len(numbers):
        raise ValueError("a gene without a finite log2 FC cannot be signed")
    by_symbol: dict[str, list[int]] = {}
    for n, s in zip(numbers, symbols, strict=True):
        if s:
            by_symbol.setdefault(s.upper(), []).append(n)
    return GeneTable(
        numbers=numbers,
        symbols=tuple(s or "" for s in symbols),
        lfc=None if finite is None else tuple(finite),
        row={n: i for i, n in enumerate(numbers)},
        by_symbol={k: tuple(sorted(v)) for k, v in by_symbol.items()},
    )


@dataclass(frozen=True)
class Selection:
    """One page's selection.

    Attributes:
        src: ``<dataset>:<contrast>``, the contrast whose raw log2 FC signs the genes.
        genes: Unversioned Ensembl numbers in ascending order, negated for log2 FC < 0.
    """

    src: str
    genes: tuple[int, ...]

    def __post_init__(self) -> None:
        """Refuse what no fragment carries, so ``encode`` and ``decode`` stay inverses."""
        numbers = [abs(g) for g in self.genes]
        if not _SOURCE.fullmatch(self.src) or not all(0 < n < 10**11 for n in numbers):
            raise ValueError(f"no fragment carries source {self.src!r} or these numbers")
        if any(a >= b for a, b in zip(numbers, numbers[1:], strict=False)):
            raise ValueError("a selection's numbers ascend without repeats")


def select(src: str, genes: Iterable[str | int], table: GeneTable) -> Selection:
    """A selection, in the one order: strip the version, attach the sign, drop repeats, sort.

    Args:
        src: ``<dataset>:<contrast>`` of the table's log2 FC.
        genes: Ensembl IDs, versioned or not, or their numbers, signed or not, as ``decode``
            gives them: a number's sign is stripped like a version, and the table signs it.
        table: The genes, with their raw log2 FC in ``src``.

    Returns:
        The selection; a versioned ID and its unversioned twin are one gene.

    Raises:
        ValueError: If the table has no log2 FC, lacks a gene, or no fragment carries
            ``src``.
    """
    if table.lfc is None:
        raise ValueError("a table without log2 FC cannot sign a selection")
    signed: dict[int, int] = {}
    for gene in genes:
        n = abs(gene) if isinstance(gene, int) else number(gene)
        if n not in table.row:
            raise ValueError(f"{ensembl(n)} is not in the gene table")
        signed[n] = -n if table.lfc[table.row[n]] < 0 else n
    return Selection(src, tuple(sorted(signed.values(), key=abs)))


def carries(fragment: str) -> bool:
    """Whether a fragment is meant as a selection; ``#results`` and other anchors are not."""
    return fragment.startswith("#v=")


def encode(selection: Selection) -> str:
    """The fragment ``#v=1&src=<dataset>:<contrast>&n=<count>&g=<numbers>&end`` of a selection.

    ``END`` closes it, so ``decode`` refuses the link cut short anywhere.

    Raises:
        SelectionError: If the fragment is longer than ``LIMIT``, so a page says so before
            it writes a link that ``decode`` refuses.
    """
    genes = ",".join(str(g) for g in selection.genes)
    count = len(selection.genes)
    return _within_limit(f"#v={VERSION}&src={selection.src}&n={count}&g={genes}{END}")


def _within_limit(fragment: str) -> str:
    """The fragment, if it is at most ``LIMIT`` UTF-16 units long, as ``location.hash`` counts."""
    # surrogatepass: a lone surrogate, which argv holds for a byte that is not UTF-8, is one
    # unit on the page too, and strict encoding would raise on it.
    length = len(fragment.encode("utf-16-le", "surrogatepass")) // 2
    if length > LIMIT:
        raise SelectionError("select.long", length=group(length), limit=group(LIMIT))
    return fragment


def _field(key: str, empty: str, found: str, **values: str) -> SelectionError:
    """A field's refusal: ``key`` quoting what the field holds, or ``empty`` if it holds none."""
    return SelectionError(key, found=found, **values) if found else SelectionError(empty, **values)


def decode(fragment: str) -> Selection:
    """The selection a fragment carries, read strictly: only what ``encode`` writes is read.

    Args:
        fragment: The fragment with its ``#``, as ``location.hash`` gives it.

    Returns:
        The selection.

    Raises:
        SelectionError: Naming the first problem: the length, the version, the end mark,
            the fields, the source, the count, a gene number, a repeat or the order.
    """
    _within_limit(fragment)
    if not (head := _HEAD.match(fragment)):
        raise SelectionError("select.form", version=str(VERSION))
    if head[1] != str(VERSION):
        raise SelectionError("select.version", found=head[1], expected=str(VERSION))
    # Before the fields, which end in END too: a link cut short is then named as cut, at
    # any point, where the field match would call it malformed. Inside the last number,
    # a cut keeps the count and its numbers well formed, so only the end mark can tell.
    if not fragment.endswith(END):
        raise SelectionError("select.cut")
    if not (fields := _FRAGMENT.fullmatch(fragment)):
        raise SelectionError("select.form", version=str(VERSION))
    src, n, g = fields.groups()
    if not _SOURCE.fullmatch(src):
        raise _field("select.source", "select.source_empty", src)
    if not _COUNT.fullmatch(n) or len(n) > _DIGITS:
        raise _field("select.n", "select.n_empty", n, digits=str(_DIGITS))
    items = g.split(",") if g else []
    # With the end mark checked, a count that differs is an edit: a gene added or removed by
    # hand. It counts the items that hold something, so _genes names an empty one as such.
    held = len(items) - items.count("")
    if n != str(held):
        raise SelectionError("select.count", expected=group(n), decoded=group(held))
    return Selection(src, _genes(items))


def _genes(items: list[str]) -> tuple[int, ...]:
    """The signed numbers of a fragment's ``g``, each checked in turn."""
    genes: list[int] = []
    seen: set[int] = set()
    for item in items:
        if not _GENE.fullmatch(item):
            raise _field("select.number", "select.number_empty", item)
        gene = int(item)
        if abs(gene) in seen:
            raise SelectionError("select.duplicate", gene=ensembl(gene))
        if genes and abs(gene) < abs(genes[-1]):
            raise SelectionError("select.order", gene=ensembl(gene))
        seen.add(abs(gene))
        genes.append(gene)
    return tuple(genes)


@dataclass(frozen=True)
class Pasted:
    """What a paste box or a list file names.

    Attributes:
        found: The numbers of the genes named, each once, in the order first named.
        ambiguous: Each symbol naming several genes, as first typed, with their IDs; none of
            its genes is added.
        missing: Each name no gene has, as first typed.
    """

    found: tuple[int, ...]
    ambiguous: tuple[tuple[str, tuple[str, ...]], ...]
    missing: tuple[str, ...]


def parse(pasted: str, table: GeneTable) -> Pasted:
    """The genes pasted text names: symbols or Ensembl IDs, any separator, any case.

    A ``#`` starts a comment to the end of its line, so a saved list's header is not read as
    names. On a line with an ID found in the table, a symbol naming that ID's gene, alone or
    with others, is its label: a saved list's ``ENSG00000202019 Y_RNA`` line is one gene, not
    an ambiguous symbol.

    Args:
        pasted: The text.
        table: The genes it may name.

    Returns:
        The genes found, the ambiguous symbols and the names not found.
    """
    found: dict[int, None] = {}
    ambiguous: dict[str, tuple[str, tuple[str, ...]]] = {}
    missing: dict[str, str] = {}
    for line in _LINE.split(pasted):
        words = _words(line.split("#", 1)[0], table)
        ids = {number(w) for w in words if _ID.fullmatch(w)} & table.row.keys()
        for word in words:
            if _ID.fullmatch(word):
                named: tuple[int, ...] = (n,) if (n := number(word)) in ids else ()
                key = f"#{n}"
            else:
                named, key = table.by_symbol.get(word.upper(), ()), word.upper()
                if ids & set(named):
                    continue
            if len(named) == 1:
                found.setdefault(named[0])
            elif named:
                ambiguous.setdefault(key, (word, tuple(ensembl(n) for n in named)))
            else:
                missing.setdefault(key, word)
    return Pasted(tuple(found), tuple(ambiguous.values()), tuple(missing.values()))


def _words(text: str, table: GeneTable) -> list[str]:
    """A line's words; one joined by slashes is its parts, unless it is a symbol (THRA1/BTR)."""
    words: list[str] = []
    for word in _WORD.findall(text):
        whole = "/" not in word or word.upper() in table.by_symbol
        words += [word] if whole else _PART.findall(word)
    return words


@dataclass(frozen=True)
class GeneList:
    """A ``--genes-from`` list file.

    Attributes:
        source: ``<dataset>:<contrast>`` from its ``# source:`` line, or None.
        time: Its ``# time:`` line as written, or None.
        genes: What its lines name, read as the paste box reads them.
    """

    source: str | None
    time: str | None
    genes: Pasted


def read_list(path: Path, table: GeneTable) -> GeneList:
    """Read a list file: the Save-list format, or any text the paste box reads.

    The Save-list format is ``# source:``, ``# time:`` and ``# count:`` lines, then one gene
    per line, its ID and symbol. Each header line is optional; a count is checked against the
    lines that name something, so a list cut short is refused, as a link is.

    Args:
        path: The file, UTF-8 with or without a byte-order mark.
        table: The genes it may name.

    Returns:
        The list.

    Raises:
        GeneListError: If the file is not UTF-8 text (a UTF-16 export, a Parquet table), a
            header line repeats, the source is not ``<dataset>:<contrast>``, or the count
            differs from the lines that name genes.
    """
    try:
        content = path.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as e:
        raise GeneListError(f"{path}: not UTF-8 text; a list file is UTF-8 text") from e
    header: dict[str, str] = {}
    lines = 0
    for line in _LINE.split(content):
        body, _, comment = line.partition("#")
        lines += bool(_WORD.search(body))
        key, colon, value = (part.strip() for part in comment.partition(":"))
        if not body.strip() and colon and key in _HEADER:
            if key in header:
                raise GeneListError(f"{path}: two '# {key}:' lines")
            header[key] = value
    source, count = header.get("source"), header.get("count")
    if source is not None and not _SOURCE.fullmatch(source):
        raise GeneListError(f"{path}: '# source: {source}' is not <dataset>:<contrast>")
    if count is not None and count != str(lines):
        raise GeneListError(f"{path}: '# count: {count}', but {lines} lines name genes")
    return GeneList(source, header.get("time"), parse(content, table))
