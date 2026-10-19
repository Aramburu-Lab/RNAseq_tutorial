"""``[[contrasts]]``, ``[de]``, ``[function]`` and ``[ncrna]``: what M5, M6 and M8 fit and draw."""

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from ctx.dataset.checks import _list, _need, _sequence, _table, _values
from ctx.dataset.sources import Resource


@dataclass(frozen=True)
class Split:
    """A seeded split of individuals into two arbitrary groups, balanced within strata."""

    seed: int
    strata: str


@dataclass(frozen=True)
class Contrast:
    """One comparison: which samples, which two levels, which design, and the n it must have.

    Attributes:
        name: Identifier, used in output names.
        variable: Column holding the compared levels; ignored for a split.
        reference: Baseline level; the LFC is tested over reference.
        tested: The other level.
        expected_n: (reference, tested) samples after every step; a mismatch stops the fit.
        design: pydeseq2 formula over ``group`` and sheet columns, e.g. ``~donor + group``.
        filters: Column to allowed values; a sample must match every one.
        paired: Keep only donors with a sample at both levels.
        split: Replace the compared levels by a seeded split of individuals.
        label: The contrast in plain words, as titles, leads, captions and axes name it; empty
            reads as its two levels (Rytis 2026-10-01: never the identifier on a page).
    """

    name: str
    variable: str
    reference: str
    tested: str
    expected_n: tuple[int, int]
    design: str = "~group"
    filters: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    paired: bool = False
    split: Split | None = None
    label: str = ""


@dataclass(frozen=True)
class PcaSpec:
    """A PCA figure and its PERMANOVA metrics.

    Attributes:
        name: Identifier, used in output and metric names.
        filters: Column to allowed values.
        factors: Sheet columns to test; the first colours the points.
        symbol: Column setting the point shape, or None.
        strata: Permute labels only within this column (the individual), or None.
        one_per: Keep one sample per value of this column (the first run), or None.
    """

    name: str
    filters: Mapping[str, tuple[str, ...]]
    factors: tuple[str, ...]
    symbol: str | None = None
    strata: str | None = None
    one_per: str | None = None


@dataclass(frozen=True)
class DeSpec:
    """What the M5 page draws beyond each contrast's tables."""

    featured: str
    pca: tuple[PcaSpec, ...] = ()
    lessons: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    scatters: tuple[tuple[str, str], ...] = ()
    box_contrast: str | None = None
    box_genes: tuple[str, ...] = ()


@dataclass(frozen=True)
class FunctionSpec:
    """Which gene sets M6 tests, their licence line, and what the page features."""

    gene_sets: str
    licence: str
    featured_set: str
    running: tuple[str, str]


@dataclass(frozen=True)
class NcrnaSpec:
    """What M8 builds: the contrast it stratifies, the pairs its select offers, its own rules.

    Attributes:
        featured: The contrast whose M5 results and samples M8 reads.
        pairs: How many lncRNA and neighbour pairs the page's select offers.
        sanity: M8's pre-registered rules, a ``[sanity]``-shaped table with its
            ``registered_on`` time; apart from ``[sanity]``, whose record is fixed
            (tests/test_sanity_record.py), which finds its registration by any line starting
            ``registered = ``, so this one is named otherwise.
        mirna: The miRNA contrast of a GDC dataset (``[ncrna.mirna]``), fitted on the mature
            miRNA counts gather wrote; None for a dataset without small-RNA libraries.
        mirnas: The mature miRNAs, by miRBase name, whose targets the page tests and pairs.
    """

    featured: str
    pairs: int
    sanity: Mapping[str, Any] = field(default_factory=dict)
    mirna: Contrast | None = None
    mirnas: tuple[str, ...] = ()


def _split(raw: Any, where: str, problems: list[str]) -> Split | None:
    """A contrast's split table, seed and strata; absent or empty is no split."""
    if not raw:
        return None
    table = _table(raw, f"{where}.split", problems)
    seed = _need(table, "seed", int, f"{where}.split", problems) if table else None
    strata = _need(table, "strata", str, f"{where}.split", problems) if table else None
    return Split(seed, strata) if seed is not None and strata is not None else None


def _contrasts(raw: list[Any], problems: list[str]) -> tuple[Contrast, ...]:
    """Parse [[contrasts]]."""
    out = []
    for i, entry in enumerate(raw):
        where = f"contrasts[{i}]"
        fields = [_need(entry, key, str, where, problems)
                  for key in ("name", "variable", "reference", "tested")]  # fmt: skip
        n = _need(entry, "n", list, where, problems)
        filters = _values(entry.get("filters", {}), f"{where}.filters", problems)
        split = _split(entry.get("split"), where, problems)
        if n is not None and (len(n) != 2 or not all(isinstance(x, int) for x in n)):  # noqa: PLR2004 - a pair
            problems.append(f"{where}.n: expected [reference n, tested n], found {n!r}")
        if None in fields or n is None:
            continue
        name, variable, reference, tested = (str(f) for f in fields)
        out.append(Contrast(
            name, variable, reference, tested, expected_n=(n[0], n[-1]),
            design=entry.get("design", "~group"), filters=filters,
            paired=bool(entry.get("paired", False)), split=split,
            label=str(entry.get("label", "")),
        ))  # fmt: skip
    names = [c.name for c in out]
    problems += [
        f"contrasts: {n!r} appears twice" for n in sorted(set(names)) if names.count(n) > 1
    ]
    return tuple(out)


def _pcas(raw: list[Any], problems: list[str]) -> tuple[PcaSpec, ...]:
    """Parse [[de.pca]]: a name and factors each, the rest optional."""
    out = []
    for i, entry in enumerate(raw):
        where = f"de.pca[{i}]"
        found = len(problems)
        p = _table(entry, where, problems)
        filters = _values(p.get("filters", {}), f"{where}.filters", problems)
        if "name" not in p:
            problems.append(f"{where}.name: expected str, found None")
        factors = _sequence(p.get("factors"), f"{where}.factors", problems)
        # Any problem refuses the file, so a broken entry need not be built.
        if len(problems) == found:
            out.append(PcaSpec(p["name"], filters, tuple(factors), p.get("symbol"),
                               p.get("strata"), p.get("one_per")))  # fmt: skip
    return tuple(out)


def _de(raw: Mapping[str, Any], names: set[str], problems: list[str]) -> DeSpec:
    """Parse [de] and check it names only declared contrasts."""
    lessons = _table(raw.get("lessons", {}), "de.lessons", problems)
    scatters = _list(raw.get("scatters", []), "de.scatters", problems)
    if not all(isinstance(pair, list) and len(pair) == 2 for pair in scatters):  # noqa: PLR2004 - pairs
        problems.append(f"de.scatters: expected [[contrast, contrast], ...], found {scatters!r}")
        scatters = []
    spec = DeSpec(
        featured=raw.get("featured", ""),
        pca=_pcas(_list(raw.get("pca", []), "de.pca", problems), problems),
        lessons={k: tuple(_list(v, f"de.lessons.{k}", problems)) for k, v in lessons.items()},
        scatters=tuple((a, b) for a, b in scatters),
        box_contrast=raw.get("box_contrast"),
        box_genes=tuple(_sequence(raw.get("box_genes", []), "de.box_genes", problems)),
    )
    used = [spec.featured, *(c for v in spec.lessons.values() for c in v)]
    used += [c for pair in spec.scatters for c in pair] + (
        [spec.box_contrast] if spec.box_contrast else []
    )
    problems += [
        f"de: {c!r} is not a declared contrast" for c in dict.fromkeys(used) if c not in names
    ]
    return spec


def _function(raw: Mapping[str, Any], names: set[str], resources: tuple[Resource, ...],
              problems: list[str]) -> FunctionSpec | None:  # fmt: skip
    """Parse [function]."""
    fields = [_need(raw, k, str, "function", problems)
              for k in ("gene_sets", "licence", "featured_set")]  # fmt: skip
    running = raw.get("running", [])
    if not isinstance(running, list) or len(running) != 2:  # noqa: PLR2004 - a pair
        problems.append("function.running: expected [contrast, gene set]")
        return None
    if fields[0] is not None and fields[0] not in {r.name for r in resources}:
        problems.append(f"function.gene_sets: no [[resources]] entry named {fields[0]!r}")
    if running[0] not in names:
        problems.append(f"function.running: {running[0]!r} is not a declared contrast")
    if None in fields:
        return None
    gene_sets, licence, featured = (str(f) for f in fields)
    return FunctionSpec(gene_sets, licence, featured, running=(running[0], running[1]))


def _ncrna(raw: Mapping[str, Any], names: set[str], problems: list[str]) -> NcrnaSpec | None:
    """Parse [ncrna]."""
    featured = _need(raw, "featured", str, "ncrna", problems)
    pairs = _need(raw, "pairs", int, "ncrna", problems)
    if featured is not None and featured not in names:
        problems.append(f"ncrna.featured: {featured!r} is not a declared contrast")
    if pairs is not None and pairs < 1:
        problems.append(f"ncrna.pairs: {pairs}, at least 1")
    sanity = raw.get("sanity", {})
    if sanity and not isinstance(sanity.get("registered_on"), str):
        problems.append("ncrna.sanity.registered_on: the time the rules were registered")
    mirna = None
    if "mirna" in raw:
        parsed = _contrasts([raw["mirna"]], problems)
        mirna = parsed[0] if parsed else None
    mirnas = raw.get("mirnas", [])
    if not isinstance(mirnas, list) or not all(isinstance(m, str) for m in mirnas):
        problems.append("ncrna.mirnas: a list of mature miRNA names")
        mirnas = []
    if mirnas and mirna is None:
        problems.append("ncrna.mirnas: needs [ncrna.mirna], the contrast they come from")
    if featured is None or pairs is None:
        return None
    return NcrnaSpec(featured, pairs, sanity, mirna, tuple(mirnas))
