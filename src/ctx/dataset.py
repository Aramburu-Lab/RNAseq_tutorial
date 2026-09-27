"""The one reader of a dataset file (``conf/datasets/<slug>.toml``).

A dataset file says where a cohort's metadata comes from and how to read it: the SDRF join,
which columns become curation columns and how their values are spelled, the batch rule, and
the counts the curated sheet must reproduce. Every schema error is reported at once.
"""

import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Curation columns every dataset with an SDRF must fill (harmonise.CURATION).
REQUIRED_COLUMNS = ("donor", "cell_type", "condition")
TABLES = {
    "dataset", "resources", "guard", "sdrf", "columns", "derived", "order", "expect",
    "contrasts", "de", "function", "sanity",
}  # fmt: skip


class DatasetError(ValueError):
    """A dataset file that cannot be used; the message lists every problem."""


@dataclass(frozen=True)
class Resource:
    """A file ``ctx setup`` fetches: its mirrors, size, checksum and place under the root."""

    name: str
    urls: tuple[str, ...]
    bytes: int
    sha256: str
    dest: str


@dataclass(frozen=True)
class Column:
    """A curated column: the SDRF header it comes from, its unit, meaning and value map.

    ``values`` maps each SDRF value, whitespace-collapsed, to its curated spelling; None
    passes values through unchanged.
    """

    name: str
    source: str
    unit: str
    meaning: str
    values: Mapping[str, str] | None


@dataclass(frozen=True)
class Sdrf:
    """How SDRF rows meet the fetchngs runs: the join pair and the cross-checked pairs."""

    resource: str
    join: tuple[str, str]
    cross_check: tuple[tuple[str, str], ...]


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
class Dataset:
    """A parsed dataset file."""

    path: Path
    slug: str
    info: Mapping[str, str]
    resources: tuple[Resource, ...]
    sdrf: Sdrf | None
    columns: tuple[Column, ...]
    batch_bins: Mapping[str, tuple[int, ...]]
    order: Mapping[str, tuple[str, ...]]
    expect_runs: int | None
    expect_counts: Mapping[str, Mapping[str, int]]
    contrasts: tuple[Contrast, ...] = ()
    de: DeSpec | None = None
    function: FunctionSpec | None = None
    sanity: Mapping[str, Any] = field(default_factory=dict)

    def resource(self, name: str) -> Resource:
        """The resource with this name."""
        return next(r for r in self.resources if r.name == name)


def _need(table: Mapping[str, Any], key: str, kind: type, where: str, problems: list[str]) -> Any:
    """Read one required key of a given type, recording a problem instead of raising."""
    value = table.get(key)
    # bool is an int to isinstance; a count written as true is still a mistake.
    if not isinstance(value, kind) or (isinstance(value, bool) and kind is not bool):
        problems.append(f"{where}.{key}: expected {kind.__name__}, found {value!r}")
        return None
    return value


def _resources(raw: list[Any], problems: list[str]) -> tuple[Resource, ...]:
    """Parse [[resources]]."""
    out = []
    for i, entry in enumerate(raw):
        where = f"resources[{i}]"
        fields = [
            _need(entry, "name", str, where, problems),
            _need(entry, "urls", list, where, problems),
            _need(entry, "bytes", int, where, problems),
            _need(entry, "sha256", str, where, problems),
            _need(entry, "dest", str, where, problems),
        ]
        if None not in fields:
            name, urls, size, digest, dest = fields
            if Path(dest).is_absolute() or ".." in Path(dest).parts:
                problems.append(f"{where}.dest: {dest!r} must stay under the root")
            out.append(Resource(name, tuple(urls), size, digest, dest))
    return tuple(out)


def _columns(raw: Mapping[str, Any], problems: list[str]) -> tuple[Column, ...]:
    """Parse [columns.<name>]."""
    out = []
    for name, entry in raw.items():
        where = f"columns.{name}"
        source = _need(entry, "source", str, where, problems)
        meaning = _need(entry, "meaning", str, where, problems)
        values = entry.get("values")
        if values is not None and not all(isinstance(v, str) for v in values.values()):
            problems.append(f"{where}.values: every mapped value must be a string")
        if source is not None and meaning is not None:
            out.append(Column(name, source, entry.get("unit", "none"), meaning, values))
    return tuple(out)


def _sdrf(raw: Mapping[str, Any], resources: tuple[Resource, ...], problems: list[str]) -> Sdrf:
    """Parse [sdrf]."""
    resource = _need(raw, "resource", str, "sdrf", problems)
    join = _need(raw, "join", list, "sdrf", problems) or ["", ""]
    pairs = raw.get("cross_check", [])
    if resource is not None and resource not in {r.name for r in resources}:
        problems.append(f"sdrf.resource: no [[resources]] entry named {resource!r}")
    if len(join) != 2 or any(len(p) != 2 for p in pairs):  # noqa: PLR2004 - pairs
        problems.append("sdrf.join and each sdrf.cross_check entry must be [sdrf header, column]")
    return Sdrf(resource or "", (join[0], join[-1]), tuple((a, b) for a, b in pairs))


def _values(raw: Mapping[str, Any]) -> dict[str, tuple[str, ...]]:
    """A filter table: column to allowed values."""
    return {k: tuple(str(x) for x in v) for k, v in raw.items()}


def _contrasts(raw: list[Any], problems: list[str]) -> tuple[Contrast, ...]:
    """Parse [[contrasts]]."""
    out = []
    for i, entry in enumerate(raw):
        where = f"contrasts[{i}]"
        fields = [_need(entry, key, str, where, problems)
                  for key in ("name", "variable", "reference", "tested")]  # fmt: skip
        n = _need(entry, "n", list, where, problems)
        split = entry.get("split")
        if n is not None and (len(n) != 2 or not all(isinstance(x, int) for x in n)):  # noqa: PLR2004 - a pair
            problems.append(f"{where}.n: expected [reference n, tested n], found {n!r}")
        if None in fields or n is None:
            continue
        name, variable, reference, tested = (str(f) for f in fields)
        out.append(Contrast(
            name, variable, reference, tested, expected_n=(n[0], n[-1]),
            design=entry.get("design", "~group"),
            filters=_values(entry.get("filters", {})), paired=bool(entry.get("paired", False)),
            split=Split(int(split["seed"]), str(split["strata"])) if split else None,
        ))  # fmt: skip
    names = [c.name for c in out]
    problems += [
        f"contrasts: {n!r} appears twice" for n in sorted(set(names)) if names.count(n) > 1
    ]
    return tuple(out)


def _de(raw: Mapping[str, Any], names: set[str], problems: list[str]) -> DeSpec:
    """Parse [de] and check it names only declared contrasts."""
    pca = tuple(
        PcaSpec(
            p["name"],
            _values(p.get("filters", {})),
            tuple(p["factors"]),
            p.get("symbol"),
            p.get("strata"),
            p.get("one_per"),
        )  # fmt: skip
        for p in raw.get("pca", [])
    )
    spec = DeSpec(
        featured=raw.get("featured", ""),
        pca=pca,
        lessons={k: tuple(v) for k, v in raw.get("lessons", {}).items()},
        scatters=tuple((a, b) for a, b in raw.get("scatters", [])),
        box_contrast=raw.get("box_contrast"),
        box_genes=tuple(raw.get("box_genes", [])),
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
    if len(running) != 2:  # noqa: PLR2004 - a pair
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


def load(path: Path) -> Dataset:
    """Read and check a dataset file.

    Args:
        path: The TOML file.

    Returns:
        The parsed dataset.

    Raises:
        DatasetError: Listing every problem found, not only the first.
    """
    raw = tomllib.loads(path.read_text())
    problems = [f"unknown table [{name}]" for name in sorted(set(raw) - TABLES)]
    info = raw.get("dataset", {})
    slug = _need(info, "slug", str, "dataset", problems)
    resources = _resources(raw.get("resources", []), problems)
    columns = _columns(raw.get("columns", {}), problems)
    sdrf = _sdrf(raw["sdrf"], resources, problems) if "sdrf" in raw else None
    named = {c.name for c in columns}
    if sdrf and (missing := [c for c in REQUIRED_COLUMNS if c not in named]):
        problems.append(f"columns: an SDRF dataset must define {missing}")
    bins = raw.get("derived", {}).get("batch", {}).get("bins", {})
    lengths = [n for levels in bins.values() for n in levels]
    if len(lengths) != len(set(lengths)):
        problems.append("derived.batch.bins: a read length sits in two bins")
    order = {k: tuple(v) for k, v in raw.get("order", {}).items()}
    problems += [f"order.{k}: no [columns.{k}]" for k in order if k not in named]
    expect = raw.get("expect", {})
    counts = expect.get("counts", {})
    known = named | {"batch", "read_length"}
    problems += [f"expect.counts.{k}: not a curated column" for k in counts if k not in known]
    contrasts = _contrasts(raw.get("contrasts", []), problems)
    names = {c.name for c in contrasts}
    de = _de(raw["de"], names, problems) if "de" in raw else None
    function = _function(raw["function"], names, resources, problems) if "function" in raw else None
    if problems:
        raise DatasetError(f"{path}:\n  " + "\n  ".join(problems))
    return Dataset(
        path=path,
        slug=slug,
        info=info,
        resources=resources,
        sdrf=sdrf,
        columns=columns,
        batch_bins={level: tuple(v) for level, v in bins.items()},
        order=order,
        expect_runs=expect.get("runs"),
        expect_counts={k: {str(v): int(n) for v, n in c.items()} for k, c in counts.items()},
        contrasts=contrasts,
        de=de,
        function=function,
        sanity=raw.get("sanity", {}),
    )
