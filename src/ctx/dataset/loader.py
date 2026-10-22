"""``load``: read a dataset file, every table through its parser, every problem at once.

The tables' parsers live one family to a module beside this one; ``load`` calls them in one
fixed order, so a broken file lists its problems in the order the tables are read.
"""

import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ctx.dataset.analyses import (
    Contrast,
    DeSpec,
    FunctionSpec,
    NcrnaSpec,
    _contrasts,
    _de,
    _function,
    _ncrna,
)
from ctx.dataset.checks import DatasetError, _list, _need, _table
from ctx.dataset.counts import CountsSpec, _counts
from ctx.dataset.gdc import ClinicalSpec, GdcSpec, _clinical, _gdc
from ctx.dataset.proteo import ProteoSpec, _proteo
from ctx.dataset.sources import (
    REQUIRED_COLUMNS,
    Column,
    Resource,
    Sdrf,
    _columns,
    _resources,
    _sdrf,
)

TABLES = {
    "dataset", "resources", "guard", "sdrf", "columns", "derived", "order", "expect",
    "contrasts", "de", "function", "sanity", "gdc", "clinical", "proteo", "ncrna", "counts",
}  # fmt: skip
# The tables each module's own command refuses to start without, by Dataset field: the command
# checks them (``Dataset.lacks``) and ``ctx run`` plans the module from them, so neither restates
# the other. ``ctx ncrna`` and a cross-study ``ctx proteo`` read M5's results, not [de].
NEEDS = {"de": ("contrasts", "de"), "function": ("de", "function"), "proteo": ("proteo",),
         "ncrna": ("ncrna",), "clinical": ("gdc", "clinical")}  # fmt: skip

# The dataset file this process read for each slug: a figure's spec names it as an input
# (ctx.figures.record), and the modules that draw do not pass it on to ctx.viz.save.
_LOADED: dict[str, Path] = {}


def loaded() -> dict[str, Path]:
    """The dataset file ``load`` last read for each slug in this process."""
    return dict(_LOADED)


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
    gdc: GdcSpec | None = None
    clinical: ClinicalSpec | None = None
    proteo: ProteoSpec | None = None
    ncrna: NcrnaSpec | None = None
    counts: CountsSpec | None = None

    @property
    def name(self) -> str:
        """The name pages show: ``[dataset] label`` where given, else the slug.

        The slug names files, so it is lowercase with underscores; the label need not be.
        """
        return str(self.info.get("label", self.slug))

    @property
    def accessions(self) -> tuple[str, ...]:
        """``[dataset] accessions``: what fetchngs downloads, in the file's order."""
        return tuple(self.info.get("accessions", ()))

    def table(self) -> Path | None:
        """The own table of ``[sdrf] file``, beside the dataset file, never the working folder."""
        return self.path.parent / self.sdrf.file if self.sdrf and self.sdrf.file else None

    def lacks(self, module: str) -> tuple[str, ...]:
        """The tables of ``NEEDS`` the module's command needs and this file does not declare."""
        return tuple(table for table in NEEDS[module] if not getattr(self, table))

    def resource(self, name: str) -> Resource:
        """The resource with this name."""
        return next(r for r in self.resources if r.name == name)

    def parent(self) -> "Dataset | None":
        """The dataset this one is a subset of (``[dataset] subset_of``), or None for a whole one.

        Its file sits beside this one, named by its slug, as ``conf/datasets`` holds both; a
        page names it by its own ``name``, its label: M3's subset note named the slug, "tier1",
        where the parent's label is "Step-1" (kit v6, 2026-10-07).

        Raises:
            DatasetError: If that file is missing, unusable or names another slug.
        """
        slug = self.info.get("subset_of")
        if not slug:
            return None
        path = self.path.with_name(f"{slug}.toml")
        if not path.is_file():
            raise DatasetError(f"{self.path.name}: subset_of {slug!r}, but no {path} beside it")
        found = load(path)
        if found.slug != slug:
            raise DatasetError(f"{path}: slug {found.slug!r}, not the subset_of {slug!r}")
        return found


def _expected_counts(raw: Any, problems: list[str]) -> dict[str, dict[str, int]]:
    """[expect.counts]: per curated column, the runs each value must have, each a whole number."""
    out = {}
    for column, table in _table(raw, "expect.counts", problems).items():
        counts = _table(table, f"expect.counts.{column}", problems)
        # bool is an int to isinstance; a count written as true is still a mistake (F.6a).
        bad = {v: n for v, n in counts.items() if not isinstance(n, int) or isinstance(n, bool)}
        problems += [f"expect.counts.{column}.{v}: expected int, found {n!r}"
                     for v, n in bad.items()]  # fmt: skip
        out[column] = {str(v): int(n) for v, n in counts.items() if v not in bad}
    return out


def _routes(raw: Mapping[str, Any], gdc: GdcSpec | None, proteo: ProteoSpec | None,
            ncrna: NcrnaSpec | None, counts: CountsSpec | None) -> list[str]:  # fmt: skip
    """Tables that need a route into gather the file does not declare, or that clash."""
    problems = []
    if counts and gdc:
        problems.append("[counts] and [gdc] are two routes into gather; keep one")
    if "clinical" in raw and "gdc" not in raw:
        problems.append("[clinical] needs [gdc]: M10 reads a GDC project's gather outputs")
    # gather reads small-RNA libraries from a GDC project alone (--mirna-dir), so a file that
    # declared them elsewhere would list miRNA figures no root can hold.
    if ncrna and ncrna.mirna and "gdc" not in raw:
        problems.append("[ncrna.mirna] needs [gdc]: gather reads small-RNA libraries from a "
                        "GDC project alone")  # fmt: skip
    if proteo and proteo.design == "tumour_normal" and not (gdc and gdc.biospecimens):
        problems.append("proteo: the tumour_normal design needs [gdc] biospecimens, which "
                        "gather reads for each aliquot's case and sample type")  # fmt: skip
    return problems


def load(path: Path) -> Dataset:
    """Read and check a dataset file.

    Args:
        path: The TOML file.

    Returns:
        The parsed dataset.

    Raises:
        DatasetError: Listing every problem found, not only the first.
    """
    try:
        raw = tomllib.loads(path.read_text())
    except tomllib.TOMLDecodeError as error:
        # tomllib's message ends "(at line L, column C)"; no other check can run on the file.
        raise DatasetError(f"{path}:\n  not valid TOML: {error}") from None
    problems = [f"unknown table [{name}]" for name in sorted(set(raw) - TABLES)]
    info = raw.get("dataset", {})
    slug = _need(info, "slug", str, "dataset", problems)
    if "label" in info:
        _need(info, "label", str, "dataset", problems)
    # The runs `ctx guide` writes into the fetchngs list (Part 1), in the file's order.
    accessions = info.get("accessions", [])
    if not isinstance(accessions, list) or not all(isinstance(a, str) for a in accessions):
        problems.append(f"dataset.accessions: expected a list of accessions, found {accessions!r}")
    resources = _resources(raw.get("resources", []), problems)
    columns = _columns(raw.get("columns", {}), problems)
    sdrf = _sdrf(raw["sdrf"], resources, problems) if "sdrf" in raw else None
    named = {c.name for c in columns}
    if sdrf and (missing := [c for c in REQUIRED_COLUMNS if c not in named]):
        problems.append(f"columns: an SDRF dataset must define {missing}")
    bins = _table(raw.get("derived", {}).get("batch", {}).get("bins", {}), "derived.batch.bins",
                  problems)  # fmt: skip
    bins = {k: _list(v, f"derived.batch.bins.{k}", problems) for k, v in bins.items()}
    lengths = [n for levels in bins.values() for n in levels]
    if len(lengths) != len(set(lengths)):
        problems.append("derived.batch.bins: a read length sits in two bins")
    order = {k: tuple(_list(v, f"order.{k}", problems))
             for k, v in raw.get("order", {}).items()}  # fmt: skip
    problems += [f"order.{k}: no [columns.{k}]" for k in order if k not in named]
    expect = raw.get("expect", {})
    counts = _expected_counts(expect.get("counts", {}), problems)
    known = named | {"batch", "read_length"}
    problems += [f"expect.counts.{k}: not a curated column" for k in counts if k not in known]
    contrasts = _contrasts(raw.get("contrasts", []), problems)
    names = {c.name for c in contrasts}
    de = _de(raw["de"], names, problems) if "de" in raw else None
    function = _function(raw["function"], names, resources, problems) if "function" in raw else None
    gdc = _gdc(raw["gdc"], problems) if "gdc" in raw else None
    clinical = _clinical(raw["clinical"], problems) if "clinical" in raw else None
    proteo = _proteo(raw["proteo"], names, resources, problems) if "proteo" in raw else None
    ncrna = _ncrna(raw["ncrna"], names, problems) if "ncrna" in raw else None
    counts_spec = _counts(raw["counts"], problems) if "counts" in raw else None
    problems += _routes(raw, gdc, proteo, ncrna, counts_spec)
    if problems:
        raise DatasetError(f"{path}:\n  " + "\n  ".join(problems))
    _LOADED[slug] = path
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
        expect_counts=counts,
        contrasts=contrasts,
        de=de,
        function=function,
        sanity=raw.get("sanity", {}),
        gdc=gdc,
        clinical=clinical,
        proteo=proteo,
        ncrna=ncrna,
        counts=counts_spec,
    )
