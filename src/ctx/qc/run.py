"""M2 QC (`ctx qc`): every sample judged by the registry, and the judgement tested.

The run reads the gather outputs, never the fetchngs layout (CONVENTIONS), measures each
stage it has inputs for, judges every metric the registry holds against its bands and its
stratum, writes one status row per run, and asks whether a metric, or a flag, differs between
the groups of a comparison M5 makes. Nothing is excluded.
"""

from collections.abc import Mapping
from dataclasses import dataclass, replace
from datetime import date
from pathlib import Path

import polars as pl

from ctx import runtime
from ctx.dataset import Dataset, DatasetError, load
from ctx.de.run import Writer, newest
from ctx.io.naming import MODULES, stem
from ctx.qc import checkpoint
from ctx.qc import registry as reg
from ctx.report import text
from ctx.report.builder import (
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    write,
)

MODULE = "qc"
PAGE = MODULE
GTF_RESOURCE = "annotation_gtf"
# The page lists this many confounding tests, smallest adjusted p first; the CSV holds all.
SHOWN_TESTS = 10
ALPHA = 0.05


@dataclass(frozen=True)
class Options:
    """What `ctx qc` was asked to do, beside the dataset.

    Attributes:
        dataset_file: The dataset file: contrasts, level orders and resources.
        registry: The registry file; None resolves it (``registry.resolve``).
        strata: Sheet columns whose shared values make a stratum, for the MAD bands.
        peers: Sheet columns whose shared values make a peer group.
        mad: The MAD multiplier; None keeps the registry's.
    """

    dataset_file: Path
    registry: Path | None = None
    strata: tuple[str, ...] = ("cell_type",)
    peers: tuple[str, ...] = ("cell_type", "condition")
    mad: float | None = None


@dataclass(frozen=True)
class Stage:
    """One stage's judgement: its values, bands and flags, all keyed by sample."""

    values: pl.DataFrame
    bands: pl.DataFrame
    flags: pl.DataFrame
    notes: tuple[str, ...] = ()


def load_registry(options: Options) -> reg.Registry:
    """The registry, with the MAD multiplier the reader chose."""
    registry = reg.load(reg.resolve(options.registry))
    return registry if options.mad is None else replace(registry, mad_multiplier=options.mad)


def gene_types(root: Path, spec: Dataset) -> tuple[pl.DataFrame | None, str]:
    """The annotation's gene types, or None and why: no such resource, or not fetched."""
    found = [r for r in spec.resources if r.name == GTF_RESOURCE]
    if not found:
        return None, text.fill(PAGE, "notes.no_annotation", resource=GTF_RESOURCE)
    path = root / found[0].dest
    if not path.is_file():
        return None, text.fill(PAGE, "notes.not_fetched", path=str(path))
    return checkpoint.gene_types(path), ""


def counts_stage(
    settings: runtime.Settings, spec: Dataset, gathered: Mapping[str, pl.DataFrame],
    registry: reg.Registry, options: Options,
) -> tuple[Stage, pl.DataFrame | None]:  # fmt: skip
    """The checkpoint on counts, and the sex check's own table."""
    sheet, counts = gathered["samplesheet"], gathered["counts"]
    if missing := [c for c in {*options.strata, *options.peers} if c not in sheet.columns]:
        raise DatasetError(f"--strata and --peers name columns the sheet lacks: {missing}")
    types, note = gene_types(settings.root, spec)
    values, sexes = checkpoint.sample_metrics(counts, sheet, gathered["genes"], types,
                                              options.peers, settings.threads)  # fmt: skip
    samples = sheet.unique("sample", keep="first", maintain_order=True)
    bands, flags = checkpoint.judge(values, samples, registry, options.strata)
    return Stage(values, bands, flags, (note,) if note else ()), sexes


def columns_of(frame: pl.DataFrame, registry: reg.Registry) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of a values table, from the registry's own words."""
    by_column = {m.column: (m.unit, m.meaning) for m in registry.metrics.values()}
    own = {c: (u, text.word(PAGE, f"columns.{c}")) for c, u in
           (("sample", "none"), ("run_accession", "none"))}  # fmt: skip
    return {c: own.get(c) or by_column[c] for c in frame.columns}


def _write_tables(writer: Writer, values: pl.DataFrame, frames: Mapping[str, pl.DataFrame],
                  registry: reg.Registry) -> None:  # fmt: skip
    """The values table, its columns described by the registry, then every other table.

    Each other table's columns, with their units, are described in ``content/qc.toml``; a
    column it does not describe stops the run.
    """
    writer.table(values, "metrics", columns_of(values, registry))
    for what, frame in frames.items():
        units = text.raw(PAGE, f"tables.{what}.units")
        writer.table(frame, what, {c: (units.get(c, "none"),
                                       text.word(PAGE, f"tables.{what}.columns.{c}"))
                                   for c in frame.columns})  # fmt: skip


def _status_table(status: pl.DataFrame) -> ResultTable:
    """The flagged runs, as the page lists them."""
    flagged = status.filter(pl.col("qc_status") != "pass")
    rows = tuple((r, s, q, reason) for r, s, q, reason, _ in flagged.rows())
    return ResultTable(text.fill(PAGE, "status.title", flagged=flagged.height, runs=status.height),
                       tuple(text.words(PAGE, "status.columns")), rows)  # fmt: skip


def _confounding_table(tests: pl.DataFrame) -> ResultTable:
    """The confounding tests with the smallest adjusted p; the CSV holds every one."""
    tested = tests.filter(pl.col("p_adj").is_not_null()).sort("p_adj").head(SHOWN_TESTS)
    rows = tuple(
        (r["contrasts"], r["metric"], r["test"], f"{r['effect']:.2f}", f"{r['p_adj']:.2g}",
         text.word(PAGE, "confounding.tied") if r["read_length_tied"] else "")
        for r in tested.iter_rows(named=True)
    )  # fmt: skip
    total = tests["p_adj"].is_not_null().sum()
    title = text.fill(PAGE, "confounding.title", shown=len(rows), tests=total)
    columns = tuple(text.words(PAGE, "confounding.columns"))
    return ResultTable(title, columns, rows, numeric=(3, 4))


def _keys(stage: Stage, status: pl.DataFrame, tests: pl.DataFrame) -> tuple[KeyNumber, ...]:
    """The numbers at the top of Results."""
    flagged = status.filter(pl.col("qc_status") != "pass")
    below = tests.filter(pl.col("p_adj") < ALPHA).height
    return (
        KeyNumber(text.word(PAGE, "keys.samples"), f"{stage.values.height:,}"),
        KeyNumber(text.word(PAGE, "keys.flagged"), f"{flagged.height:,}",
                  text.fill(PAGE, "keys.flagged_detail", runs=status.height)),
        KeyNumber(text.word(PAGE, "keys.metrics"), f"{stage.bands['metric'].n_unique():,}"),
        KeyNumber(text.word(PAGE, "keys.confounded"), f"{below:,}",
                  text.fill(PAGE, "keys.confounded_detail", alpha=ALPHA)),
    )  # fmt: skip


def run(
    settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str, dataset: str,
    options: Options,
) -> list[Path]:  # fmt: skip
    """Judge the dataset's samples and write the M2 tables and page.

    Returns:
        The page, then every table.

    Raises:
        DatasetError: If the dataset file does not match the slug, or gather has not run.
        RegistryError: If the registry cannot be found or read.
    """
    spec = load(options.dataset_file)
    if spec.slug != dataset:
        raise DatasetError(f"{options.dataset_file}: slug {spec.slug!r}, not {dataset!r}")
    registry = load_registry(options)
    gathered_dir = settings.root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    inputs = {w: newest(gathered_dir, dataset, w) for w in ("samplesheet", "counts", "genes")}
    gathered = {w: pl.read_parquet(p) for w, p in inputs.items()}
    stage, sexes = counts_stage(settings, spec, gathered, registry, options)
    sheet = gathered["samplesheet"]
    status = checkpoint.status(sheet, stage.flags, registry)
    tests = checkpoint.confounding(stage.values, stage.flags, sheet, spec.contrasts, spec.order,
                                   registry)  # fmt: skip
    day = date.today()
    releases = {"registry": str(registry.path), **{f"input:{w}": p.name for w, p in inputs.items()}}
    writer = Writer(settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset,
                    day, releases)  # fmt: skip
    frames = {"bands": stage.bands, "flags": stage.flags, "status": status, "confounding": tests}
    _write_tables(writer, stage.values, frames | ({"sex_check": sexes} if sexes is not None
                                                  else {}), registry)  # fmt: skip
    results: list[Result] = [*stage.notes, _status_table(status), _confounding_table(tests)]
    usage, budget = runtime.measure(settings.started), runtime.budget(MODULE)
    page = ModulePage(
        title=text.fill(PAGE, "title", dataset=spec.name), **text.skeleton(PAGE),
        parameters=cards, results=tuple(results), key_numbers=_keys(stage, status, tests),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=tuple(p.name for p in inputs.values()),
        ),
    )  # fmt: skip
    out = writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"
    return [write(page, out), *writer.files]
