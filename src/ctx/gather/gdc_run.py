"""M1c for a GDC project: check the folder, choose one file per case, write the tables and page.

Everything is read and checked before the first file is written (``ctx.gather.gdc``); the
outputs follow ``harmonise.run``'s pattern under ``outputs/01_gather/``.
"""

import dataclasses
import functools
import hashlib
import time
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl

from ctx import runtime
from ctx.dataset import Dataset, DatasetError, GdcSpec, load
from ctx.gather import gdc, gdc_mirna
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table
from ctx.proteo import annotate
from ctx.qc import registry
from ctx.report import text
from ctx.report.builder import (
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    ResultTable,
    write,
)

MODULE = "gather"
PAGE = "gather_gdc"
# (unit, meaning) of every sheet column, from the page file.
SHEET = "columns.samplesheet"


@dataclass(frozen=True)
class Read:
    """Everything gather read from one folder, checked, before anything is written."""

    folder: gdc.Folder
    counts: gdc.Counts
    sheet: pl.DataFrame
    annotations: pl.DataFrame
    cases: pl.DataFrame
    treatments: pl.DataFrame
    tests: pl.DataFrame
    cdr: pl.DataFrame | None
    md5_seconds: float
    aliquots: pl.DataFrame | None = None


@dataclass(frozen=True)
class MirnaRead:
    """Both miRNA kinds of the project, checked, and their sheet; nothing written yet."""

    matrices: gdc_mirna.Matrices
    sheet: pl.DataFrame
    md5_seconds: float


def read_mirna(mirna_dir: Path, isoform_dir: Path, spec: GdcSpec, star: Read) -> MirnaRead:
    """Read and check both miRNA folders, choose one library per patient, join the STAR files.

    Raises:
        gdc.GdcError: On any failed check; nothing is written.
    """
    if spec.mirna_workflow is None:
        raise gdc.GdcError("--mirna-dir needs [gdc] mirna_workflow in the dataset file")
    folders = {
        kind: gdc_mirna.read_kind(path.resolve(), spec, kind)
        for kind, path in ((gdc_mirna.MIRNA, mirna_dir), (gdc_mirna.ISOFORM, isoform_dir))
    }
    started = time.perf_counter()
    for folder in folders.values():
        gdc.verify(folder)
    seconds = time.perf_counter() - started
    found = gdc_mirna.read_matrices(folders[gdc_mirna.MIRNA], folders[gdc_mirna.ISOFORM])
    notes = gdc.annotations(star.folder.path / gdc.ANNOTATIONS, found.files)
    sheet = gdc_mirna.sheet(found, notes, star.sheet, spec)
    return MirnaRead(found, sheet.select(*columns("columns.mirna_samplesheet")), seconds)


def _mirna_tables(found: MirnaRead) -> list[tuple[str, pl.DataFrame, dict[str, tuple[str, str]]]]:
    """The miRNA sheet, both count matrices and the hairpins' loci, each with its columns."""
    m, names = found.matrices, found.sheet["sample"].to_list()
    meaning = text.raw(PAGE, "columns.mirna_counts")
    per_file = {s: ("reads", meaning["count"].format(sample=s)) for s in names}
    tables = []
    for what, key, ids, matrix in (("mirna_mature_counts", "mimat", m.mature, m.mature_counts),
                                   ("mirna_hairpin_counts", "hairpin", m.hairpins,
                                    m.hairpin_counts)):  # fmt: skip
        data = {key: list(ids)} | {s: pl.Series(s, matrix[:, j], dtype=pl.Int32)
                                   for j, s in enumerate(names)}  # fmt: skip
        tables.append((what, pl.DataFrame(data), {key: tuple(meaning[key]), **per_file}))
    return [("mirna_samplesheet", found.sheet, columns("columns.mirna_samplesheet")), *tables,
            ("mirna_hairpins", m.loci, columns("columns.mirna_hairpins"))]  # fmt: skip


def _mirna_checks(found: MirnaRead, spec: GdcSpec) -> ResultTable:
    """What the two miRNA folders hold, checked, one row per check."""
    fill = functools.partial(text.fill, PAGE)
    s, m = found.sheet, found.matrices
    share = {r: (pl.col(f"{r}_reads") / pl.col("mirna_reads") * 100).median()
             for r in gdc_mirna.REGIONS}  # fmt: skip
    medians = s.select(**share, cross=(pl.col("cross_mapped_reads") / pl.col("mirna_reads")
                                       * 100).median()).row(0, named=True)  # fmt: skip
    taken = gdc.chosen(s)
    types = ", ".join(t for _, _, t in gdc_mirna.KINDS.values())
    values = (
        fill("mirna.query_value", project=spec.project, workflow=spec.mirna_workflow, types=types),
        fill("mirna.files_value", files=s.height, seconds=found.md5_seconds),
        fill("mirna.hairpins_value", hairpins=len(m.hairpins),
             detected=int((m.hairpin_counts.sum(axis=1) > 0).sum())),
        fill("mirna.mature_value", mature=len(m.mature), isoforms=s["mature_isoforms"].median()),
        fill("mirna.sum_value", files=s.height, reads=int(s["mirna_reads"].sum())),
        fill("mirna.regions_value", **{r: medians[r] for r in gdc_mirna.REGIONS}),
        fill("mirna.cross_value", pct=medians["cross"]),
        fill("mirna.taken_value", tumours=gdc.chosen(s, spec.tumour).height,
             normals=gdc.chosen(s, spec.normal).height, taken=taken.height,
             joined=taken["rna_sample"].is_not_null().sum()),
    )  # fmt: skip
    rows = tuple(zip(text.words(PAGE, "mirna.rows"), values, strict=True))
    return ResultTable(text.word(PAGE, "mirna.title"), text.words(PAGE, "checks.columns"), rows)


def columns(key: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) per column of one table, from the page file."""
    return {name: (unit, meaning) for name, (unit, meaning) in text.raw(PAGE, key).items()}


def build_sheet(folder: gdc.Folder, counts: gdc.Counts, notes: pl.DataFrame,  # noqa: PLR0913 - each a part of the sheet
                spec: GdcSpec, metric: registry.Metric,
                libraries: pl.DataFrame | None = None) -> pl.DataFrame:  # fmt: skip
    """One row per counts file: identity, barcode parts, strand, library size and the choice.

    ``libraries`` is ``gdc.library_selection``'s record; without one the column stays empty.
    """
    stats = gdc.strand(counts.stats, metric)
    hit = notes.filter(pl.col("action") == "drops its files")["samples"]
    redacted = {s for joined in hit for s in joined.split(";")}
    parts = pl.DataFrame(
        [gdc.barcode(s, c) for s, c in folder.files.select("sample", "case").rows()],
        schema=dict.fromkeys(("vial", "portion", "plate", "site"), pl.String),
    )
    condition = (pl.when(pl.col("sample_type") == spec.tumour).then(pl.lit("tumour"))
                 .when(pl.col("sample_type") == spec.normal).then(pl.lit("normal"))
                 .otherwise(pl.lit("other")))  # fmt: skip
    frame = folder.files.hstack(parts).join(stats, on="sample", maintain_order="left")
    frame = frame.with_columns(condition=condition, batch=pl.col("plate"), donor=pl.col("case"),
                               run_accession=gdc.run_accession(), project=pl.lit(spec.project),
                               redacted=pl.col("sample").is_in(sorted(redacted)))  # fmt: skip
    frame = gdc.choose(frame, spec.protocol, spec.exclude_cases)
    if libraries is None:
        frame = frame.with_columns(library_selection=pl.lit(None, dtype=pl.String))
    else:
        frame = frame.join(libraries, on="sample", how="left", maintain_order="left")
    return frame.select(*columns(SHEET))


def read_all(path: Path, spec: Dataset) -> Read:
    """Read and check the whole folder; nothing is written.

    Raises:
        gdc.GdcError: On any failed check.
        registry.RegistryError: If the QC registry lacks the strand metric.
    """
    assert spec.gdc is not None  # noqa: S101 - run() refuses a dataset file without [gdc]
    metric = registry.load(registry.resolve(None)).metrics.get(gdc.STRAND_METRIC)
    if metric is None:
        raise registry.RegistryError(f"the QC registry has no [metrics.{gdc.STRAND_METRIC}]")
    folder = gdc.read_folder(path, spec.gdc)
    started = time.perf_counter()
    gdc.verify(folder)
    md5_seconds = time.perf_counter() - started
    cdr = None
    if spec.gdc.cdr_file and spec.gdc.cdr_type:
        cdr = gdc.cdr(path / spec.gdc.cdr_file, spec.gdc.cdr_type)
    counts = gdc.read_counts(folder, spec.gdc)
    notes = gdc.annotations(path / gdc.ANNOTATIONS, folder.files)
    sheet = build_sheet(
        folder, counts, notes, spec.gdc, metric, gdc.library_selection(path, spec.gdc)
    )
    cases, treatments, tests = gdc.clinical(path / gdc.CASES)
    aliquots = gdc.biospecimens(path, spec.gdc) if spec.gdc.biospecimens else None
    return Read(folder, counts, sheet, notes, cases, treatments, tests, cdr, md5_seconds,
                aliquots)  # fmt: skip


def read_projects(paths: list[Path], spec: Dataset) -> list[tuple[GdcSpec, Read]]:
    """Each folder read with its own project's spec, in the dataset file's project order.

    A folder's project is the one its file query names (``gdc.project_of``), so the folders may
    be given in any order; each declared project needs exactly one.

    Raises:
        gdc.GdcError: If a folder is of an undeclared project, two are of one, or one is
            missing (``gdc.by_project``), every problem named.
    """
    assert spec.gdc is not None  # noqa: S101 - run() refuses a dataset file without [gdc]
    by, problems = gdc.by_project(paths, spec.gdc.projects)
    if problems:
        raise gdc.GdcError("; ".join(problems))
    found = []
    for project in spec.gdc.projects:
        one = spec.gdc.one(project)
        found.append((one, read_all(by[project], dataclasses.replace(spec, gdc=one))))
    return found


def combine(found: list[Read]) -> Read:
    """One sheet and one matrix over every project's folder, files in project then aliquot order.

    ``run_accession`` is numbered again over the combined sheet, so no two files share one; the
    choice of one file per case is each folder's own, as no case is in two projects.

    Raises:
        gdc.GdcError: If the projects' counts files list different genes, or share a case.
    """
    if len(found) == 1:
        return found[0]
    first = found[0]
    if any(not r.counts.genes.equals(first.counts.genes) for r in found[1:]):
        raise gdc.GdcError("the projects' counts files list different genes")
    cases = [set(r.sheet["case"]) for r in found]
    if len(set().union(*cases)) < sum(len(c) for c in cases):
        raise gdc.GdcError("a case has files in two projects")
    files = pl.concat([r.folder.files for r in found])
    folder = gdc.Folder(first.folder.path, files, first.folder.release,
                        sum(r.folder.dropped_site for r in found))  # fmt: skip
    counts = dataclasses.replace(first.counts,
                                 matrix=np.hstack([r.counts.matrix for r in found]),
                                 stats=pl.concat([r.counts.stats for r in found]))  # fmt: skip
    sheet = pl.concat([r.sheet for r in found]).with_columns(run_accession=gdc.run_accession())

    def stacked(name: str) -> pl.DataFrame:
        """One clinical or annotation table over every project; a column one lacks is null."""
        frames: list[pl.DataFrame] = [getattr(r, name) for r in found]
        return pl.concat(frames, how="diagonal_relaxed")

    return Read(folder, counts, sheet, stacked("annotations"), stacked("cases"),
                stacked("treatments"), stacked("tests"), None,
                sum(r.md5_seconds for r in found))  # fmt: skip


def _flow(sheet: pl.DataFrame) -> ResultTable:
    """Files, cases and each drop of the rule, per sample type, in the rule's order."""
    word = functools.partial(text.word, PAGE)
    kinds = [k for k in ("tumour", "normal", "other") if k in set(sheet["condition"])]
    cells: dict[str, list[str]] = {k: [] for k in kinds}
    for kind in kinds:
        part = sheet.filter(pl.col("condition") == kind)
        kept = gdc.chosen(part)["case"].n_unique()
        reasons = Counter(part["drop_reason"])
        candidates = part.filter(pl.col("drop_reason").is_in(["", "not_first"]))
        ties = candidates.group_by("case").len().filter(pl.col("len") > 1).height
        cells[kind] = [f"{part.height:,}", f"{part['case'].n_unique():,}",
                       *(f"{reasons.get(r, 0):,}"
                         for r in ("excluded", "ffpe", "protocol", "redacted")),
                       f"{kept:,}", f"{part['case'].n_unique() - kept:,}", f"{ties:,}"]  # fmt: skip
    steps = text.words(PAGE, "flow.steps")
    rows = tuple((step, *(cells[k][i] for k in kinds)) for i, step in enumerate(steps))
    heads = (word("flow.step"), *(word(f"flow.{k}") for k in kinds))
    return ResultTable(word("flow.title"), heads, rows, numeric=tuple(range(1, len(heads))))


def _shares(sheet: pl.DataFrame) -> str:
    """Files per strand call with the range of their shares."""
    by = sheet.group_by("strand_call").agg(
        pl.len(),
        pl.col("first_strand_share").min().alias("lo"),
        pl.col("first_strand_share").max().alias("hi"),
    )
    return "; ".join(text.fill(PAGE, "checks.strand_value", call=c, n=n, lo=lo, hi=hi)
                     for c, n, lo, hi in by.sort("strand_call").rows())  # fmt: skip


def _no_feature(sheet: pl.DataFrame) -> str | None:
    """Median N_noFeature share of FFPE files against the files kept, where there are FFPE."""
    share = pl.col("n_nofeature") / (pl.col("n_nofeature") + pl.col("library_size")) * 100
    ffpe = sheet.filter(pl.col("preservation") == gdc.FFPE)
    if ffpe.is_empty():
        return None
    kept = gdc.chosen(sheet)
    return text.fill(PAGE, "checks.no_feature_value", ffpe=ffpe.select(share.median()).item(),
                     n=ffpe.height, kept=kept.select(share.median()).item())  # fmt: skip


def _checks(found: Read, spec: GdcSpec, several: bool = False) -> ResultTable:
    """What the folder holds, checked, one row per check; titled by project when ``several``.

    A dataset joining projects shows one table per project, and two tables under one title
    read as one table printed twice (S.4, the lung page).
    """
    fill, word = functools.partial(text.fill, PAGE), functools.partial(text.word, PAGE)
    files, genes = found.sheet.height, found.counts.genes.height
    sizes = found.sheet["library_size"]
    rows = [
        (word("checks.release"), fill("checks.release_value", release=spec.release)),
        (word("checks.query"), fill("checks.query_value", project=spec.project,
                                    workflow=spec.workflow)),
        (word("checks.manifest"), fill("checks.manifest_value", files=files)),
        (word("checks.md5"), fill("checks.md5_value", files=files, seconds=found.md5_seconds)),
        (word("checks.model"), fill("checks.model_value", model=spec.gene_model, files=files,
                                    rows=found.counts.rows, n=len(gdc.N_ROWS))),
        (word("checks.genes"), fill("checks.genes_value", par_y=found.counts.par_y, genes=genes)),
        (word("checks.stored"), fill("checks.stored_value", genes=genes, files=files,
                                     bytes=found.counts.matrix.nbytes)),
        (word("checks.strand"), _shares(found.sheet)),
        (word("checks.library"), fill("checks.library_value", column=spec.protocol,
                                      median=sizes.median(), low=sizes.min(), high=sizes.max())),
    ]  # fmt: skip
    if (no_feature := _no_feature(found.sheet)) is not None:
        rows.append((word("checks.no_feature"), no_feature))
    if spec.primary_site:
        rows.append((word("checks.site"), fill("checks.site_value", site=spec.primary_site,
                                               n=found.folder.dropped_site)))  # fmt: skip
    clinical = fill("checks.clinical_value", cases=found.cases.height,
                    with_counts=found.sheet["case"].unique().is_in(found.cases["case"]).sum(),
                    all=found.sheet["case"].n_unique())  # fmt: skip
    if found.cdr is not None:
        clinical += fill("checks.cdr_value", rows=found.cdr.height, type=spec.cdr_type)
    rows.append((word("checks.clinical"), clinical))
    title = fill("checks.title_project", project=spec.project) if several else word("checks.title")
    return ResultTable(title, text.words(PAGE, "checks.columns"), tuple(rows))


def _pairs(sheet: pl.DataFrame) -> str:
    """Cases with a chosen tumour and a chosen normal, of those with files of both."""

    def cases(frame: pl.DataFrame, kind: str) -> set[str]:
        """The cases with a file of one condition in a frame."""
        return set(frame.filter(pl.col("condition") == kind)["case"])

    kept = cases(gdc.chosen(sheet), "tumour") & cases(gdc.chosen(sheet), "normal")
    both = cases(sheet, "tumour") & cases(sheet, "normal")
    return text.fill(PAGE, "flow.pairs", kept=len(kept), both=len(both))


def _annotations(notes: pl.DataFrame) -> ResultTable:
    """Annotations of the cases' biospecimens by category: how many, on how many, the action."""
    by = (notes.group_by("classification", "category", "action")
          .agg(pl.len().alias("n"), pl.col("case").n_unique().alias("cases"),
               pl.col("files").sum().alias("files"))
          .sort("classification", "category", "action"))  # fmt: skip
    rows = tuple((c, k, f"{n:,}", f"{cases:,}", f"{f:,}", a) for c, k, a, n, cases, f in by.rows())
    word = functools.partial(text.word, PAGE)
    empty = ((word("annotations.none"), *("",) * 5),)
    return ResultTable(word("annotations.title"), text.words(PAGE, "annotations.columns"),
                       rows or empty, numeric=(2, 3, 4))  # fmt: skip


def _keys(found: Read, spec: GdcSpec) -> tuple[KeyNumber, ...]:
    """Files, cases, tumours taken, release."""
    word = functools.partial(text.word, PAGE)
    tumours = gdc.chosen(found.sheet, spec.tumour).height
    return (
        KeyNumber(word("keys.files"), f"{found.sheet.height:,}"),
        KeyNumber(word("keys.cases"), f"{found.sheet['case'].n_unique():,}"),
        KeyNumber(word("keys.tumours"), f"{tumours:,}", word("keys.tumours_detail")),
        KeyNumber(word("keys.release"), spec.release.split(" - ")[0]),
    )


def _matrix(found: Read) -> pl.DataFrame:
    """Genes by files, int32, gene_id then one column per aliquot in sheet order."""
    names = found.folder.files["sample"].to_list()
    assert names == found.sheet["sample"].to_list()  # noqa: S101 - both in aliquot order
    data = {"gene_id": found.counts.genes["gene_id"]}
    data |= {name: pl.Series(name, found.counts.matrix[:, j], dtype=pl.Int32)
             for j, name in enumerate(names)}  # fmt: skip
    return pl.DataFrame(data)


def _write(  # noqa: PLR0913 - each a part of the write
    found: Read, out: Path, dataset: str, day: date, spec: Dataset,
    releases: Mapping[str, str], *, mirna: MirnaRead | None = None,
) -> list[Path]:  # fmt: skip
    """Every table, each with its README; the miRNA tables too where both folders were read."""
    assert spec.gdc is not None  # noqa: S101 - run() checked
    matrix = _matrix(found)
    per_file = {s: ("reads", text.fill(PAGE, "columns.count", protocol=spec.gdc.protocol,
                                       sample=s))
                for s in matrix.columns[1:]}  # fmt: skip
    # The sheet's own columns, then those the dataset file declares (ctx.proteo.annotate).
    declared = columns(SHEET) | {c.name: (c.unit, c.meaning) for c in spec.columns}
    tables = [
        ("samplesheet", found.sheet, declared),
        ("counts", matrix, {"gene_id": columns("columns.genes")["gene_id"], **per_file}),
        ("genes", found.counts.genes, columns("columns.genes")),
        ("annotations", found.annotations, columns("columns.annotations")),
        ("cases", found.cases, columns("columns.cases")),
        ("treatments", found.treatments, columns("columns.treatments")),
        ("molecular_tests", found.tests, columns("columns.molecular_tests")),
    ]
    if found.cdr is not None:
        tables.append(("cdr", found.cdr, {c: ("none", text.fill(PAGE, "columns.cdr", column=c))
                                          for c in found.cdr.columns}))  # fmt: skip
    if found.aliquots is not None:
        tables.append(("aliquots", found.aliquots, columns("columns.aliquots")))
    tables += _mirna_tables(mirna) if mirna is not None else []
    files: list[Path] = []
    for name, frame, described in tables:
        files += write_table(frame, out, stem(MODULE, f"{dataset}_{name}", 1, day=day),
                             columns=described, data_releases=releases)  # fmt: skip
    return files


def _releases(path: Path, spec: Dataset, one: GdcSpec, prefix: str = "") -> dict[str, str]:
    """The releases and files every table derives from, with the SHA-256 of each record read.

    ``one`` is the folder's own project; ``prefix`` names it before each folder-specific key
    where the sheet joins several projects.
    """
    named = [gdc.METADATA, gdc.MANIFEST, gdc.STATUS, gdc.CASES, gdc.ANNOTATIONS]
    named += [n for n in (one.cdr_file, one.biospecimens) if n]
    return {
        "GDC": one.release, f"{prefix}GDC folder": str(path), f"{prefix}project": one.project,
        "workflow": one.workflow, "gene model": one.gene_model,
        "dataset file": str(spec.path),
        f"sha256:{spec.path.name}": hashlib.sha256(spec.path.read_bytes()).hexdigest(),
        **{f"sha256:{prefix}{n}": hashlib.sha256((path / n).read_bytes()).hexdigest()
           for n in named},
    }  # fmt: skip


def _mirna_releases(found: MirnaRead) -> dict[str, str]:
    """The miRNA folders and the SHA-256 of each one's metadata, manifest and status."""
    out = {}
    for kind, path in ((gdc_mirna.MIRNA, found.matrices.paths[0]),
                       (gdc_mirna.ISOFORM, found.matrices.paths[1])):  # fmt: skip
        metadata, manifest, _ = gdc_mirna.KINDS[kind]
        out[f"GDC {kind} folder"] = str(path)
        out |= {f"sha256:{kind}/{n}": hashlib.sha256((path / n).read_bytes()).hexdigest()
                for n in (metadata, manifest, gdc.STATUS)}  # fmt: skip
    return out


def run(  # noqa: PLR0913 - one call from the CLI, every argument is a named option or context
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset: str,
    *,
    gdc_dirs: list[Path],
    dataset_file: Path,
    mirna_dir: Path | None = None,
    isoform_dir: Path | None = None,
) -> list[Path]:
    """Check a GDC folder, choose one file per case and sample type, write the tables and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: Slug naming the dataset in output names.
        gdc_dirs: The GDC project folders as the fetch wrote them, one per declared project.
        dataset_file: The dataset file, with a ``[gdc]`` table.
        mirna_dir: The project's miRNA folder as the fetch wrote it; given with isoform_dir.
        isoform_dir: The project's isoform folder, the same aliquots' files.

    Returns:
        The page, then the tables.

    Raises:
        DatasetError: If the dataset file has no ``[gdc]`` or names another slug, or miRNA
            folders are given for several projects.
        gdc.GdcError: On any failed check of the folder; nothing is written.
    """
    day = date.today()
    page_stem = stem(MODULE, f"{dataset}_overview", 1, day=day)
    spec = load(dataset_file)
    if spec.slug != dataset:
        raise DatasetError(f"{dataset_file} is for {spec.slug!r}, not --dataset {dataset!r}")
    if spec.gdc is None:
        raise DatasetError(f"{dataset_file}: no [gdc] table; not a GDC dataset")
    if (mirna_dir is None) != (isoform_dir is None):
        raise DatasetError("--mirna-dir and --isoform-dir are given together or not at all")
    if mirna_dir is not None and len(spec.gdc.projects) > 1:
        raise DatasetError("--mirna-dir reads one project; this dataset file names several")
    projects = read_projects(gdc_dirs, spec)
    found = combine([r for _, r in projects])
    if spec.columns:
        found = dataclasses.replace(
            found, sheet=annotate.sheet_columns(found.sheet, spec, settings.root)
        )
    mirna = None
    if mirna_dir is not None and isoform_dir is not None:
        mirna = read_mirna(mirna_dir, isoform_dir, spec.gdc, found)
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    several = len(projects) > 1
    releases: dict[str, str] = {}
    for one, read in projects:
        releases |= _releases(read.folder.path, spec, one, f"{one.project} " if several else "")
    if mirna is not None:
        releases |= _mirna_releases(mirna)
    files = _write(found, out, dataset, day, spec, releases, mirna=mirna)
    word = functools.partial(text.word, PAGE)
    results: tuple[ResultTable | str, ...] = (*(_checks(r, one, several) for one, r in projects),
                                              _flow(found.sheet), _pairs(found.sheet),
                                              _annotations(found.annotations))  # fmt: skip
    if mirna is not None:
        results = (*results, _mirna_checks(mirna, spec.gdc))
    usage, budget = runtime.measure(settings.started), runtime.budget(MODULE)
    # The page's own input names the TCGA-CDR workbook; a folder read without one (CPTAC) has none.
    variant = None if found.cdr is not None else "variants.no_cdr"
    inputs = word(f"{variant}.inputs" if variant else "inputs")
    page = ModulePage(
        title=text.fill(PAGE, "title", dataset=spec.name), **text.skeleton(PAGE, variant),
        parameters=cards, results=results, key_numbers=_keys(found, spec.gdc),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=(inputs,),
        ),
    )  # fmt: skip
    return [write(page, out / f"{page_stem}.html"), *files]
