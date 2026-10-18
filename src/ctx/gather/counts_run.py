"""M1c for a count matrix: check it at the door, curate its samples, write gather's tables.

The sheet, counts and genes tables are the ones every route of gather writes, so qc, de,
function and report read them unchanged. Everything is read and checked before the first
file is written; a failed check writes the problem table alone (``ctx.gather.counts``).
"""

import functools
import gzip
import hashlib
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl

from ctx import runtime
from ctx.dataset import Dataset, DatasetError, load
from ctx.gather import counts, counts_rsem
from ctx.gather.curation import PROBLEM_COLUMNS, CurationError, Problem, curate, read_table
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table
from ctx.qc.checkpoint import gene_types
from ctx.report import text
from ctx.report.builder import KeyNumber, ModulePage, ParameterCard, Reproduction, ResultTable
from ctx.report.builder import write as write_page

MODULE = "gather"
PAGE = "gather_counts"
GTF_RESOURCE = "annotation_gtf"
# The sheet columns this route writes before the curated ones, in order.
FIXED = ("sample", "run_accession", counts.SOURCE_COLUMN)
# Rows left out shown by name on the page, the most counts first; the table lists them all.
LARGEST = 10


@dataclass(frozen=True)
class Checked:
    """Everything read and checked; nothing written yet."""

    raw: counts.Raw
    mapped: counts.Mapped
    sheet: pl.DataFrame
    releases: dict[str, str]
    rounding: counts_rsem.Rounding | None = None


def _resource_path(root: Path, spec: Dataset, name: str) -> Path | None:
    """Where ``ctx setup`` put a resource the dataset file declares, or None if it has none.

    Raises:
        DatasetError: If it is declared but not fetched.
    """
    if not any(r.name == name for r in spec.resources):
        return None
    path = root / spec.resource(name).dest
    if not path.is_file():
        raise DatasetError(f"{path} is missing; run ctx setup --dataset-file {spec.path}")
    return path


def gtf_release(gtf: Path) -> str:
    """The GTF's own description line, which names its release, or "not stated"."""
    with gzip.open(gtf, "rt") if gtf.suffix == ".gz" else gtf.open() as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if line.startswith("##description:"):
                return line.removeprefix("##description:").strip()
    return "not stated"


def _runs(samples: tuple[str, ...]) -> pl.DataFrame:
    """One run per sample column, in fetchngs's shape where curation reads it.

    A matrix holds no reads, so its layout and base counts are unknown: curation's read
    length is then null and no read-length batch can be derived.
    """
    names = list(samples)
    return pl.DataFrame({
        "run_accession": names, "sample": names, "experiment_accession": names,
        "library_layout": [""] * len(names),
        "base_count": pl.Series([None] * len(names), dtype=pl.Int64),
        "read_count": pl.Series([None] * len(names), dtype=pl.Int64),
    })  # fmt: skip


def _table_text(root: Path, spec: Dataset, releases: dict[str, str]) -> str:
    """The sample table as curation reads it: the own table, else the fetched SDRF.

    Raises:
        DatasetError: If the dataset file has no ``[sdrf]`` or its table is not there.
    """
    if spec.sdrf is None:
        raise DatasetError(f"{spec.path}: no [sdrf] table; the samples need a donor, cell type "
                           "and condition from your own table")  # fmt: skip
    own = spec.table()
    path = own or root / spec.resource(spec.sdrf.resource).dest
    if not path.is_file():
        raise DatasetError(f"{path} is missing; write it beside {spec.path.name}, one row per "
                           "sample column of the matrix")  # fmt: skip
    body, digest = read_table(path) if own else (path.read_text(), "")
    releases[f"sha256:{path.name}"] = digest or hashlib.sha256(path.read_bytes()).hexdigest()
    return body


def _curated(
    raw: counts.Raw, body: str, spec: Dataset
) -> tuple[pl.DataFrame | None, list[Problem]]:
    """The curated sheet of the matrix's samples, or the curation problems."""
    try:
        curated = curate(_runs(raw.samples), body, spec)
    except CurationError as error:
        return None, list(error.problems)
    return curated.drop("read_length"), []


def _entrez(raw: counts.Raw, root: Path, spec: Dataset, gtf: Path
            ) -> tuple[dict[str, set[str]] | None, list[Problem]]:  # fmt: skip
    """GENCODE's Entrez table read through the GTF, where the matrix holds Entrez IDs."""
    n = sum(1 for i in raw.ids if counts.kind_of(i) == "entrez")
    path = _resource_path(root, spec, counts.ENTREZ_RESOURCE)
    if n and path is None:
        detail = (f"{n:,} rows have Entrez gene IDs; add the {counts.ENTREZ_RESOURCE} "
                  "resource, GENCODE's metadata.EntrezGene of the GTF's release")  # fmt: skip
        return None, [("gene ids", "entrez", detail)]
    return (counts.entrez_genes(path, gtf) if n and path else None), []


def _hgnc(raw: counts.Raw, root: Path, spec: Dataset, releases: dict[str, str]
          ) -> dict[str, set[str]] | None:  # fmt: skip
    """HGNC's Ensembl IDs per NCBI Gene ID, where the matrix holds them and the file is declared."""
    path = _resource_path(root, spec, counts.HGNC_RESOURCE)
    if path is None or not any(counts.kind_of(i) == "entrez" for i in raw.ids):
        return None
    releases["NCBI Gene IDs GENCODE does not link"] = f"HGNC, {path.name}"
    return counts.hgnc_ensembl(path)


def _digest(path: Path) -> str:
    """SHA-256 of a file, or of a folder's ``sha256sum`` lines for its files, sorted by name."""
    if not path.is_dir():
        return _sha256(path)
    lines = "".join(f"{_sha256(f)}  {f.name}\n" for f in sorted(path.iterdir()) if f.is_file())
    return hashlib.sha256(lines.encode()).hexdigest()


def check(root: Path, spec: Dataset, matrix: Path) -> Checked:
    """Read and check the matrix, map its genes and curate its samples; nothing is written.

    Order: the annotation is found; the matrix is read; its values, repeated IDs and samples
    are checked; its IDs placed on the annotation; its samples joined to the table both ways.

    Raises:
        DatasetError: If the dataset file lacks ``[counts]`` or a resource is not fetched.
        counts.CountsError: With every problem the matrix and the curation found.
    """
    assert spec.counts is not None  # noqa: S101 - run() refuses a dataset file without one
    gtf = _resource_path(root, spec, GTF_RESOURCE)
    if gtf is None:
        raise DatasetError(f"{spec.path}: no {GTF_RESOURCE} resource; gene IDs map through it")
    releases = {
        "count matrix": matrix.name, f"sha256:{matrix.name}": _digest(matrix),
        "counts source": spec.counts.source, "counts url": spec.counts.url or "not given",
        "annotation the source names": spec.counts.annotation or "not given",
        "annotation": gtf_release(gtf),
        "dataset file": str(spec.path), f"sha256:{spec.path.name}": _sha256(spec.path),
    }  # fmt: skip
    kind = spec.counts.kind
    raw = counts_rsem.read(matrix, spec.counts)
    problems = [*raw.problems, *counts_rsem.kind_problems(raw, kind),
                *counts.value_problems(raw, kind), *counts.duplicate_problems(raw)]  # fmt: skip
    # Rounded only after the values passed as read: rounded log2 values would look whole.
    rounding = None
    if kind == counts.RSEM_EXPECTED:
        raw, rounding = counts_rsem.rounded(raw)
        releases["values"] = text.fill(PAGE, "rounding.readme", **vars(rounding))
    entrez, missing = _entrez(raw, root, spec, gtf)
    mapped = counts.place(raw, gene_types(gtf), entrez, _hgnc(raw, root, spec, releases))
    if mapped.counts.height == 0:
        missing.append(("gene ids", "all", "no row maps to one gene of the annotation"))
    curated, curation = _curated(raw, _table_text(root, spec, releases), spec)
    problems += missing + curation
    if problems or curated is None:
        raise counts.CountsError(problems)
    sheet = curated.with_columns(sample=pl.col("run_accession"),
                                 **{counts.SOURCE_COLUMN: pl.lit(matrix.name)})  # fmt: skip
    rest = [c for c in sheet.columns if c not in FIXED]
    return Checked(raw, mapped, sheet.select(*FIXED, *rest), releases, rounding)


def _sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def columns(key: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) per column of one table, from the page file."""
    return {name: (unit, meaning) for name, (unit, meaning) in text.raw(PAGE, key).items()}


def _tables(
    found: Checked, spec: Dataset
) -> list[tuple[str, pl.DataFrame, dict[str, tuple[str, str]]]]:
    """(description, frame, column dictionary) of every table this route writes."""
    meaning = "count.meaning" if found.rounding is None else "count.meaning_rsem"
    per_sample = {s: ("reads", text.fill(PAGE, meaning, sample=s))
                  for s in found.mapped.counts.columns[1:]}  # fmt: skip
    genes = columns("columns.genes")
    sheet = columns("columns.samplesheet") | {c.name: (c.unit, c.meaning) for c in spec.columns}
    return [
        ("samplesheet", found.sheet, {c: sheet[c] for c in found.sheet.columns}),
        ("counts", found.mapped.counts, {"gene_id": genes["gene_id"], **per_sample}),
        ("genes", found.mapped.genes, genes),
        ("genes_left_out", found.mapped.left_out, columns("columns.left_out")),
    ]


def _checks(found: Checked) -> ResultTable:
    """What the matrix held, checked, one row per check."""
    fill, word = functools.partial(text.fill, PAGE), functools.partial(text.word, PAGE)
    raw, mapped = found.raw, found.mapped
    sizes = mapped.counts.select(pl.exclude("gene_id")).sum().row(0)
    left = mapped.left_out
    share = float(left["total"].sum()) / max(float(np.nansum(raw.values)), 1.0) * 100
    kinds = ", ".join(f"{k} {n:,}" for k, n in sorted(mapped.kinds.items()))
    v = mapped.versions
    summary = ", ".join(name for name, _ in raw.summary) or word("checks.none")
    rows = (
        (word("checks.matrix"), fill("checks.matrix_value", rows=len(raw.ids),
                                     samples=len(raw.samples), format=raw.format)),
        (word("checks.values"), fill("checks.values_value", values=raw.values.size)
         if found.rounding is None else fill("checks.values_rsem", values=raw.values.size)),
        (word("checks.rounding"), word("checks.rounding_none") if found.rounding is None
         else fill("checks.rounding_value", **vars(found.rounding))),
        (word("checks.summary"), summary),
        (word("checks.kinds"), kinds),
        (word("checks.mapped"), fill("checks.mapped_value", genes=mapped.counts.height,
                                     left=left.height, share=share)),
        (word("checks.hgnc"), fill("checks.hgnc_value", rows=mapped.through_hgnc)),
        (word("checks.versions"), fill("checks.versions_value", **v)),
        (word("checks.annotation"), fill("checks.annotation_value",
                                         gtf=found.releases["annotation"],
                                         named=found.releases["annotation the source names"])),
        (word("checks.library"), fill("checks.library_value", median=float(np.median(sizes)),
                                      low=min(sizes), high=max(sizes))),
    )  # fmt: skip
    return ResultTable(word("checks.title"), text.words(PAGE, "checks.columns"), rows)


def _reasons(left: pl.DataFrame) -> ResultTable:
    """The rows left out, by reason: how many, and the counts they held."""
    word = functools.partial(text.word, PAGE)
    by = left.group_by("kind", "reason").agg(pl.len(), pl.col("total").sum())
    by = by.sort("kind", "reason")
    rows = tuple((k, w, f"{n:,}", f"{t:,.0f}") for k, w, n, t in by.rows())
    empty = ((word("left.none"), "", "", ""),)
    return ResultTable(word("left.title"), text.words(PAGE, "left.columns"), rows or empty,
                       numeric=(2, 3))  # fmt: skip


def _largest(left: pl.DataFrame) -> ResultTable:
    """The rows left out that held the most counts, so a reader sees which genes are lost."""
    word = functools.partial(text.word, PAGE)
    top = left.sort("total", "source_id", descending=[True, False]).head(LARGEST)
    rows = tuple((i, k, w, f"{t:,.0f}") for i, k, w, t in
                 top.select("source_id", "kind", "reason", "total").rows())  # fmt: skip
    empty = ((word("left.none"), "", "", ""),)
    title = text.fill(PAGE, "largest.title", n=LARGEST)
    return ResultTable(title, text.words(PAGE, "largest.columns"), rows or empty, numeric=(3,))


def _write_problems(error: counts.CountsError, out: Path, dataset: str) -> None:
    """The problem table beside where the sheet would have gone; no sheet is written."""
    frame = pl.DataFrame(error.problems, schema=list(PROBLEM_COLUMNS), orient="row")
    write_table(frame, out, stem(MODULE, f"{dataset}_counts_problems", 1, day=date.today()),
                columns=PROBLEM_COLUMNS, data_releases={})  # fmt: skip


def run(  # noqa: PLR0913 - one call from the CLI, every argument is a named option or context
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset: str,
    *,
    matrix: Path,
    dataset_file: Path,
) -> list[Path]:
    """Check a count matrix and its sample table, then write the sheet, counts, genes and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: Slug naming the dataset in output names.
        matrix: The count matrix, genes by samples.
        dataset_file: The dataset file, with ``[counts]``, ``[sdrf]`` and the annotation.

    Returns:
        The page, then the tables.

    Raises:
        DatasetError: If the dataset file is unusable for this route.
        counts.CountsError: After writing the problem table, on any failed check.
    """
    day = date.today()
    page_stem = stem(MODULE, f"{dataset}_overview", 1, day=day)
    spec = load(dataset_file)
    if spec.slug != dataset:
        raise DatasetError(f"{dataset_file} is for {spec.slug!r}, not --dataset {dataset!r}")
    if spec.counts is None:
        raise DatasetError(f"{dataset_file}: no [counts] table naming the matrix's source")
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    if not matrix.is_file() and not matrix.is_dir():
        raise DatasetError(f"--counts {matrix} is not a file or a folder of RSEM files")
    try:
        found = check(settings.root, spec, matrix.absolute())
    except counts.CountsError as error:
        _write_problems(error, out, dataset)
        raise
    files: list[Path] = []
    for name, frame, described in _tables(found, spec):
        files += write_table(frame, out, stem(MODULE, f"{dataset}_{name}", 1, day=day),
                             columns=described, data_releases=found.releases)  # fmt: skip
    usage, budget = runtime.measure(settings.started), runtime.budget(MODULE)
    word = functools.partial(text.word, PAGE)
    page = ModulePage(
        title=text.fill(PAGE, "title", dataset=spec.name), **text.skeleton(PAGE),
        parameters=cards, results=(_checks(found), _reasons(found.mapped.left_out),
                 _largest(found.mapped.left_out)),
        key_numbers=(
            KeyNumber(word("keys.samples"), f"{len(found.raw.samples):,}"),
            KeyNumber(word("keys.genes"), f"{found.mapped.counts.height:,}"),
            KeyNumber(word("keys.left_out"), f"{found.mapped.left_out.height:,}",
                      word("keys.left_out_detail")),
        ),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=(word("inputs"),),
        ),
    )  # fmt: skip
    return [write_page(page, out / f"{page_stem}.html"), *files]
