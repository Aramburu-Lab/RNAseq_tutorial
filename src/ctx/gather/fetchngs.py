"""Read one nf-core/fetchngs outdir, in the layout of the pinned release as its code writes it.

See D-006 for the layout facts and their sources. Nothing here writes into the outdir.
"""

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl

# The one definition. An upgrade is one deliberate change: this, a new fixture, and the
# versions file, which 1.13.0 renames.
FETCHNGS_VERSION = "1.12.0"
VERSIONS_FILE = "nf_core_fetchngs_software_mqc_versions.yml"
STATUSES = ("verified", "mismatch", "no_md5", "missing", "not_downloaded")


class FetchngsError(ValueError):
    """An outdir that cannot be used as it stands; the message lists every offender."""


@dataclass(frozen=True)
class Outdir:
    """What one fetchngs outdir holds, checked.

    Attributes:
        runs: One row per run: every samplesheet column as text, FASTQ paths relative to the
            outdir, a status per file, and ``run_index``, ``n_runs``, ``single_end``,
            ``ena_files``.
        unlisted_fastq: FASTQ files in ``fastq/`` that no row names.
        provenance: fetchngs version and commit, Nextflow version, and the SHA-256 of every
            pipeline_info file read.
    """

    runs: pl.DataFrame
    unlisted_fastq: tuple[str, ...]
    provenance: dict[str, str]


def _sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _version(outdir: Path) -> dict[str, str]:
    """Refuse any release but the pinned one; the versions file names it with its commit."""
    path = outdir / "pipeline_info" / VERSIONS_FILE
    if not path.is_file():
        raise FetchngsError(f"{path} is missing: not a fetchngs {FETCHNGS_VERSION} outdir")
    text = path.read_text()
    found = re.search(r"nf-core/fetchngs: v(\d+\.\d+\.\d+)-g([0-9a-f]+)", text)
    if found is None or found.group(1) != FETCHNGS_VERSION:
        seen = found.group(1) if found else "unreadable"
        raise FetchngsError(f"fetchngs {seen} wrote {outdir}; ctx reads {FETCHNGS_VERSION} only")
    nextflow = re.search(r"Nextflow: (\S+)", text)
    return {
        "fetchngs": found.group(1),
        "fetchngs_commit": found.group(2),
        "nextflow": nextflow.group(1) if nextflow else "not recorded",
        f"sha256:pipeline_info/{VERSIONS_FILE}": _sha256(path),
    }


def _params(outdir: Path) -> tuple[dict[str, object], dict[str, str]]:
    """The one params file; two mean a -resume into the same outdir, which mixes runs."""
    found = sorted((outdir / "pipeline_info").glob("params_*.json"))
    if len(found) != 1:
        raise FetchngsError(f"{outdir}: {len(found)} params_*.json files; expected one run")
    return json.loads(found[0].read_text()), {
        f"sha256:pipeline_info/{found[0].name}": _sha256(found[0])
    }


def _sheet(outdir: Path) -> pl.DataFrame:
    """samplesheet.csv as text, so accessions and md5 strings are never reinterpreted."""
    sheet = pl.read_csv(outdir / "samplesheet" / "samplesheet.csv", infer_schema=False)
    if "strandedness" not in sheet.columns:
        raise FetchngsError(
            f"{outdir}: no strandedness column; rerun with --nf_core_pipeline rnaseq"
        )
    return sheet


def _completeness(outdir: Path, params: dict[str, object], sheet: pl.DataFrame) -> list[str]:
    """Every input ID has metadata, and metadata and samplesheet list the same runs."""
    problems = []
    name = Path(str(params.get("input", ""))).name
    listed = outdir / name
    if not name or not listed.is_file():
        return [f"input list {name or '(none)'} is not in the outdir; copy it there"]
    ids = [line.strip() for line in listed.read_text().splitlines() if line.strip()]
    missing = [i for i in ids if not (outdir / "metadata" / f"{i}.runinfo_ftp.tsv").is_file()]
    problems += [f"input {i}: no metadata/{i}.runinfo_ftp.tsv" for i in missing]
    in_metadata = set()
    for path in (outdir / "metadata").glob("*.runinfo_ftp.tsv"):
        in_metadata |= set(pl.read_csv(path, separator="\t", infer_schema=False)["run_accession"])
    in_sheet = set(sheet["run_accession"])
    problems += [
        f"run {r}: in metadata/ but not in samplesheet.csv" for r in sorted(in_metadata - in_sheet)
    ]
    problems += [
        f"run {r}: in samplesheet.csv but not in metadata/" for r in sorted(in_sheet - in_metadata)
    ]
    return problems


def _mappings(outdir: Path, sheet: pl.DataFrame) -> list[str]:
    """id_mappings.csv must say what the samplesheet says: the same runs, the same values."""
    mappings = pl.read_csv(outdir / "samplesheet" / "id_mappings.csv", infer_schema=False)
    shared = [c for c in mappings.columns if c in sheet.columns]
    joined = mappings.select(shared).join(
        sheet.select(shared), on="run_accession", how="full", suffix="_sheet", nulls_equal=True
    )
    problems = []
    for row in joined.iter_rows(named=True):
        run = row["run_accession"] or row["run_accession_sheet"]
        differing = [c for c in shared if c != "run_accession" and row[c] != row[f"{c}_sheet"]]
        if row["run_accession"] is None or row["run_accession_sheet"] is None:
            problems.append(f"run {run}: in only one of id_mappings.csv and samplesheet.csv")
        elif differing:
            problems.append(f"run {run}: id_mappings.csv disagrees on {differing}")
    return problems


def _file_status(outdir: Path, name: str) -> str:
    """Status of one downloaded FASTQ, by file name inside fastq/, never by the sheet's path."""
    path = outdir / "fastq" / name
    if not path.is_file():
        return "missing"
    md5 = outdir / "fastq" / "md5" / f"{name}.md5"
    if not md5.is_file():
        return "no_md5"  # the sra-tools path writes none
    expected = md5.read_text().split()[0]
    return "verified" if hashlib.md5(path.read_bytes()).hexdigest() == expected else "mismatch"


def _files(outdir: Path, sheet: pl.DataFrame, downloaded: bool) -> pl.DataFrame:
    """Relative FASTQ paths and a status per file."""
    rows = []
    columns = ("fastq_1", "fastq_2", "library_layout")
    for fastq_1, fastq_2, layout in sheet.select(columns).iter_rows():
        out: dict[str, str] = {}
        for key, value in (("1", fastq_1), ("2", fastq_2)):
            name = Path(value).name if value else ""
            if not downloaded:
                # A run with no ENA link has no name here either; its layout says what is owed.
                owed = key == "1" or layout == "PAIRED"
                out |= {f"fastq_{key}": "", f"status_{key}": "not_downloaded" if owed else ""}
            else:
                out |= {
                    f"fastq_{key}": f"fastq/{name}" if name else "",
                    f"status_{key}": _file_status(outdir, name) if name else "",
                }
        rows.append(out)
    return pl.DataFrame(
        rows, schema=dict.fromkeys(("fastq_1", "status_1", "fastq_2", "status_2"), pl.String)
    )


def _runs(sheet: pl.DataFrame, files: pl.DataFrame) -> pl.DataFrame:
    """One row per run with its replicate structure and ENA file count."""
    number = pl.col("run_accession").str.extract(r"(\d+)$").cast(pl.Int64)
    return (
        sheet.drop("fastq_1", "fastq_2")
        .hstack(files)
        .with_columns(
            number.rank("ordinal").over("sample").cast(pl.Int32).alias("run_index"),
            pl.len().over("sample").cast(pl.Int32).alias("n_runs"),
            (pl.col("library_layout") == "SINGLE").alias("single_end"),
            pl.col("fastq_ftp")
            .fill_null("")
            .str.split(";")
            .list.eval(pl.element().filter(pl.element() != ""))
            .list.len()
            .cast(pl.Int32)
            .alias("ena_files"),
        )
    )


def read(outdir: Path) -> Outdir:
    """Check and read one fetchngs outdir, in a fixed order.

    Order: release version; samplesheet as text with its strandedness column; a single params
    file; completeness of input IDs, metadata and samplesheet; agreement of id_mappings.csv;
    FASTQ files by name and md5. The first three stop at once. The rest are collected and
    raised together, so one run lists every offender.

    Args:
        outdir: A fetchngs outdir, as fetched.

    Returns:
        The checked runs, unlisted FASTQ files and provenance.

    Raises:
        FetchngsError: If the outdir cannot be used as it stands.
    """
    provenance = _version(outdir)
    sheet = _sheet(outdir)
    params, params_hash = _params(outdir)
    provenance |= params_hash
    problems = _completeness(outdir, params, sheet) + _mappings(outdir, sheet)
    downloaded = not params.get("skip_fastq_download", False)
    files = _files(outdir, sheet, downloaded)
    for row in sheet.select("run_accession").hstack(files).iter_rows(named=True):
        problems += [
            f"run {row['run_accession']}: {row[f'fastq_{k}'] or 'fastq'} is {row[f'status_{k}']}"
            for k in ("1", "2")
            if row[f"status_{k}"] in {"missing", "mismatch"}
        ]
    if problems:
        raise FetchngsError(f"{outdir}:\n  " + "\n  ".join(problems))
    named = {Path(p).name for c in ("fastq_1", "fastq_2") for p in files[c] if p}
    on_disk = {p.name for p in (outdir / "fastq").glob("*.fastq.gz")} if downloaded else set()
    return Outdir(_runs(sheet, files), tuple(sorted(on_disk - named)), provenance)
