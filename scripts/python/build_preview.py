#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side; builds the course preview, every module on made-up data until it is built
# ///
r"""Build the course preview: one page per module M0 to M10 and the hub, on one made-up dataset.

Built modules run as a student runs them, through ``ctx``, on the test suite's synthetic
design (tests/synth.py): counts with planted differential genes, a batch and three planted QC
problems. A module not built yet draws its page from ``scripts/python/preview_<module>.py``,
whose ``page(root, dataset, day)`` writes ``outputs/NN_<module>/<day>_<module>_<dataset>_
overview_v1.html`` from made-up data and the page text in ``content/<module>.toml``. Every
preview page says, in its first result, that its numbers are made up.

Usage:
    pixi run python scripts/python/build_preview.py --root ~/ctx_preview/synthetic
    pixi run python scripts/python/build_preview.py --root <dir> --only process
"""

from __future__ import annotations

import argparse
import importlib
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / "tests"), str(Path(__file__).parent)]
CTX = str(Path(sys.executable).parent / "ctx")
DATASET = "desynth"
# Course order; each is a MODULES token with no command yet (D-044).
PREVIEWS = ("process", "quant", "proteo", "ncrna", "genecard", "clinical")
# The first result of every preview page: text.word("ui", "preview.note"); M3, built but
# short of reads and a genome here, says why in content/process.toml's preview.note.
NOTE = ("ui", "preview.note")
# The course's batch names (D-043) in place of the test data's short and long.
BATCH = {"short": "Reads ≤ 76 nt (min 35 nt)", "long": "Reads > 76 nt (max 500 nt)"}


def ctx(root: Path, *args: str) -> None:
    """Run one ctx command in root, timed; stop the build on a failure."""
    start = time.perf_counter()
    done = subprocess.run([CTX, "--root", str(root), *args], cwd=root, capture_output=True,
                          text=True, check=False)  # fmt: skip
    print(f"{time.perf_counter() - start:6.1f} s  exit {done.returncode}  ctx {' '.join(args)}")
    if done.returncode:
        sys.exit(done.stderr[-3000:])


def synthetic(root: Path) -> Path:
    """The made-up gather outputs, gene types and gene sets; return the dataset file."""
    import polars as pl  # noqa: PLC0415
    from synth import write_gene_sets  # noqa: PLC0415
    from test_qc_run import plant, with_annotation  # noqa: PLC0415

    from ctx.dataset import load  # noqa: PLC0415
    from ctx.report.text import word  # noqa: PLC0415

    plant(root)
    spec = with_annotation(root, root)
    sheet = next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))
    frame = pl.read_parquet(sheet)
    frame.with_columns(batch=pl.col("batch").replace(BATCH)).write_parquet(sheet)
    text = spec.read_text()
    for old, new in BATCH.items():
        text = text.replace(f'batch = ["{old}"]', f'batch = ["{new}"]')
    label = word("ui", "preview.dataset")
    spec.write_text(text.replace("[dataset]\n", f'[dataset]\nlabel = "{label}"\n', 1))
    write_gene_sets(root / load(spec).resource("sets").dest)
    return spec


# Read length per made-up batch, so gather derives the batch as it does for a real study.
LENGTHS = {BATCH["short"]: 76, BATCH["long"]: 101}
FIXTURE = REPO / "tests" / "fixtures" / "fetchngs" / "metadata_only"
SDRF_COLUMNS = {
    f"Characteristics[{header}]": column
    for header, column in (("individual", "donor"), ("cell type", "cell_type"),
                           ("condition", "condition"), ("time point", "time_point"),
                           ("sex", "sex"), ("donor sex", "donor_sex"), ("tissue", "tissue"))
}  # fmt: skip


def _rows(path: Path, delimiter: str) -> tuple[list[str], dict[str, str]]:
    """A fixture table's header and first row, the template of every made-up row."""
    import csv  # noqa: PLC0415

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or []), next(reader)


def _write(path: Path, header: list[str], rows: list[dict[str, str]], delimiter: str) -> None:
    """One table in the fixture's own format."""
    import csv  # noqa: PLC0415

    path.parent.mkdir(parents=True, exist_ok=True)
    quoting = csv.QUOTE_ALL if delimiter == "," else csv.QUOTE_MINIMAL
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, header, delimiter=delimiter, quoting=quoting,
                                lineterminator="\n")  # fmt: skip
        writer.writeheader()
        writer.writerows(rows)


def made_up_download(root: Path, spec: Path) -> Path:
    """A metadata-only fetchngs outdir and an SDRF for the made-up samples; extend the spec.

    Written from the planted sheet, so gather curates exactly the samples, groups and sexes the
    counts carry. Accessions are made up (XRR, XRX, XRS, XRP); read counts are the counted reads
    over a made-up assigned share, drawn per run from a seeded generator.
    """
    import hashlib  # noqa: PLC0415
    import json  # noqa: PLC0415
    import shutil  # noqa: PLC0415

    import numpy as np  # noqa: PLC0415
    import polars as pl  # noqa: PLC0415

    gathered = root / "outputs" / "01_gather"
    sheet = pl.read_parquet(next(gathered.glob("*_samplesheet_v1.parquet")))
    counts = pl.read_parquet(next(gathered.glob("*_counts_v1.parquet")))
    counted = {c: int(counts[c].sum()) for c in counts.columns[1:]}
    assigned = np.random.default_rng(7).uniform(0.7, 0.9, sheet.height)
    outdir = root / "data" / "raw" / DATASET
    head_meta, meta = _rows(next((FIXTURE / "metadata").glob("*.tsv")), "\t")
    head_sheet, template = _rows(FIXTURE / "samplesheet" / "samplesheet.csv", ",")
    head_ids, id_row = _rows(FIXTURE / "samplesheet" / "id_mappings.csv", ",")
    sheet_rows, id_rows, sdrf = [], [], []
    for i, r in enumerate(sheet.iter_rows(named=True)):
        run, n = r["run_accession"], i + 1
        reads = int(counted[r["sample"]] / assigned[i])
        length = LENGTHS[r["batch"]]
        fields = {
            "run_accession": run, "experiment_accession": f"XRX{n:06d}",
            "sample_accession": f"XRS{n:06d}", "secondary_sample_accession": f"XRS{n:06d}",
            "study_accession": "XRP000001", "secondary_study_accession": "XRP000001",
            "submission_accession": "XRA000001", "run_alias": run, "experiment_alias": run,
            "sample_alias": r["sample"], "study_alias": "made-up", "library_layout": "SINGLE",
            "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC",
            "library_selection": "cDNA", "library_name": r["sample"],
            "instrument_model": "Illumina NovaSeq 6000", "instrument_platform": "ILLUMINA",
            "base_count": str(reads * length), "read_count": str(reads), "tax_id": "9606",
            "scientific_name": "Homo sapiens", "sample_title": r["sample"],
            "experiment_title": "made-up run", "study_title": "made-up study",
            "sample_description": "made up by scripts/python/build_preview.py", "fastq_md5": "",
            "fastq_bytes": "", "fastq_ftp": "", "fastq_galaxy": "", "fastq_aspera": "",
            "fastq_1": "", "fastq_2": "", "md5_1": "", "md5_2": "", "single_end": "true",
        }  # fmt: skip
        row = {**meta, **{k: v for k, v in fields.items() if k in head_meta}, "id": run}
        _write(outdir / "metadata" / f"{run}.runinfo_ftp.tsv", head_meta, [row], "\t")
        sheet_rows.append({**template, **{k: v for k, v in fields.items() if k in head_sheet},
                           "sample": r["sample"], "fastq_1": f"made-up/{run}.fastq.gz",
                           "fastq_2": "", "strandedness": "auto"})  # fmt: skip
        id_rows.append({**id_row, **{k: v for k, v in fields.items() if k in head_ids},
                        "sample": r["sample"]})  # fmt: skip
        sdrf.append({"Source Name": r["sample"], "Comment[ENA_RUN]": run,
                     **{k: str(r.get(c) or "") for k, c in SDRF_COLUMNS.items()}})  # fmt: skip
    _write(outdir / "samplesheet" / "samplesheet.csv", head_sheet, sheet_rows, ",")
    _write(outdir / "samplesheet" / "id_mappings.csv", head_ids, id_rows, ",")
    (outdir / "ids.csv").write_text("\n".join(sheet["run_accession"]) + "\n")
    for name in (
        "samplesheet/multiqc_config.yml",
        "pipeline_info/nf_core_fetchngs_software_mqc_versions.yml",
    ):
        (outdir / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(FIXTURE / name, outdir / name)
    (params,) = (FIXTURE / "pipeline_info").glob("params_*.json")
    recorded = json.loads(params.read_text()) | {
        "input": str(outdir / "ids.csv"),
        "outdir": str(outdir),
    }
    (outdir / "pipeline_info" / params.name).write_text(json.dumps(recorded, indent=4))
    sdrf_path = root / "resources" / "datasets" / DATASET / "made_up.sdrf.txt"
    header = ["Source Name", *SDRF_COLUMNS, "Comment[ENA_RUN]"]
    _write(sdrf_path, header, sdrf, "\t")
    digest = hashlib.sha256(sdrf_path.read_bytes()).hexdigest()
    columns = "\n".join(f'[columns.{c}]\nsource = "{k}"\nmeaning = "Made up."\n'
                        for k, c in SDRF_COLUMNS.items())  # fmt: skip
    bins = ", ".join(f'"{b}" = [{n}]' for b, n in LENGTHS.items())
    spec.write_text(
        spec.read_text()
        + f"""
[[resources]]
name = "sdrf"
urls = ["http://127.0.0.1/made_up.sdrf.txt"]
bytes = {sdrf_path.stat().st_size}
sha256 = "{digest}"
dest = "{sdrf_path.relative_to(root).as_posix()}"

[sdrf]
resource = "sdrf"
join = ["Comment[ENA_RUN]", "run_accession"]
cross_check = [["Source Name", "sample_title"]]

{columns}
[derived.batch]
bins = {{ {bins} }}
"""
    )
    return outdir


def main() -> None:
    """Build every page, or only the preview modules named."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--only", nargs="*", choices=PREVIEWS, help="preview modules only")
    args = parser.parse_args()
    root, day = args.root.expanduser().resolve(), date.today()
    if args.only is None:
        if root.exists():
            sys.exit(f"{root} exists; move it aside first")
        root.mkdir(parents=True)
        spec = synthetic(root)
        outdir = made_up_download(root, spec)
        ctx(root, "setup")
        ctx(
            root,
            "gather",
            "--dataset",
            DATASET,
            "--fetchngs-dir",
            str(outdir),
            "--dataset-file",
            str(spec),
        )
        for command in ("qc", "de", "function"):
            ctx(root, command, "--dataset", DATASET, "--dataset-file", str(spec))
    for module in args.only or PREVIEWS:
        start = time.perf_counter()
        page = importlib.import_module(f"preview_{module}").page(root, DATASET, day)
        print(f"{time.perf_counter() - start:6.1f} s  preview {module}: {page.name}")
    if args.only is None:
        ctx(root, "report", "--dataset", DATASET, "--dataset-file", str(spec))


if __name__ == "__main__":
    main()
