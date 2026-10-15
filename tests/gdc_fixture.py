"""A GDC project folder shaped as the fetch writes one, made up, for tests/test_gdc.py and M10.

Barcodes start with FAKE, not a cohort's prefix, so the cohort guard (D-020) passes; the reader
takes nothing from barcode shape but the TCGA parts it parses where they fit. Nothing here comes
from a real dataset; the layout and field names are GDC's, as the BRCA folder holds them.
"""

import csv
import hashlib
import io
import json
import random
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ctx.gather import gdc

PROJECT = "FAKE-BRCA"
RELEASE = "Data Release 99.0 - January 1, 2099"
WORKFLOW = "STAR - Counts"
MODEL = "GENCODE v36"
HEADER = ("gene_id", "gene_name", "gene_type", "unstranded", "stranded_first", "stranded_second",
          "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded")  # fmt: skip
# Twelve genes, one _PAR_Y copy of the last X gene; GENCODE-shaped versioned IDs.
GENES = tuple(
    (f"ENSG{n:011d}.{v}", f"G{n}", "protein_coding")
    for n, v in zip(range(1, 13), (3, 1, 7, 2, 5, 1, 4, 9, 2, 6, 1, 3), strict=True)
)
PAR_Y = (f"{GENES[-1][0]}_PAR_Y", GENES[-1][1], GENES[-1][2])
METADATA_COLUMNS = (*gdc.FIELDS.values(), "id")
SPEC = f"""[dataset]
slug = "fake_brca"
label = "Made-up breast"

[gdc]
project = "{PROJECT}"
release = "{RELEASE}"
workflow = "{WORKFLOW}"
gene_model = "{MODEL}"
protocol = "{{protocol}}"
tumour = "Primary Tumor"
normal = "Solid Tissue Normal"
{{cdr}}
[guard]
ids = ["{PROJECT}"]
patterns = []
"""


@dataclass(frozen=True)
class Aliquot:
    """One counts file: its aliquot barcode, sample type, preservation and strand.

    ``strand`` is unstranded, forward or reverse, or a (first, second) pair of strand sums,
    which the file's last gene then carries alone.
    """

    sample: str
    sample_type: str = "Primary Tumor"
    preservation: str = "OCT"
    strand: str | tuple[int, int] = "unstranded"

    @property
    def case(self) -> str:
        """The case barcode the aliquot's begins with."""
        return "-".join(self.sample.split("-")[:3])


@dataclass
class Folder:
    """What ``build`` wrote, for tests that damage one part of it."""

    path: Path
    aliquots: tuple[Aliquot, ...]
    counts: dict[str, list[tuple[str, int, int, int]]] = field(default_factory=dict)


TYPES = {"01": "Primary Tumor", "11": "Solid Tissue Normal", "06": "Metastatic"}


def aliquot(
    case: int, code: str = "01", vial: str = "A", plate: str = "A001", **kw: Any
) -> Aliquot:
    """An aliquot of made-up case ``FAKE-A1-<case>``, its sample type that of its code."""
    kw.setdefault("sample_type", TYPES[code])
    return Aliquot(f"FAKE-A1-{case:04d}-{code}{vial}-11R-{plate}-07", **kw)


def _rows(item: Aliquot, rng: random.Random) -> list[tuple[str, int, int, int]]:
    """(gene, unstranded, first, second) per gene row, _PAR_Y row included, zero."""
    rows = []
    for gene_id, _, _ in GENES:
        total = rng.randint(50, 5000)
        if item.strand == "unstranded" or isinstance(item.strand, tuple):
            first = total // 2 + rng.randint(-total // 50, total // 50)
            second = total - first
        elif item.strand == "reverse":
            # A stranded column holds fewer reads than unstranded: a read on two genes of
            # opposite strands is ambiguous unstranded and counted stranded, and the reverse.
            first, second = total // 25, total - total // 10
        else:
            first, second = total - total // 10, total // 25
        rows.append((gene_id, total, first, second))
    if isinstance(item.strand, tuple):
        rows = [(g, u, 0, 0) for g, u, _, _ in rows[:-1]] + [
            (rows[-1][0], rows[-1][1], *item.strand)
        ]
    rows.insert(len(rows), (PAR_Y[0], 0, 0, 0))
    return rows


def counts_text(rows: Sequence[tuple[str, int, int, int]], model: str = MODEL) -> str:
    """A STAR - Counts file: the gene-model line, the header, four N_ rows, then the genes."""
    names = {g: (n, t) for g, n, t in (*GENES, PAR_Y)}
    lines = [f"# gene-model: {model}", "\t".join(HEADER)]
    lines += [
        f"{n}\t\t\t{v}\t{v}\t{v}\t\t\t" for n, v in zip(gdc.N_ROWS, (11, 22, 33, 44), strict=True)
    ]
    for gene_id, u, f, s in rows:
        name, kind = names[gene_id]
        lines.append(f"{gene_id}\t{name}\t{kind}\t{u}\t{f}\t{s}\t1.0\t1.0\t1.0")
    return "\n".join(lines) + "\n"


def _tsv(rows: Iterable[Mapping[str, Any]], columns: Sequence[str], crlf: bool) -> str:
    """Tab-separated text with a header, lines ended as GDC's API ends them."""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(columns), delimiter="\t",
                            lineterminator="\r\n" if crlf else "\n")  # fmt: skip
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def _query(project: str, workflow: str | None = None) -> str:
    """A fetch query file naming the project and, for files, the workflow."""
    content = [{"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}}]
    if workflow:
        content.append({"op": "in", "content": {"field": "analysis.workflow_type",
                                                "value": [workflow]}})  # fmt: skip
    return json.dumps({"filters": {"op": "and", "content": content}})


def case_record(case: str, **extra: Any) -> dict[str, Any]:
    """One GDC case as the cases endpoint returns it with the fetch's expand fields."""
    record = {
        "submitter_id": case, "case_id": f"uuid-{case}", "primary_site": "Breast",
        "disease_type": "Ductal and Lobular Neoplasms",
        "demographic": {"sex_at_birth": extra.pop("sex", "female"), "age_at_index": 60,
                        "vital_status": "Alive"},
        "diagnoses": [{"diagnosis_is_primary_disease": True,
                       "primary_diagnosis": "Infiltrating duct carcinoma, NOS",
                       "treatments": extra.pop("treatments", [])}],
        "follow_ups": [{"molecular_tests": extra.pop("tests", [])}],
    }  # fmt: skip
    return record | extra


def build(  # noqa: PLR0913 - every part of a folder a test may change
    path: Path,
    aliquots: Sequence[Aliquot],
    *,
    seed: int = 0,
    crlf: bool = True,
    shuffle: bool = True,
    annotations: Sequence[Mapping[str, Any]] = (),
    cases: Mapping[str, Mapping[str, Any]] | None = None,
    release: str = RELEASE,
    project: str = PROJECT,
) -> Folder:
    """Write a GDC folder holding one counts file per aliquot, every record consistent.

    Args:
        path: The folder, created.
        aliquots: The counts files.
        seed: Seeds the counts.
        crlf: End the metadata and manifest lines with CRLF, as GDC's API does.
        shuffle: Put the metadata columns in a seeded order, not GDC's.
        annotations: Annotation hits as the annotations endpoint returns them.
        cases: Extra fields per case for ``case_record``; every case gets a record.
        release: The release the status file and the log name.
        project: The project the file and case queries name.

    Returns:
        What was written.
    """
    rng = random.Random(seed)
    path.mkdir(parents=True)
    folder = Folder(path, tuple(aliquots))
    meta, log = [], []
    for i, item in enumerate(aliquots):
        file_id, name = (
            f"{i:08d}-0000-4000-8000-000000000000",
            f"{i:04d}.rna_seq.augmented_star_gene_counts.tsv",
        )
        rows = _rows(item, rng)
        folder.counts[item.sample] = rows
        body = counts_text(rows).encode()
        (path / file_id).mkdir()
        (path / file_id / name).write_bytes(body)
        values = {
            "file_id": file_id, "file_name": name, "md5sum": hashlib.md5(body).hexdigest(),  # noqa: S324 - GDC's checksum
            "file_size": str(len(body)), "cases.0.case_id": f"uuid-{item.case}",
            "cases.0.submitter_id": item.case, "cases.0.samples.0.submitter_id": item.sample[:16],
            "cases.0.samples.0.sample_type": item.sample_type,
            "cases.0.samples.0.preservation_method": item.preservation,
            "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id": item.sample,
            "id": file_id,
        }  # fmt: skip
        meta.append(values)
        log.append(f"{file_id}/{name}")
    columns = list(METADATA_COLUMNS)
    if shuffle:
        random.Random(seed + 1).shuffle(columns)
    (path / gdc.METADATA).write_text(_tsv(meta, columns, crlf), newline="")
    manifest = [{"id": m["file_id"], "filename": m["file_name"], "md5": m["md5sum"],
                 "size": m["file_size"], "state": "released"} for m in meta]  # fmt: skip
    (path / gdc.MANIFEST).write_text(_tsv(manifest, ("id", "filename", "md5", "size", "state"),
                                          crlf), newline="")  # fmt: skip
    status = {"data_release": release, "status": "OK"}
    (path / gdc.STATUS).write_text(json.dumps(status))
    (path / f"{gdc.METADATA}{gdc.QUERY}").write_text(_query(project, WORKFLOW))
    (path / f"{gdc.CASES}{gdc.QUERY}").write_text(_query(project))
    every = sorted({a.case for a in aliquots})
    extra = cases or {}
    hits = [case_record(c, **dict(extra.get(c, {}))) for c in every]
    (path / gdc.CASES).write_text(json.dumps({"data": {"hits": hits}, "warnings": {}}))
    (path / gdc.ANNOTATIONS).write_text(json.dumps({"data": {"hits": list(annotations)}}))
    logged = [gdc.STATUS, gdc.MANIFEST, gdc.METADATA, gdc.CASES, gdc.ANNOTATIONS, *log]
    retrieval = [{"retrieved_utc": "2099-01-01T00:00:00Z", "gdc_release": release, "file": f,
                  "source": "made up"} for f in logged]  # fmt: skip
    (path / gdc.RETRIEVAL).write_text(_tsv(retrieval, ("retrieved_utc", "gdc_release", "file",
                                                       "source"), False))  # fmt: skip
    return folder


def spec_file(path: Path, protocol: str = "unstranded", cdr: str | None = None) -> Path:
    """A dataset file for the made-up project; ``cdr`` names a CDR workbook in the folder."""
    lines = f'cdr_file = "{cdr}"\ncdr_type = "BRCA"\n' if cdr else ""
    path.write_text(SPEC.format(protocol=protocol, cdr=lines))
    return path


def write_cdr(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    """A TCGA-CDR workbook: the unnamed row-number column, then the columns, in a TCGA-CDR sheet.

    Raises:
        ImportError: Without openpyxl, which writes it.
    """
    import openpyxl  # noqa: PLC0415 - only a CDR test needs it

    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = "TCGA-CDR"
    sheet.append([None, *columns])
    for i, row in enumerate(rows, 1):
        sheet.append([i, *(row.get(c) for c in columns)])
    book.create_sheet("Notes").append(["Recommended use of the endpoints:"])
    book.save(path)
