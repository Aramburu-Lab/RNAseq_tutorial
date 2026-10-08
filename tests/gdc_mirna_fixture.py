"""A made-up GDC project's miRNA and isoform folders, their file rows real, for test_ncrna_tcga.

The four file pairs in tests/fixtures/gdc_mirna are GDC's own files of one TCGA project (Data
Release 46.0), cut to seven hairpins: let-7a's three, miR-21, miR-141, miR-200c and miR-210;
nothing else in them was changed, so each isoform file's rows still sum, hairpin by hairpin,
to its miRNA file's counts. The folders around them are made up as tests/gdc_fixture.py makes
the STAR folder: FAKE barcodes, so the cohort guard (D-020) passes, and GDC's field names.
"""

import hashlib
import json
from collections.abc import Callable, Sequence
from pathlib import Path

from gdc_fixture import METADATA_COLUMNS, PROJECT, RELEASE, Aliquot, _tsv
from gdc_fixture import spec_file as star_spec

from ctx.gather import gdc, gdc_mirna

FIXTURE = Path(__file__).parent / "fixtures" / "gdc_mirna"
WORKFLOW = "BCGSC miRNA Profiling"
FILES = ("a1", "a2", "a3", "a4")
SUFFIX = {gdc_mirna.MIRNA: "mirbase21.mirnas.quantification.txt",
          gdc_mirna.ISOFORM: "mirbase21.isoforms.quantification.txt"}  # fmt: skip


def mirna_aliquot(item: Aliquot) -> Aliquot:
    """The small-RNA library of a STAR aliquot's sample: same sample barcode, another aliquot."""
    parts = item.sample.split("-")
    return Aliquot("-".join([*parts[:4], "21R", "M001", "13"]), item.sample_type,
                   item.preservation)  # fmt: skip


def _query(kind: str, project: str, workflow: str) -> str:
    """A fetch query naming the project, the workflow and the kind's data type."""
    data_type = gdc_mirna.KINDS[kind][2]
    content = [{"op": "in", "content": {"field": f, "value": [v]}} for f, v in (
        ("cases.project.project_id", project), ("analysis.workflow_type", workflow),
        ("data_type", data_type))]  # fmt: skip
    return json.dumps({"filters": {"op": "and", "content": content}})


def build(  # noqa: PLR0913 - every part a test may change
    path: Path,
    kind: str,
    aliquots: Sequence[Aliquot],
    files: Sequence[str] = FILES,
    *,
    workflow: str = WORKFLOW,
    edit: Callable[[str, str], str] | None = None,
) -> Path:
    """Write one kind's folder: one fixture file per aliquot, every record consistent.

    Args:
        path: The folder, created.
        kind: ``gdc_mirna.MIRNA`` or ``gdc_mirna.ISOFORM``.
        aliquots: The libraries, one per file in ``files`` order.
        files: Which fixture files, a1 to a4, each aliquot gets.
        workflow: The workflow the query names.
        edit: Changes a file's text before it is written, given (fixture name, text); the
            manifest's md5 is taken after it, so only the content check can catch the change.

    Returns:
        The folder.
    """
    metadata, manifest, _ = gdc_mirna.KINDS[kind]
    path.mkdir(parents=True)
    meta, log = [], []
    for i, (item, name) in enumerate(zip(aliquots, files, strict=True)):
        file_id = f"{i:08d}-1111-4000-8000-{'0' if kind == gdc_mirna.MIRNA else '1'}" + "0" * 11
        file_name = f"{i:04d}.{SUFFIX[kind]}"
        body = (FIXTURE / f"{name}.{SUFFIX[kind]}").read_text()
        body = edit(name, body) if edit else body
        data = body.encode()
        (path / file_id).mkdir()
        (path / file_id / file_name).write_bytes(data)
        meta.append({
            "file_id": file_id, "file_name": file_name,
            "md5sum": hashlib.md5(data).hexdigest(),  # noqa: S324 - GDC's checksum
            "file_size": str(len(data)), "cases.0.case_id": f"uuid-{item.case}",
            "cases.0.submitter_id": item.case, "cases.0.samples.0.submitter_id": item.sample[:16],
            "cases.0.samples.0.sample_type": item.sample_type,
            "cases.0.samples.0.preservation_method": item.preservation,
            "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id": item.sample,
            "id": file_id,
        })  # fmt: skip
        log.append(f"{file_id}/{file_name}")
    (path / metadata).write_text(_tsv(meta, METADATA_COLUMNS, True), newline="")
    rows = [{"id": m["file_id"], "filename": m["file_name"], "md5": m["md5sum"],
             "size": m["file_size"], "state": "released"} for m in meta]  # fmt: skip
    (path / manifest).write_text(_tsv(rows, ("id", "filename", "md5", "size", "state"), True),
                                 newline="")  # fmt: skip
    (path / gdc.STATUS).write_text(json.dumps({"data_release": RELEASE, "status": "OK"}))
    (path / f"{metadata}{gdc.QUERY}").write_text(_query(kind, PROJECT, workflow))
    retrieval = [{"retrieved_utc": "2099-01-01T00:00:00Z", "gdc_release": RELEASE, "file": f,
                  "source": "made up"} for f in (gdc.STATUS, manifest, metadata, *log)]  # fmt: skip
    (path / gdc.RETRIEVAL).write_text(_tsv(retrieval, ("retrieved_utc", "gdc_release", "file",
                                                       "source"), False))  # fmt: skip
    return path


def spec_file(path: Path) -> Path:
    """The made-up project's dataset file, its [gdc] naming the miRNA workflow."""
    star_spec(path)
    line = 'normal = "Solid Tissue Normal"\n'
    text = path.read_text().replace(line, f'{line}mirna_workflow = "{WORKFLOW}"\n')
    path.write_text(text)
    return path


def rows(name: str, kind: str) -> list[list[str]]:
    """A fixture file's rows, header first, split on tabs."""
    lines = (FIXTURE / f"{name}.{SUFFIX[kind]}").read_text().splitlines()
    return [line.split("\t") for line in lines]
