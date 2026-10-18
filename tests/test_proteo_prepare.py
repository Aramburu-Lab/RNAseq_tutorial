"""R3: the PRIDE and starter-table helper, the template's [proteo] block, and the guide's M7 part.

The block is pasted as the guide says, into the template copied as the guide copies it, and
``ctx proteo`` runs on it. The PRIDE calls go to a made-up opener. Each test names the mutation
it catches. No real accession is typed here.
"""

import csv
import hashlib
import io
import json
import shlex
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pytest
from clinical_root import ctx
from proteo_study_root import SAMPLES, gather, maxquant_text, protein
from test_guide import parse

from ctx import pipeline
from ctx.dataset import load
from ctx.proteo import guide, prepare
from ctx.report import guide as page

ROOT = Path(__file__).resolve().parents[1]
FOUR = ROOT / "conf" / "datasets" / "four_runs.toml"
BODY = b"made-up protein table\n"


class Opener:
    """A made-up network: PRIDE's listing as JSON, and the file's bytes."""

    def __init__(self, body: bytes, sha1: str):
        self.body, self.sha1, self.calls = body, sha1, []

    def __call__(self, url: str, timeout: int):
        self.calls.append(url)
        if "/ws/archive/" in url:
            entry = {
                "fileName": "proteinGroups.txt",
                "fileSizeBytes": len(BODY),
                "checksum": self.sha1,
                "publicFileLocations": [{"value": "ftp://ftp.pride.ebi.ac.uk/a/proteinGroups.txt"}],
            }
            return io.BytesIO(json.dumps([entry]).encode())
        return io.BytesIO(self.body)


def test_fetch_checks_size_and_sha1(tmp_path):
    """Mutations caught: the FTP URL kept where HTTPS is served; a body of the wrong SHA-1 kept;
    a partial file left behind; the SHA-256 in the entry not the file's."""
    good = Opener(BODY, hashlib.sha1(BODY).hexdigest())  # noqa: S324 - PRIDE's own checksum
    (entry,) = prepare.listing("PXDTEST", good)
    assert entry["url"] == "https://ftp.pride.ebi.ac.uk/a/proteinGroups.txt"
    bad = Opener(b"made-up protein tablX\n", entry["sha1"])
    with pytest.raises(OSError, match="SHA-1"):
        prepare.fetch(entry, tmp_path, bad)
    assert list(tmp_path.iterdir()) == []
    path = prepare.fetch(entry, tmp_path, good)
    block = prepare.resource_block("protein_table", entry["url"], path, "data/raw/x/p.txt")
    assert (
        f'sha256 = "{hashlib.sha256(BODY).hexdigest()}"' in block
        and f"bytes = {len(BODY)}" in block
    )


def test_template_block_runs_ctx_proteo(tmp_path):
    """The guide's route on a made-up study: the starter sample table from the export, the
    template's [proteo] block pasted, ctx proteo. Mutations caught: the block's keys out of
    step with ctx.dataset; the starter missing a sample column; a zip member not found."""
    rng = np.random.default_rng(11)
    root = tmp_path / "root"
    samples, rna = gather(root, rng)
    columns = [f"ms_{j:02d}" for j in range(SAMPLES)]
    archive = root / "data" / "raw" / "fake_study" / "txt.zip"
    archive.parent.mkdir(parents=True)
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(
            "combined/txt/proteinGroups.txt", maxquant_text(protein(rna, rng), columns, rng)
        )
    spec_path = tmp_path / "fake_study.toml"
    shutil.copy(pipeline.TEMPLATE, spec_path)
    table = tmp_path / "fake_study_proteo.csv"
    assert (
        prepare.main(
            [
                "samples",
                "--format",
                "maxquant",
                "--table",
                str(archive),
                "--member",
                "proteinGroups.txt",
                "--out",
                str(table),
            ]
        )
        == 0
    )
    rows = list(csv.reader(table.open()))
    assert rows[0] == list(prepare.COLUMNS) and [r[0] for r in rows[1:]] == columns
    table.write_text(
        "protein_sample,rna_sample,batch\n"
        + "".join(f"{c},{samples[j]},run{j}\n" for j, c in enumerate(columns))
    )
    text = spec_path.read_text().replace('slug = "mine"', 'slug = "fake_study"')
    text += pipeline.block("[proteo]", slug="fake_study")
    text += (
        '\n[[resources]]\nname = "protein_table"\nurls = ["https://example.org/txt.zip"]\n'
        f'bytes = 1\nsha256 = "{"0" * 64}"\ndest = "data/raw/fake_study/txt.zip"\n'
    )
    spec_path.write_text(text)
    spec = load(spec_path)
    gmt = root / spec.resource("msigdb_hallmark").dest
    gmt.parent.mkdir(parents=True, exist_ok=True)
    gmt.write_text("HALLMARK_OXIDATIVE_PHOSPHORYLATION\tx\tG5\tG6\tG7\tG8\n")
    done = ctx(
        "--root", str(root), "proteo", "--dataset", "fake_study", "--dataset-file", str(spec_path)
    )
    assert any(line.endswith("_overview_v1.html") for line in done.stdout.splitlines())


def test_guide_part_lines_parse():
    """Every line of the M7 part is one its own parser takes, and the block is the template's.
    Mutations caught: a prepare option renamed in the parser alone; a ctx option renamed; the
    block typed into the guide instead of read from the template."""
    spec = load(FOUR)
    found = guide.lines(spec, "conf/datasets/four_runs.toml")
    for line in (*found["list"], *found["fetch"], *found["samples"]):
        if line.startswith(guide.PREPARE):
            args = shlex.split(line.removeprefix(guide.PREPARE))
            prepare.parser().parse_args(args)
    assert {parse(line).info_name for line in found["run"][1:]} == {"setup", "proteo"}
    assert "\n".join(found["block"]) + "\n" == pipeline.block("[proteo]", slug=spec.slug)
    groups = page.blocks(spec, FOUR)
    # Part 4, before the TCGA part (Part 5), which the guide ends on.
    assert [g.title for g in groups][-2:] == [
        "Part 4: protein of the same samples, for M7",
        "Part 5: any TCGA project",
    ]
