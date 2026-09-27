"""The committed nf-core/rnaseq fixture is real, pinned, complete, and holds only COPIED_FILES."""

import csv
import hashlib
import re
import subprocess
from pathlib import Path

from ctx.gather.rnaseq import COPIED_FILES, RNASEQ_VERSION, VERSIONS_FILE

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "rnaseq"
BUDGET = 2_000_000
# Sample names of nf-core's test samplesheet at test-datasets 626c8fab (samplesheet/v3.10).
TEST_SAMPLES = {
    "WT_REP1",
    "WT_REP2",
    "RAP1_UNINDUCED_REP1",
    "RAP1_UNINDUCED_REP2",
    "RAP1_IAA_30M_REP1",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_manifest():
    lines = (FIXTURE / "MANIFEST.sha256").read_text().splitlines()
    entries = {name: digest for digest, name in (line.split("  ", 1) for line in lines)}
    on_disk = sorted(
        p.relative_to(FIXTURE).as_posix()
        for p in FIXTURE.rglob("*")
        if p.is_file() and p.name != "MANIFEST.sha256"
    )
    assert entries and sorted(entries) == on_disk
    assert {name: sha256(FIXTURE / name) for name in entries} == entries
    assert 0 < sum(p.stat().st_size for p in FIXTURE.rglob("*") if p.is_file()) <= BUDGET
    files = [p.relative_to(ROOT).as_posix() for p in FIXTURE.rglob("*") if p.is_file()]
    tracked = subprocess.run(
        ["git", "ls-files", "--", *files], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.split()
    assert sorted(tracked) == sorted(files)


def test_holds_copied_files_only():
    data = {p for p in FIXTURE.rglob("*") if p.is_file()} - {
        FIXTURE / "README.md",
        FIXTURE / "MANIFEST.sha256",
    }
    matched = {p for pattern in COPIED_FILES for p in FIXTURE.glob(pattern)}
    assert all(any(FIXTURE.glob(pattern)) for pattern in COPIED_FILES)
    assert data == matched


def test_is_pinned_release():
    text = (FIXTURE / VERSIONS_FILE).read_text()
    assert re.search(rf"nf-core/rnaseq: v{re.escape(RNASEQ_VERSION)}\b", text)
    assert re.search(r"^\s+salmon: 2\.7\.0$", text, re.MULTILINE)


def test_layout():
    path = FIXTURE / "salmon" / "salmon.merged.gene_counts_length_scaled.tsv"
    with path.open(newline="") as handle:
        header, *rows = list(csv.reader(handle, delimiter="\t"))
    assert header[:2] == ["gene_id", "gene_name"]
    assert len(header[2:]) == len(set(header[2:])) and set(header[2:]) == TEST_SAMPLES
    values = [float(v) for row in rows for v in row[2:]]
    assert rows and min(values) >= 0
    assert any(v != int(v) for v in values), "length-scaled counts should include fractions"
