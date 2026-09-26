"""The committed nf-core/fetchngs fixtures are real, pinned, complete, and exempt only as stated.

Consistency between files of one run proves nothing about whether the run happened (a copy
made by hand could agree with itself), so the proof of reality is external: the network tests
against ENA and a rebuild, plus the README provenance.
"""

import csv
import hashlib
import importlib.util
import json
import re
import subprocess
import urllib.request
from pathlib import Path

import pytest

from ctx.gather.fetchngs import FETCHNGS_VERSION

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "fetchngs"
MODES = ("download", "metadata_only")
BUDGET = 2_000_000
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&fields=run_accession,fastq_md5,fastq_bytes&format=json"


@pytest.fixture(scope="module")
def builder():
    spec = importlib.util.spec_from_file_location(
        "make_fetchngs_fixture", ROOT / "scripts" / "python" / "make_fetchngs_fixture.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest(mode):
    lines = (FIXTURES / mode / "MANIFEST.sha256").read_text().splitlines()
    return {name: digest for digest, name in (line.split("  ", 1) for line in lines)}


def sheet(mode):
    with (FIXTURES / mode / "samplesheet" / "samplesheet.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


@pytest.mark.parametrize("mode", MODES)
def test_manifest(mode):
    root = FIXTURES / mode
    on_disk = sorted(
        p.relative_to(root).as_posix()
        for p in root.rglob("*")
        if p.is_file() and p.name != "MANIFEST.sha256"
    )
    entries = manifest(mode)
    assert len(entries) > 0
    assert sorted(entries) == on_disk
    assert {name: sha256(root / name) for name in entries} == entries
    files = [p.relative_to(ROOT).as_posix() for p in root.rglob("*") if p.is_file()]
    tracked = subprocess.run(
        ["git", "ls-files", "--", *files], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.split()
    assert sorted(tracked) == sorted(files)


def test_budget():
    total = sum(p.stat().st_size for p in FIXTURES.rglob("*") if p.is_file())
    assert 0 < total <= BUDGET


def test_exempt_input_is_nf_core_list(builder):
    # The cohort scan skips metadata_only/ (D20); this is what keeps that exemption narrow.
    assert sha256(FIXTURES / "metadata_only" / "ids.csv") == builder.NF_CORE_LIST_SHA256
    download_ids = (FIXTURES / "download" / "ids.csv").read_text().split()
    assert download_ids == list(builder.DOWNLOAD_IDS)


@pytest.mark.parametrize("mode", MODES)
def test_is_pinned_release(mode):
    text = (
        FIXTURES / mode / "pipeline_info" / "nf_core_fetchngs_software_mqc_versions.yml"
    ).read_text()
    found = re.search(r"nf-core/fetchngs: v(\d+\.\d+\.\d+)-g([0-9a-f]{7})", text)
    assert found and found.group(1) == FETCHNGS_VERSION
    (params,) = (FIXTURES / mode / "pipeline_info").glob("params_*.json")
    recorded = json.loads(params.read_text())
    assert recorded["nf_core_pipeline"] == "rnaseq"
    assert recorded["skip_fastq_download"] is (mode == "metadata_only")


def test_md5_consistency():
    """Consistency only: every md5 file agrees with its FASTQ and with the sheet's ENA md5."""
    root = FIXTURES / "download"
    md5_files = sorted((root / "fastq" / "md5").glob("*.md5"))
    assert md5_files
    by_name = {}
    for row in sheet("download"):
        names = [Path(row[k]).name for k in ("fastq_1", "fastq_2") if row[k]]
        ena = row["fastq_md5"].split(";")[-len(names) :] if row["fastq_md5"] else []
        by_name |= dict(zip(names, ena, strict=False))
    for path in md5_files:
        digest, name = path.read_text().split()
        assert hashlib.md5((root / "fastq" / name).read_bytes()).hexdigest() == digest
        assert by_name[name] == digest


@pytest.mark.network
def test_matches_ena():
    root = FIXTURES / "download"
    for row in sheet("download"):
        with urllib.request.urlopen(ENA.format(row["run_accession"]), timeout=60) as response:
            (ena,) = json.load(response)
        names = [Path(row[k]).name for k in ("fastq_1", "fastq_2") if row[k]]
        if not ena["fastq_md5"]:
            # No ENA FASTQ: fetchngs fell back to sra-tools, which writes no md5 file.
            assert not any((root / "fastq" / "md5" / f"{n}.md5").exists() for n in names)
            continue
        expected = ena["fastq_md5"].split(";")[-len(names) :]
        assert [
            hashlib.md5((root / "fastq" / n).read_bytes()).hexdigest() for n in names
        ] == expected


@pytest.mark.network
@pytest.mark.parametrize("mode", MODES)
def test_rebuilds(builder, tmp_path, mode):
    out = tmp_path / mode
    assert (
        builder.main(["--mode", mode, "--out", str(out), "--scratch", str(tmp_path / "scratch")])
        == 0
    )
    compared = [
        p.relative_to(FIXTURES / mode)
        for p in (FIXTURES / mode).rglob("*")
        if p.is_file()
        and p.parts[len((FIXTURES / mode).parts)] in {"samplesheet", "metadata", "fastq"}
    ]
    assert compared
    assert {p: (out / p).read_bytes() for p in compared} == {
        p: (FIXTURES / mode / p).read_bytes() for p in compared
    }
