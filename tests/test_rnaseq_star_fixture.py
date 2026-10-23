"""The committed aligned nf-core/rnaseq fixture is real, pinned, and holds the copy list only.

Every file M2's read stage and M3 read from an outdir aligned with STAR (D-116), and nothing
else: no BAM, no HTML report, no file of a folder the pages never open.
"""

import hashlib
import json
import re
import subprocess
from pathlib import Path

from ctx.gather.rnaseq import ALIGNED_FILES, COPIED_FILES, RNASEQ_VERSION, VERSIONS_FILE

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "rnaseq_star"
# Measured 2026-10-23: 416,986 bytes; the cap leaves room for one more table, not for a BAM.
BUDGET = 600_000
SAMPLES = {"WT_REP1", "WT_REP2", "RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2",
           "RAP1_IAA_30M_REP1"}  # fmt: skip


def _data() -> set[Path]:
    """Every file of the fixture but its README and manifest."""
    return {p for p in FIXTURE.rglob("*") if p.is_file()} - {
        FIXTURE / "README.md", FIXTURE / "MANIFEST.sha256"}  # fmt: skip


def test_manifest_and_budget():
    lines = (FIXTURE / "MANIFEST.sha256").read_text().splitlines()
    entries = {name: digest for digest, name in (line.split("  ", 1) for line in lines)}
    on_disk = sorted(p.relative_to(FIXTURE).as_posix() for p in FIXTURE.rglob("*")
                     if p.is_file() and p.name != "MANIFEST.sha256")  # fmt: skip
    assert entries and sorted(entries) == on_disk
    assert {n: hashlib.sha256((FIXTURE / n).read_bytes()).hexdigest() for n in entries} == entries
    assert 0 < sum(p.stat().st_size for p in FIXTURE.rglob("*") if p.is_file()) <= BUDGET
    files = [p.relative_to(ROOT).as_posix() for p in FIXTURE.rglob("*") if p.is_file()]
    tracked = subprocess.run(["git", "ls-files", "--", *files], cwd=ROOT, check=True,
                             capture_output=True, text=True).stdout.split()  # fmt: skip
    assert sorted(tracked) == sorted(files)


def test_holds_the_copy_list_and_nothing_else():
    """Every pattern of COPIED_FILES and ALIGNED_FILES matches, and no other file is here."""
    patterns = (*COPIED_FILES, *ALIGNED_FILES)
    assert all(any(FIXTURE.glob(p)) for p in patterns), [p for p in patterns
                                                         if not any(FIXTURE.glob(p))]  # fmt: skip
    assert _data() == {p for pattern in patterns for p in FIXTURE.glob(pattern)}
    assert not [p for p in _data() if p.suffix in {".bam", ".bai", ".html"}]


def test_is_the_pinned_release_aligned_with_star():
    assert re.search(rf"nf-core/rnaseq: v{re.escape(RNASEQ_VERSION)}\b",
                     (FIXTURE / VERSIONS_FILE).read_text())  # fmt: skip
    for path in sorted((FIXTURE / "pipeline_info").glob("params_*.json")):
        params = json.loads(path.read_text())
        assert (params["aligner"], params["pseudo_aligner"], params["skip_alignment"]) == (
            "star_salmon", "salmon", False), path.name  # fmt: skip


def test_both_matrices_hold_the_test_samples():
    for relative in ("salmon/salmon.merged.gene_counts_length_scaled.tsv",
                     "star_salmon/salmon.merged.gene_counts_length_scaled.tsv"):  # fmt: skip
        header = (FIXTURE / relative).read_text().splitlines()[0].split("\t")
        assert header[:2] == ["gene_id", "gene_name"] and set(header[2:]) == SAMPLES, relative
