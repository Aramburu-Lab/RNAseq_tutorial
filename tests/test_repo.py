"""Repository layout after adopting codebase_template (D-002), and dataset independence."""

import re
import subprocess
from pathlib import Path

import pytest

from ctx.gather.fetchngs import FETCHNGS_VERSION

ROOT = Path(__file__).resolve().parents[1]
PRECEDENCE = (
    "Where `CONVENTIONS.md` conflicts with `AGENTS.md`, `CONVENTIONS.md` wins; "
    "each divergence is a D-entry in `knowledge/decisions.md`."
)
RULE_FILES = ("CONVENTIONS.md", "AGENTS.md", "CLAUDE.md")
COHORT_ID = re.compile(r"E-MTAB-\d+|GSE\d+|TCGA-[A-Z]+|PXD\d+|PDC\d+")
# D20: nf-core/fetchngs's own test list holds GEO IDs; this one directory is exempt, and
# tests/test_fetchngs_fixture.py pins its input list to nf-core's byte for byte.
NF_CORE_FIXTURE = ROOT / "tests" / "fixtures" / "fetchngs" / "metadata_only"


def test_layout():
    assert [p.name for p in ROOT.iterdir() if p.name.lower() == "todo.md"] == ["TODO.md"]
    # Order matters: the project rules are imported after the template's, and override them.
    imports = re.findall(r"^@(\S+)$", (ROOT / "CLAUDE.md").read_text(), flags=re.MULTILINE)
    assert imports == ["AGENTS.md", "CONVENTIONS.md"]
    assert all((ROOT / name).is_file() for name in imports)
    for name in RULE_FILES:
        text = (ROOT / name).read_text()
        assert not re.search(r"(?<!knowledge/)lessons\.md", text), f"{name}: bare lessons.md"
        assert "todo.md" not in text, f"{name}: lowercase todo.md"
        assert PRECEDENCE in text, f"{name}: precedence sentence missing"


@pytest.mark.parametrize(
    ("path", "ignored"),
    [
        ("data/x.csv", True),
        ("resources/manifest.json", True),
        ("tests/fixtures/chr22/.DS_Store", True),
        ("tests/fixtures/x.key", True),
        ("results/samplesheet/samplesheet.csv", True),
        ("tests/fixtures/x.fa", False),
        ("tests/fixtures/x.fq", False),
        ("tests/fixtures/x.gz", False),
        ("tests/fixtures/fetchngs/download/fastq/x_1.fastq.gz", False),
        # Proves /data/ is anchored: a package subdirectory called data stays tracked.
        ("src/ctx/data/x.py", False),
    ],
)
def test_gitignore_rules(path, ignored):
    probe = subprocess.run(["git", "check-ignore", "-q", "--no-index", path], cwd=ROOT, check=False)
    assert probe.returncode == (0 if ignored else 1)


def test_no_cohort_identifiers():
    files = [
        p
        for top in ("src", "tests", "scripts")
        for p in (ROOT / top).rglob("*")
        if p.is_file() and "__pycache__" not in p.parts and not p.is_relative_to(NF_CORE_FIXTURE)
    ] + [ROOT / "pyproject.toml", ROOT / "pixi.toml"]
    assert len(files) > 2
    hits = [
        f"{p.relative_to(ROOT)}: {m.group()}"
        for p in files
        for m in COHORT_ID.finditer(p.read_text(errors="ignore"))
    ]
    assert not hits, "Phase 1 must not depend on a dataset"


def test_fetchngs_version_single_source():
    sources = [ROOT / "TODO.md", ROOT / "CONVENTIONS.md"] + [
        p for top in ("scripts", "src") for p in (ROOT / top).rglob("*.py")
    ]
    pinned = [
        (p.name, v)
        for p in sources
        for v in re.findall(r"nf-core/fetchngs -r (\d+\.\d+\.\d+)", p.read_text())
    ]
    assert pinned, "no pinned fetchngs command found"
    assert all(v == FETCHNGS_VERSION for _, v in pinned), pinned
