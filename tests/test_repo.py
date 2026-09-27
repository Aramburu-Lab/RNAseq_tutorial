"""Repository layout after adopting codebase_template (D-002), and dataset independence."""

import csv
import io
import re
import subprocess
import tomllib
import urllib.request
from pathlib import Path

import pytest

from ctx.gather.fetchngs import FETCHNGS_VERSION
from ctx.gather.rnaseq import RNASEQ_VERSION

ROOT = Path(__file__).resolve().parents[1]
PRECEDENCE = (
    "Where `CONVENTIONS.md` conflicts with `AGENTS.md`, `CONVENTIONS.md` wins; "
    "each divergence is a D-entry in `knowledge/decisions.md`."
)
RULE_FILES = ("CONVENTIONS.md", "AGENTS.md", "CLAUDE.md")
COHORT_ID = re.compile(r"E-MTAB-\d+|GSE\d+|TCGA-[A-Z]+|PXD\d+|PDC\d+")
# D20: nf-core/fetchngs's own test list holds GEO IDs; this one directory is exempt, and
# tests/test_fetchngs_fixture.py pins its input list to nf-core's byte for byte.
NF_CORE_FIXTURE = Path("tests", "fixtures", "fetchngs", "metadata_only")
DATASETS = Path("conf", "datasets")


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


def guard_patterns(root):
    """COHORT_ID plus every id and pattern a dataset file declares in its [guard] table."""
    patterns = [COHORT_ID]
    for path in sorted((root / DATASETS).glob("*.toml")):
        declared = tomllib.loads(path.read_text()).get("guard")
        if not declared:
            raise ValueError(f"{path.name}: no [guard] table; a dataset file must declare its ids")
        patterns += [re.compile(rf"\b{re.escape(i)}\b") for i in declared["ids"]]
        patterns += [re.compile(p) for p in declared["patterns"]]
    return patterns


def cohort_hits(root):
    """Cohort identifiers outside the dataset files, as 'path: token'."""
    patterns = guard_patterns(root)
    files = [
        p
        for top in ("src", "tests", "scripts")
        for p in (root / top).rglob("*")
        if p.is_file()
        and "__pycache__" not in p.parts
        and not p.is_relative_to(root / NF_CORE_FIXTURE)
    ] + [root / "pyproject.toml", root / "pixi.toml"]
    return [
        f"{p.relative_to(root)}: {m.group()}"
        for p in files
        if p.is_file()
        for pattern in patterns
        for m in pattern.finditer(p.read_text(errors="ignore"))
    ]


def test_no_cohort_identifiers():
    assert len(list((ROOT / DATASETS).glob("*.toml"))) > 0
    hits = cohort_hits(ROOT)
    assert not hits, "cohort identifiers belong in conf/datasets/*.toml only"


def test_guard_reads_dataset_files(tmp_path):
    (tmp_path / DATASETS).mkdir(parents=True)
    (tmp_path / DATASETS / "x.toml").write_text(
        "[guard]\nids = [\"PRJXX1\"]\npatterns = ['\\bZZ\\d{3}\\b']\n"
    )
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "a.py").write_text("PRJXX1 ZZ123 ZZ1234\n")
    assert cohort_hits(tmp_path) == ["src/a.py: PRJXX1", "src/a.py: ZZ123"]
    (tmp_path / DATASETS / "y.toml").write_text('[dataset]\nslug = "y"\n')
    with pytest.raises(ValueError, match=r"y.toml: no \[guard\]"):
        cohort_hits(tmp_path)


@pytest.mark.network
def test_guard_patterns_cover_study():
    fields = "run_accession,experiment_accession,sample_accession,secondary_sample_accession"
    for path in sorted((ROOT / DATASETS).glob("*.toml")):
        declared = tomllib.loads(path.read_text())
        patterns = guard_patterns(ROOT)[1:]
        url = (
            "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&format=tsv"
            f"&accession={declared['dataset']['ena_study']}&fields={fields}"
        )
        with urllib.request.urlopen(url, timeout=120) as response:
            rows = list(csv.DictReader(io.StringIO(response.read().decode()), delimiter="\t"))
        values = {r[f] for r in rows for f in fields.split(",")}
        sdrf = [r for r in declared.get("resources", []) if r["name"] == "sdrf"]
        for resource in sdrf:
            with urllib.request.urlopen(resource["urls"][0], timeout=120) as response:
                table = list(csv.reader(io.StringIO(response.read().decode()), delimiter="\t"))
            column = table[0].index("Characteristics[individual]")
            values |= {row[column] for row in table[1:]}
        assert rows and len(values) > len(rows)
        uncovered = sorted(v for v in values if not any(p.search(v) for p in patterns))
        assert not uncovered, f"{path.name}: {len(uncovered)} not guarded, e.g. {uncovered[:5]}"
    fixtures = [p for p in (ROOT / "tests" / "fixtures").rglob("*") if p.is_file()]
    found = [
        f"{p.relative_to(ROOT)}: {m.group()}"
        for p in fixtures
        for pattern in guard_patterns(ROOT)[1:]
        for m in pattern.finditer(p.read_text(errors="ignore"))
    ]
    assert not found, found


@pytest.mark.parametrize(
    ("pipeline", "pinned"), [("fetchngs", FETCHNGS_VERSION), ("rnaseq", RNASEQ_VERSION)]
)
def test_version_single_source(pipeline, pinned):
    sources = [ROOT / "TODO.md", ROOT / "CONVENTIONS.md"] + [
        p for top in ("scripts", "src") for p in (ROOT / top).rglob("*.py")
    ]
    found = [
        (p.name, v)
        for p in sources
        for v in re.findall(rf"nf-core/{pipeline} -r (\d+\.\d+\.\d+)", p.read_text())
    ]
    assert found, f"no pinned {pipeline} command found"
    assert all(v == pinned for _, v in found), found
