"""Repository layout after adopting codebase_template (D-002), and dataset independence."""

import ast
import csv
import fnmatch
import io
import re
import subprocess
import tomllib
import urllib.request
from pathlib import Path

import pytest

from ctx.dataset import load
from ctx.gather.fetchngs import FETCHNGS_VERSION
from ctx.gather.rnaseq import COPIED_FILES, RNASEQ_VERSION

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


def test_subsets_curate_as_their_parent():
    """A dataset declaring ``subset_of`` reads its parent's SDRF the parent's way (D-042).

    Only [dataset] slug, label and subset_of, [expect], and resources it never uses may
    differ; a column, value map, bin or guard edited in one file alone fails here.
    """
    subsets = {
        path: raw
        for path in sorted((ROOT / DATASETS).glob("*.toml"))
        if "subset_of" in (raw := tomllib.loads(path.read_text()))["dataset"]
    }
    assert ROOT / DATASETS / "tier1_teach.toml" in subsets
    own = ("slug", "label", "subset_of")
    for path, sub in subsets.items():
        parent = tomllib.loads(
            (ROOT / DATASETS / f"{sub['dataset']['subset_of']}.toml").read_text()
        )
        assert set(sub) == {"dataset", "resources", "guard", "sdrf", "columns", "order",
                            "derived", "expect"}, path.name  # fmt: skip
        for table in ("guard", "sdrf", "columns", "order", "derived"):
            assert sub[table] == parent[table], (path.name, table)
        named = {r["name"]: r for r in parent["resources"]}
        assert all(r == named.get(r["name"]) for r in sub["resources"]), path.name
        assert {k: v for k, v in sub["dataset"].items() if k not in own} == {
            k: v for k, v in parent["dataset"].items() if k not in own
        }, path.name
        assert load(path).slug == sub["dataset"]["slug"]


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


def test_fetchngs_runs_monochrome():
    """Every place that runs fetchngs, or tells a reader how to, turns colours off (D-006).

    fetchngs 1.12.0's nf-validation 1.1.3 reads the undefined params.monochromeLogs when
    colours are on, and Nextflow warns. A guard, not a proof; what escapes it is named below.
    """
    todo_text = (ROOT / "TODO.md").read_text()
    todo = todo_text.splitlines()
    conventions = (ROOT / "CONVENTIONS.md").read_text().splitlines()
    config = "nfcore_fetchngs/nextflow.config"
    # The M1b course command, inside its backticks: the sentence after it names the flag too.
    (course,) = [line for line in todo if line.startswith("- Course command")]
    (command,) = re.findall(r"`(nextflow run nf-core/fetchngs[^`]*)`", course)
    assert "--monochrome_logs true" in command, command
    (rule,) = [line for line in conventions if line.startswith("- Raw reads enter only through")]
    assert "--monochrome_logs true" in rule and "monochrome_logs = true" in rule, rule
    # Scripts: every file under scripts/ whose code (comment lines removed) names
    # nf-core/fetchngs must carry the flag's tokens in that code. Escapes: a pipeline name
    # built from variables, and the tokens placed in a string outside the command.
    flag = re.compile(r'--monochrome_logs"?[\s,]+"?true')
    code = {}
    for path in (ROOT / "scripts").rglob("*"):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        try:
            lines = path.read_text().splitlines()
        except UnicodeDecodeError:  # not a script: compiled or binary
            continue
        code[path] = "\n".join(s for s in lines if not s.lstrip().startswith("#"))
    scripts = {p: text for p, text in code.items() if "nf-core/fetchngs" in text}
    assert scripts, "no script names nf-core/fetchngs"
    assert all(flag.search(text) for text in scripts.values()), list(scripts)
    # The Kebnekaise launcher takes no extra flags and starts Nextflow in
    # <project_dir>/nfcore_fetchngs, so its runs carry the setting in that nextflow.config.
    # Prose routes: lines naming the launcher with -profile singularity. Escapes: a route
    # written without that profile string.
    launcher = ("SHARED_run_fetchngs.sbatch", "-profile singularity")
    routes = [line for line in [*todo, *conventions] if all(k in line for k in launcher)]
    assert routes, "no prose route through the launcher"
    for line in routes:
        assert "monochrome_logs = true" in line and config in line, line
    # Code: every fenced TODO block that calls the launcher writes that config, by printf or
    # by a heredoc (D-035), with the setting, before the first call. Escape: a write by any
    # other means (cp, tee) that drops the setting.
    call = re.compile(r"SHARED_run_fetchngs\.sbatch\s")
    heredoc = re.compile(
        r"^[^\n]*" + re.escape(config) + r"[^\n]*<<'(\w+)'\n(.*?)^\1$", re.M | re.S
    )
    blocks = todo_text.split("```")[1::2]
    calls = [b for b in blocks if call.search(b)]
    assert calls, "no fenced block calls the launcher"
    for block in calls:
        writes = [
            (block.index(line), line)
            for line in block.splitlines()
            if "printf" in line and config in line
        ]
        writes += [(m.start(), m.group(0)) for m in heredoc.finditer(block)]
        assert writes and all("monochrome_logs = true" in w for _, w in writes), block[:200]
        first_call = call.search(block)
        assert first_call and min(writes)[0] < first_call.start(), "config after call"


def test_rsync_matches_copied_files():
    """Step 7's rsync copies exactly COPIED_FILES (D-022): what leaves the cluster is one list.

    A pattern the rsync lacks is never copied, and ctx fails on the Mac after the cluster run;
    an include with no pattern copies a file nothing reads.
    """
    (command,) = [
        line
        for line in (ROOT / "TODO.md").read_text().splitlines()
        if line.startswith("rsync ") and "tier1_rnaseq" in line
    ]
    includes = [i for i in re.findall(r"--include='([^']+)'", command) if i != "*/"]
    assert "--exclude='*'" in command
    for pattern in COPIED_FILES:
        assert any(i == pattern or fnmatch.fnmatch(pattern, i) for i in includes), pattern
    for include in includes:
        assert any(p == include or fnmatch.fnmatch(p, include) for p in COPIED_FILES), include


def eager_params(config: str) -> set[str]:
    """Params a Nextflow config reads as it is parsed: outside every ``key: { }`` closure."""
    return {
        name
        for line in config.splitlines()
        for name in re.findall(r"params\.(\w+)", re.sub(r"\w+:\s*\{.*\}", "", line))
    }


def test_rnaseq_override_reads_only_given_params():
    """Every param the publish override reads when parsed is set in the run's params (D-025).

    Under NXF_SYNTAX_PARSER=v2 a -c config cannot see a param only the pipeline sets:
    ``mode: params.publish_dir_mode`` came out empty and published symlinks into work/.
    """
    fixture = (ROOT / "scripts" / "python" / "make_rnaseq_fixture.py").read_text()
    names = {node.targets[0].id: node.value for node in ast.parse(fixture).body
             if isinstance(node, ast.Assign) and hasattr(node.targets[0], "id")}  # fmt: skip
    publish, params = ast.literal_eval(names["PUBLISH"]), ast.literal_eval(names["PARAMS"])
    assert eager_params(publish) == {"publish_dir_mode"}
    assert eager_params(publish) <= set(params)
    todo = (ROOT / "TODO.md").read_text()
    (given,) = re.findall(r'cat > "\$P/rnaseq_params.yml" <<EOF\n(.*?)\nEOF', todo, re.S)
    assert eager_params(publish) <= set(re.findall(r"^(\w+):", given, re.M))


def test_publish_override_single_source():
    """The cluster run and the fixture publish the Trim Galore JSON with one override text (D-025).

    Without it nf-core/rnaseq 3.27.0 publishes only the .txt reports, which MultiQC 1.35 skips for
    Trim Galore 2.x, and the read fates of ctx qc would be empty.
    """
    source = (ROOT / "scripts" / "python" / "make_rnaseq_fixture.py").read_text()
    (publish,) = [
        node.value.value
        for node in ast.parse(source).body
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "PUBLISH"
    ]
    todo = (ROOT / "TODO.md").read_text()
    assert f"""cat >> "$P/rnaseq_kebnekaise.config" <<'EOF'\n{publish}EOF""" in todo
    assert '-c "$P/rnaseq_kebnekaise.config"' in todo
