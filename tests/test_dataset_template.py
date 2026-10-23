"""U.2 and U.8: the dataset-file template a student copies, and the guide's worked example.

No accession is typed here: the worked example's come from
conf/datasets/Non-cancerous-blood-four-runs.toml, the fixture's from its own samplesheet. Each
test names the mutation it catches.
"""

import shutil
import tomllib
from pathlib import Path

import polars as pl

from ctx import pipeline
from ctx.dataset import load
from ctx.gather.curation import read_table
from ctx.gather.fetchngs import read
from ctx.gather.harmonise import curated_runs

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "conf" / "dataset_template.toml"
STEP1 = ROOT / "conf" / "datasets" / "Non-cancerous-blood.toml"
FOUR = ROOT / "conf" / "datasets" / "Non-cancerous-blood-four-runs.toml"
FIXTURE = ROOT / "tests" / "fixtures" / "fetchngs" / "metadata_only"


def _resources(path: Path) -> dict[str, dict]:
    return {r["name"]: r for r in tomllib.loads(path.read_text()).get("resources", [])}


def test_template_is_the_one_pipeline_reads():
    """ctx run --check fills its lines from this file. Catches the path drifting from it."""
    assert pipeline.TEMPLATE == TEMPLATE


def test_template_resources_are_step1s():
    """The gene sets and annotation a new study starts with are Step-1's, entry for entry, and
    so are the worked example's. Catches a URL, size or hash edited in one copy alone."""
    step1 = _resources(STEP1)
    for path in (TEMPLATE, FOUR):
        found = _resources(path)
        assert set(found) == {"msigdb_hallmark", "annotation_gtf"}, path.name
        assert all(found[n] == step1[n] for n in found), path.name


def test_template_function_block_is_step1s():
    """The [function] lines --check prints are Step-1's gene sets, licence and featured set.
    Catches the licence or the set edited in the template alone."""
    filled = tomllib.loads(pipeline.block("[function]", contrast="c"))["function"]
    step1 = tomllib.loads(STEP1.read_text())["function"]
    assert {k: v for k, v in filled.items() if k != "running"} == {
        k: v for k, v in step1.items() if k != "running"
    }
    assert filled["running"] == ["c", step1["running"][1]]


def test_filled_template_curates_the_fixture(tmp_path):
    """The template copied as the guide copies it, its slug and table filled in, curates the
    fetchngs fixture's runs. Catches a column source or the join in the template that no
    table started from gather's sheet can satisfy."""
    spec_path = tmp_path / "mine.toml"
    shutil.copy(TEMPLATE, spec_path)
    runs = read(FIXTURE).runs
    table = runs.select("run_accession").with_columns(
        donor=pl.format("d{}", pl.int_range(pl.len())), cell_type=pl.lit("monocyte"),
        condition=pl.when(pl.int_range(pl.len()) < 7).then(pl.lit("control"))
        .otherwise(pl.lit("case")),
    )  # fmt: skip
    table.write_csv(tmp_path / "mine.csv")
    spec = load(spec_path)
    assert (spec.slug, spec.accessions, spec.table()) == ("mine", (), tmp_path / "mine.csv")
    curated = curated_runs(tmp_path / "root", spec, read(FIXTURE), {})
    assert curated.height == runs.height and set(curated["condition"]) == {"control", "case"}


def test_four_run_table_is_its_accessions():
    """The worked example's table lists exactly the runs its guide downloads, and Step-1's
    guard covers it. Catches an accession added to one list alone."""
    spec = load(FOUR)
    body, _ = read_table(spec.table())
    rows = [line.split("\t") for line in body.splitlines()]
    assert rows[0][0] == "run_accession" and len(spec.accessions) == 4
    assert [r[0] for r in rows[1:]] == list(spec.accessions)
    step1 = tomllib.loads(STEP1.read_text())
    assert tomllib.loads(FOUR.read_text())["guard"] == step1["guard"]
    assert spec.contrasts[0].expected_n == (2, 2)
