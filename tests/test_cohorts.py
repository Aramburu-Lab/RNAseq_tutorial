"""652: cohort bundles, conf/cohorts/<cancer>.toml, and the global --cohort (ctx.cohorts).

Each test names its mutation.
"""

import tomllib
from pathlib import Path

import pytest
from conftest import ctx_run
from test_own_table import FIXTURE, dataset, table

from ctx import cohorts
from ctx.gather.fetchngs import read

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "conf" / "datasets"


def test_every_bundle_loads_and_each_dataset_file_has_one_cohort():
    """The bundles in the repository load, and every dataset file is listed by exactly one.
    Catches a dataset file added with no cohort, and one listed twice, which would make its
    figures' IDs depend on which bundle is read first."""
    found = cohorts.bundles()
    # ovary and pleura hold the R2.2 proofs (wo13-u), merged after the bundles were written;
    # myeloma the R3 proof (R3.3); pancreas the core edition's pancreas study (ED.8).
    assert set(found) == {"myeloid", "breast", "kidney", "lung", "headneck", "ovary", "pleura",
                          "myeloma", "pancreas"}  # fmt: skip
    listed = [e.slug for c in found.values() for e in c.datasets]
    assert sorted(listed) == sorted(p.stem for p in DATASETS.glob("*.toml"))
    assert len(listed) == len(set(listed))
    for cohort in found.values():
        for entry in cohort.datasets:
            assert cohorts.cohort_of(entry.slug) == cohort.name


BAD = """
[cohort]
name = "Other"
label = ""
colour = "red"

[extra]

[[datasets]]
slug = "tier1"
step = 4
modules = ["de", "de", "plot"]
laptop_cap = 0
note = " "
size = 3

[[datasets]]
slug = "tier1"
step = 1
modules = ["de"]
laptop_cap = "all"
note = "Again."

[[datasets]]
slug = "renamed"
step = 1
modules = ["de"]
laptop_cap = 150
note = "Its file declares another slug."

[[datasets]]
slug = "absent"
step = 3
modules = []
laptop_cap = "all"
note = "No file."
"""


def test_a_bad_bundle_lists_every_problem(tmp_path):
    """One read names each problem of a bundle. Catches a loader that stops at the first, and
    each check: unknown table and key, a name that is not the file's one-word stem, a blank
    label and note, a Step out of range, an unknown or doubled module, a cap that is neither
    all nor positive, a slug listed twice, one whose file declares another slug, one with none."""
    folder = tmp_path / "cohorts"
    folder.mkdir()
    (tmp_path / "datasets").mkdir()
    (tmp_path / "datasets" / "tier1.toml").write_bytes((DATASETS / "tier1.toml").read_bytes())
    (tmp_path / "datasets" / "renamed.toml").write_text('[dataset]\nslug = "else"\n')
    (folder / "other.toml").write_text(BAD)
    with pytest.raises(cohorts.CohortError) as error:
        cohorts.load(folder / "other.toml")
    said = str(error.value)
    for part in ("unknown table [extra]", "cohort.colour: unknown key", "cohort.name: 'Other'",
                 "cohort.label: expected non-blank", "datasets[0].step: expected one of",
                 "'plot' is not one of", "a module listed twice", "datasets[0].laptop_cap",
                 "datasets[0].note", "datasets[0].size: unknown key",
                 "datasets: 'tier1' listed twice", "renamed.toml declares slug 'else'",
                 "datasets[3].slug: no dataset file", "datasets[3].modules"):  # fmt: skip
        assert part in said, (part, said)
    assert "datasets[2].laptop_cap" not in said


def test_cohort_of_takes_the_choice_then_the_listing_then_own(tmp_path, monkeypatch):
    """A chosen cohort names every figure; else the bundle listing the slug; else own; two
    bundles listing it with none chosen is refused. Catches the listing preferred over the
    choice, and a dataset no bundle lists left without an ID."""
    monkeypatch.delenv(cohorts.COHORT_ENV, raising=False)
    assert cohorts.cohort_of("tcga_kirc") == "kidney"
    assert cohorts.cohort_of("mine") == cohorts.OWN
    monkeypatch.setenv(cohorts.COHORT_ENV, "breast")
    assert cohorts.cohort_of("tcga_kirc") == "breast"
    monkeypatch.delenv(cohorts.COHORT_ENV)
    folder = tmp_path / "cohorts"
    folder.mkdir()
    for name in ("one", "two"):
        (folder / f"{name}.toml").write_text('[[datasets]]\nslug = "shared"\n')
    monkeypatch.setenv(cohorts.FOLDER_ENV, str(folder))
    with pytest.raises(cohorts.CohortError, match="listed by one, two"):
        cohorts.cohort_of("shared")


def test_cohort_option_names_figures_and_refuses_an_unknown_cohort(tmp_path):
    """--cohort myeloid puts myeloid in the figure IDs gather records; a cohort no bundle has stops
    the command before anything is written, naming the known ones. Catches the option read
    but not passed to the figures, and an unknown cohort silently drawn as own."""
    root, conf = tmp_path / "root", tmp_path / "conf"
    root.mkdir()
    spec = dataset(conf, "mine.csv", table(read(FIXTURE).runs))
    args = ("gather", "--dataset", "mine", "--dataset-file", str(spec), "--fetchngs-dir",
            str(FIXTURE))  # fmt: skip
    done = ctx_run("--root", str(root), "--cohort", "nope", *args, cwd=root)
    assert done.returncode == 1 and "no cohort 'nope'" in done.stderr, done.stderr
    assert "myeloid" in done.stderr and not (root / "outputs").exists()
    done = ctx_run("--root", str(root), "--cohort", "myeloid", *args, cwd=root)
    assert done.returncode == 0, done.stderr[-3000:]
    (readme,) = (root / "outputs" / "01_gather").glob("*_mine_composition_v1_README.md")
    assert "`myeloid_m1_composition`" in readme.read_text()


def test_bundles_declare_no_cap_the_modules_do_not_apply():
    """Every cap is all: no module subsamples on the laptop profile yet, so a number here
    would state a cap that does not happen. Catches a planned cap written as if applied."""
    for path in sorted((ROOT / "conf" / "cohorts").glob("*.toml")):
        for entry in tomllib.loads(path.read_text())["datasets"]:
            assert entry["laptop_cap"] == cohorts.ALL, (path.name, entry["slug"])
