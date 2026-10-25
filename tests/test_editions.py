"""ED.1 and ED.2: editions, conf/editions/<name>.toml, and the global --edition (ctx.editions).

Each test names the mutation it catches. Ten were planted and each caught (D-110), among them
a caller bypassing the rule (the hub, START.html, the kit's folder filter, the CLI, the guide,
ctx run) and the full edition losing a module (``clinical`` out of ``cohorts.RUN_MODULES``,
which full's ``"all"`` reads, or a list in full.toml without it).
"""

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest
from conftest import ctx_run
from test_kit import NEW, SCRIPT, a_root, nothing_written, page

from ctx import cohorts, editions, pipeline, runtime
from ctx.dataset import load
from ctx.io.readme_companion import write_table
from ctx.report import course, guide, hub, text
from ctx.report.builder import Group, unpacked
from ctx.report.entry import onward

ROOT = Path(__file__).resolve().parents[1]
TODAY = "2026-10-22"
DATASETS = ROOT / "conf" / "datasets"
# The course's ten modules, typed here so a module lost from the code or from full.toml fails.
EVERY = ("gather", "qc", "process", "quant", "de", "function", "proteo", "ncrna", "genecard",
         "clinical")  # fmt: skip
CORE = ("gather", "qc", "process", "quant", "de", "function", "ncrna", "genecard")
CORE_DATASETS = ("Non-cancerous-blood", "Non-cancerous-blood-four-runs", "Breast-tissue-BRCA",
                 "BRCA-immune-cells-macrophages", "BRCA-immune-cells-monocytes", "Lung-tissue-LUAD",
                 "Pancreas-tissue-PAAD")  # fmt: skip
# Every course dataset file the core edition leaves out, on 2026-10-22. A dataset file added to
# conf/datasets lands here or in core.toml.
FULL_ONLY = {"cptac_brca", "cptac_ccrcc", "r2_gse130605", "r3_pxd038437", "tcga_brca",
             "tcga_hnsc", "tcga_kirc", "tcga_luad", "tcga_lung", "tcga_lusc", "u_gse327166",
             "u_gse327906"}  # fmt: skip
# Hub, START.html and teacher titles of the core modules, M0 first.
CORE_NUMBERS = ["M0", "M1", "M2", "M3", "M4", "M5", "M6", "M8", "M9"]


@pytest.fixture
def core(monkeypatch):
    monkeypatch.setenv(editions.EDITION_ENV, "core")
    return editions.active()


@pytest.fixture
def full(monkeypatch):
    monkeypatch.delenv(editions.EDITION_ENV, raising=False)
    return editions.active()


def test_the_two_editions(full):
    """full is the default, every module and every dataset file; core is Rytis's 2026-10-22 list.
    Catches full losing a module (in full.toml or in cohorts.RUN_MODULES), the default moved off
    full, and core gaining M7 or M10."""
    assert full.name == "full" and editions.DEFAULT == "full"
    assert full.modules == EVERY == cohorts.RUN_MODULES
    assert full.datasets is None
    assert full.dataset_files() == tuple(sorted(DATASETS.glob("*.toml")))
    core = editions.named("core")
    assert core.modules == CORE and core.datasets == CORE_DATASETS
    assert not core.has_module("proteo") and not core.has_module("clinical")
    assert all(core.has_module(m) for m in ("setup", "report", *CORE))


def test_every_course_dataset_is_core_or_named_full_only():
    """A dataset file added to conf/datasets must be placed: listed in core.toml, or named in
    FULL_ONLY here. Catches the ED.4 lung file landing while core still refuses its slug."""
    core = editions.named("core")
    outside = {p.stem for p in DATASETS.glob("*.toml")} - set(core.datasets or ())
    assert outside == FULL_ONLY, (
        f"{sorted(outside ^ FULL_ONLY)}: add a core dataset's slug to datasets in "
        "conf/editions/core.toml, or a full-only one to FULL_ONLY in this test"
    )


BAD = """
[edition]
name = "Other"
label = ""
colour = "red"
modules = ["qc", "de", "de", "plot"]
datasets = ["Non-cancerous-blood", "absent", "renamed"]

[extra]
"""


def test_a_bad_edition_lists_every_problem(tmp_path):
    """One read names each problem. Catches a loader that stops at the first, and each check:
    unknown table and key, a name that is not the file's one-word stem, a blank label and notes,
    an unknown or doubled module, gather missing, a slug with no file or another slug."""
    folder, datasets = tmp_path / "editions", tmp_path / "datasets"
    folder.mkdir()
    datasets.mkdir()
    (datasets / "Non-cancerous-blood.toml").write_text('[dataset]\nslug = "Non-cancerous-blood"\n')
    (datasets / "renamed.toml").write_text('[dataset]\nslug = "other"\n')
    (folder / "bad.toml").write_text(BAD)
    with pytest.raises(editions.EditionError) as caught:
        editions.load(folder / "bad.toml")
    said = str(caught.value)
    for problem in ("unknown table [extra]", "edition.colour: unknown key",
                    "edition.name: 'Other'", "edition.label: expected non-blank",
                    "edition.notes: expected non-blank", "'plot' is not one of",
                    "a name listed twice", "gather is missing", "no dataset file",
                    "renamed.toml declares slug 'other'"):  # fmt: skip
        assert problem in said, problem
    (folder / "order.toml").write_text(
        '[edition]\nname = "order"\nlabel = "L"\nnotes = "N"\nmodules = ["de", "gather"]\n'
        'datasets = "all"\n'
    )
    with pytest.raises(editions.EditionError, match="not in course order"):
        editions.load(folder / "order.toml")


def test_an_unknown_edition_is_refused_and_no_folder_means_full(tmp_path, monkeypatch):
    """A name no file has is refused with the names there are; an installed package away from
    the repository (no editions folder) runs full and nothing else. Catches a typo read as full."""
    with pytest.raises(editions.EditionError, match="known: core, full"):
        editions.named("cor")
    monkeypatch.setenv(editions.FOLDER_ENV, str(tmp_path / "none"))
    assert editions.named("full").modules == EVERY
    with pytest.raises(editions.EditionError):
        editions.named("core")


@pytest.mark.parametrize("module", ["clinical", "proteo"])
def test_core_refuses_m7_and_m10_by_name_and_writes_nothing(tmp_path, module):
    """`ctx --edition core clinical|proteo` exits 1 naming the module, and so does CTX_EDITION;
    nothing is written. Catches the CLI running a module outside the edition."""
    spec = str(DATASETS / "Non-cancerous-blood.toml")
    args = (module, "--dataset", "Non-cancerous-blood", "--dataset-file", spec)
    for flag, env in ((("--edition", "core"), None), ((), {editions.EDITION_ENV: "core"})):
        done = ctx_run("--root", str(tmp_path), *flag, *args, cwd=tmp_path, env=env)
        assert done.returncode == 1
        assert f"ctx {module} is not in the core edition" in done.stderr
        assert not (tmp_path / "outputs").exists()


def test_core_refuses_a_course_dataset_outside_it_and_keeps_the_readers_own(core, tmp_path):
    """A course dataset core does not list is refused by slug or by its dataset file; a slug with
    no course file is the reader's own and runs. Catches the own-data routes closed by core."""
    said = core.refusal("de", "tcga_brca")
    assert said and "tcga_brca" in said and str(core.path) in said
    assert core.refusal("de", core.course_slug(DATASETS / "tcga_kirc.toml"))
    assert core.course_slug(tmp_path / "tcga_kirc.toml") is None
    assert core.refusal("gather", "my_lab_study") is None
    assert core.refusal("de", "Non-cancerous-blood") is None
    done = ctx_run("--root", str(tmp_path), "--edition", "core", "guide", "--dataset-file",
                   str(DATASETS / "tcga_brca.toml"), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 1 and "outside the core edition" in done.stderr
    done = ctx_run("--root", str(tmp_path), "--edition", "core", "figure", "breast_m10_km",
                   "--params", str(tmp_path / "p.toml"), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 1 and "ctx clinical is not in the core edition" in done.stderr


def all_module_pages(root: Path, slug: str) -> Path:
    """One overview page of every module, and M0's, for one dataset."""
    out = root / "outputs"
    page(out / "00_setup" / "2026-10-22_setup_environment_v1.html", "M0 Setup")
    for module, number in [(m, n) for m, n in hub.MODULES if m not in hub.NOT_MODULES][1:]:
        page(out / f"{number:02d}_{module}" / f"2026-10-22_{module}_{slug}_overview_v1.html",
             f"{hub.TITLES[module]}: Study")  # fmt: skip
    return root


def _hub_modules(root: Path) -> list[str]:
    """The modules the hub links, in its order, from the folders its links name."""
    settings = runtime.Settings("laptop", root, 1, None, 0.0)
    html = unpacked(hub.run(settings, (), "study")[0].read_text())
    return re.findall(r'class="module-name" href="\.\./\d{2}_([a-z_]+)/', html)


def test_the_hub_lists_exactly_the_editions_modules(tmp_path, monkeypatch):
    """With a page of every module under the root, core's hub lists M0 to M6, M8 and M9, its
    lite page with no key numbers (D-111); full's lists M0 to M10. Mutation caught: hub.run
    without has_module."""
    root = all_module_pages(tmp_path, "study")
    monkeypatch.setenv(editions.EDITION_ENV, "core")
    assert _hub_modules(root) == ["setup", *CORE]
    hub_page = next((root / "outputs" / "09_report").glob("*_hub_v1.html")).read_text()
    assert "07_proteo" not in hub_page and "12_clinical" not in hub_page
    assert 'class="keys"' not in unpacked(hub_page)
    monkeypatch.delenv(editions.EDITION_ENV)
    assert _hub_modules(root) == ["setup", *EVERY]


def test_start_and_the_teacher_guide_list_exactly_the_editions_modules(tmp_path, core):
    """START.html's rows and the teacher guide's table under core are its modules. Catches
    course.page or teacher.gathered listing every module."""
    from ctx.report import teacher  # noqa: PLC0415

    root = all_module_pages(tmp_path, "study")
    home = course.page(tmp_path, [(root, True)])
    assert [r.number for r in home.rows] == CORE_NUMBERS
    assert [k.value for k in home.key_numbers] == ["9 of 9", "8 of 8"]
    assert [m for m, _ in teacher.gathered([root])] == ["setup", *CORE]


def _part_titles(groups: tuple[Group, ...]) -> list[str]:
    return [g.title for g in groups]


def test_the_guide_shows_the_m7_and_tcga_parts_only_where_the_edition_holds_them(monkeypatch):
    """Under core the guide keeps every own-data part (1b, 1c) and drops the M7 part and the
    TCGA part, and its method stops naming TCGA; full keeps both, numbered as before. Catches
    either part shown outside its edition, and full renumbered."""
    spec = load(guide.DEFAULT)
    monkeypatch.delenv(editions.EDITION_ENV, raising=False)
    full = _part_titles(guide.blocks(spec, guide.DEFAULT))
    protein = text.fill("guide", "mass_spec.group", n=4)
    tcga = text.fill("guide", "gdc.group", n=5)
    assert full[-2:] == [protein, tcga] and guide.has_tcga()
    monkeypatch.setenv(editions.EDITION_ENV, "core")
    core = _part_titles(guide.blocks(spec, guide.DEFAULT))
    assert core == full[:-2]
    assert text.fill("guide", "part1b.group", n=1) in core
    assert text.fill("guide", "part1c.group", n=1) in core
    assert not guide.has_tcga()
    assert "TCGA" not in text.raw("guide", "variants.no_tcga.method")


def test_ctx_run_and_the_gene_list_links_leave_out_what_the_edition_leaves_out(core, tmp_path):
    """ctx run under core plans no M10 and names why; M5's onward links skip M10; figure kinds
    offered read the same rule. Catches pipeline.plan, route_modules or entry.onward bypassing
    the edition."""
    spec = load(DATASETS / "tcga_brca.toml")
    options = pipeline.Options("tcga_brca", spec.path, gdc_dirs=(tmp_path,))
    planned, left_out = pipeline.plan(tmp_path, spec, options)
    assert "clinical" not in [m for m, _ in planned]
    assert text.fill("run", "outside_edition", module="clinical", edition="core") in left_out
    assert "clinical" not in [m for m, _ in onward(spec)]
    assert pipeline.route_modules(spec) & {"clinical", "proteo"} == set()


def test_full_plans_and_links_as_before(full, tmp_path):
    """The same calls under full keep M10. Catches the edition rule leaking into full."""
    spec = load(DATASETS / "tcga_brca.toml")
    options = pipeline.Options("tcga_brca", spec.path, gdc_dirs=(tmp_path,))
    assert "clinical" in [m for m, _ in pipeline.plan(tmp_path, spec, options)[0]]
    assert "clinical" in [m for m, _ in onward(spec)]
    assert {"clinical", "proteo"} <= pipeline.route_modules(spec)


def _kit(tmp: Path, name: str, roots: list[str], *extra: str
         ) -> tuple[Path, subprocess.CompletedProcess]:  # fmt: skip
    out = tmp / name
    env = {k: v for k, v in os.environ.items() if k != editions.EDITION_ENV}
    done = subprocess.run([sys.executable, str(SCRIPT), "--out", str(out), *roots, *extra],
                          capture_output=True, text=True, check=False, env=env)  # fmt: skip
    return out, done


def _every_module_root(tmp: Path, edition: str | None = None) -> Path:
    """A root with an M7 and an M10 page beside the course's, none linked from the hub.

    ``edition`` stamps every page as built under it; None leaves them as pages built before
    ED.2b, recording none.
    """
    root = a_root(tmp / "real", "Breast tissue", "Breast-tissue-BRCA")
    out = root / "outputs"
    page(out / "07_proteo" / "2026-09-30_proteo_Breast-tissue-BRCA_overview_v1.html", "M7: Breast")
    page(out / "12_clinical" / "2026-09-30_clinical_Breast-tissue-BRCA_overview_v1.html", "M10: B")
    # M10 reads a gene list from a link, so a kit that ships its folder needs its entry page.
    page(out / "12_clinical" / "index.html", "M10 Clinical and survival")
    if edition is not None:
        _stamp(out.rglob("*.html"), edition)
    return root


def _stamp(pages, edition: str) -> None:
    """Record ``edition`` on each page as the builder does, on its ``<html>`` tag."""
    for p in pages:
        html = p.read_text().replace("<html>", f'<html {editions.PAGE_ATTRIBUTE}="{edition}">', 1)
        p.write_text(html)


def test_a_core_kit_ships_no_m7_or_m10_page_and_records_its_edition(tmp_path):
    """build_kit --edition core leaves the M7 and M10 folders out, lists only core rows on
    START.html, and states core in MANIFEST.json and README.md. Catches select() shipping a
    folder outside the edition, and an unrecorded edition."""
    root = _every_module_root(tmp_path, "core")
    kit, done = _kit(tmp_path, "kit", [f"--real=step={root}"], "--edition", "core")
    assert done.returncode == 0, done.stderr
    shipped = {p.relative_to(kit).as_posix() for p in kit.rglob("*") if p.is_file()}
    assert not [p for p in shipped if "/07_proteo/" in p or "/12_clinical/" in p]
    assert json.loads((kit / "MANIFEST.json").read_text())["edition"] == "core"
    core = editions.named("core")
    line = text.fill("course", "readme.edition", label=core.label, name="core",
                     modules=", ".join(core.modules))  # fmt: skip
    assert line in (kit / "README.md").read_text()
    home = unpacked((kit / course.INDEX).read_text())
    assert hub.TITLES["proteo"].partition(" ")[2] not in home
    assert hub.TITLES["clinical"].partition(" ")[2] not in home
    assert hub.TITLES["ncrna"].partition(" ")[2] in home


def test_a_core_kit_refuses_a_root_outside_it_and_a_link_to_m7(tmp_path):
    """A root holding a course dataset core leaves out, or a shipped page linking an M7 page,
    stops the build with nothing written. Catches build() and select() bypassing the edition."""
    outside = a_root(tmp_path / "kidney", "Kidney", "tcga_kirc")
    kit, done = _kit(tmp_path, "kit", [f"--real=kidney={outside}"], "--edition", "core")
    assert done.returncode != 0 and "tcga_kirc is a course dataset outside" in done.stderr
    assert not kit.exists() and not kit.with_name("kit.partial").exists()
    root = _every_module_root(tmp_path / "linked")
    hub_page = next((root / "outputs" / "09_report").glob("*_hub_v1.html"))
    page(hub_page, "Course hub: Breast tissue",
         "../07_proteo/2026-09-30_proteo_Breast-tissue-BRCA_overview_v1.html")  # fmt: skip
    _stamp((root / "outputs").rglob("*.html"), "core")
    kit, done = _kit(tmp_path, "kit2", [f"--real=step={root}"], "--edition", "core")
    assert done.returncode != 0 and "outside the core edition" in done.stderr
    assert not kit.exists()


def _hashes(kit: Path) -> dict[str, str]:
    return {p.relative_to(kit).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in kit.rglob("*") if p.is_file()}  # fmt: skip


def test_the_full_kit_is_the_default_kit_byte_for_byte(tmp_path):
    """build_kit with no --edition and with --edition full write the same files, byte for byte,
    M7 and M10 included, and record full. Catches a default other than full, and full losing a
    module folder."""
    root = _every_module_root(tmp_path)
    plain, done = _kit(tmp_path, "plain", [f"--real=step={root}"])
    assert done.returncode == 0, done.stderr
    named, done = _kit(tmp_path, "named", [f"--real=step={root}"], "--edition", "full")
    assert done.returncode == 0, done.stderr
    assert _hashes(plain) == _hashes(named)
    assert any("/07_proteo/" in p for p in _hashes(plain))
    assert any("/12_clinical/" in p for p in _hashes(plain))
    assert json.loads((plain / "MANIFEST.json").read_text())["edition"] == "full"


def test_a_page_records_the_edition_it_was_built_under(tmp_path, monkeypatch):
    """A hub written under core records core on its <html> tag, and so do the glossary it
    writes and a table's README; under full, full. Catches builder.render not recording the
    edition, or recording a fixed one."""
    root = all_module_pages(tmp_path, "study")
    settings = runtime.Settings("laptop", root, 1, None, 0.0)
    for name in ("core", "full"):
        monkeypatch.setenv(editions.EDITION_ENV, name)
        html = hub.run(settings, (), "study")[0].read_text()
        assert editions.recorded(html) == name
        assert f'<html lang="en" {editions.PAGE_ATTRIBUTE}="{name}">' in html
        (glossary_page,) = (root / "outputs" / "09_report").glob("*_report_glossary_v1.html")
        assert editions.recorded(glossary_page.read_text()) == name
        table = write_table(pl.DataFrame({"n": [1]}), tmp_path / name, f"{TODAY}_de_x_v1",
                            columns={"n": ("count", "One.")}, data_releases={})  # fmt: skip
        assert f"| Edition | {name} |" in table[-1].read_text()
    assert editions.recorded("<html><head></head></html>") is None


def test_a_core_kit_refuses_a_page_built_under_full_naming_it(tmp_path, monkeypatch):
    """A core root whose hub was written again under full, by the builder, stops a core kit
    with nothing written, naming the hub and its edition; a page recording none is refused for
    core, with the command that rebuilds it. Mutation caught: build_kit not calling
    other_edition (the kit not reading the edition), which shipped the full hub."""
    root = _every_module_root(tmp_path, "core")
    monkeypatch.delenv(editions.EDITION_ENV, raising=False)
    settings = runtime.Settings("laptop", root, 1, None, 0.0)
    hub_page = hub.run(settings, (), "Breast-tissue-BRCA")[0]
    kit, done = _kit(tmp_path, "kit", [f"--real=step={root}"], "--edition", "core")
    assert done.returncode != 0 and nothing_written(kit)
    assert f"{hub_page.resolve()}: built under the full edition, not core" in done.stderr
    old = _every_module_root(tmp_path / "old")
    kit, done = _kit(tmp_path, "kit2", [f"--real=step={old}"], "--edition", "core")
    assert done.returncode != 0 and nothing_written(kit)
    gather = next((old / "outputs" / "01_gather").glob(f"{NEW}_*overview_v1.html")).resolve()
    assert f"{gather}: records no edition" in done.stderr
    assert "rebuild it under core (ctx --edition core)" in done.stderr


def test_a_full_kit_of_pages_built_before_ed2b_builds_and_ships_them_unchanged(tmp_path):
    """Pages recording no edition are full's: a full kit of them builds, and ships each page
    byte for byte. Catches the edition check refusing old full kits."""
    root = _every_module_root(tmp_path)
    kit, done = _kit(tmp_path, "kit", [f"--real=step={root}"])
    assert done.returncode == 0, done.stderr
    shipped = sorted((kit / "step" / "outputs").rglob("*.html"))
    assert shipped and all(editions.recorded(p.read_text()) is None for p in shipped)
    for p in shipped:
        assert p.read_bytes() == (root / "outputs" / p.relative_to(kit / "step" / "outputs")
                                  ).read_bytes()  # fmt: skip


def test_the_core_hub_keeps_a_left_out_module_in_its_place(tmp_path, monkeypatch):
    """The core hub went from M6 to M8 (Rytis 2026-10-25). Now M7 and M10 stand in their places,
    muted and linking nowhere. Catches the rows dropped, a row that links, and full showing
    them. The core hub, a lite page, shows no run count to read (D-111)."""
    root = all_module_pages(tmp_path, "study")
    settings = runtime.Settings("laptop", root, 1, None, 0.0)
    monkeypatch.setenv(editions.EDITION_ENV, "core")
    html = unpacked(hub.run(settings, (), "study")[0].read_text())
    kickers = re.findall(r'class="kicker">(M\d+)<', html)
    assert kickers == ["M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10"]
    out = re.findall(
        r'<li class="left-out"><span class="kicker">(M\d+)</span><span class="module-name">', html
    )
    assert out == ["M7", "M10"]
    monkeypatch.delenv(editions.EDITION_ENV)
    full = unpacked(hub.run(settings, (), "study")[0].read_text())
    assert 'class="left-out"' not in full


def test_an_empty_lite_page_is_left_out_of_the_hub_and_start(tmp_path, monkeypatch):
    """A core module page with nothing but a check table is marked ``data-empty``; the hub and
    START.html leave it out, and under full the same page is listed. Catches hub.run or
    course.page linking an empty page, and a page with a figure marked empty."""
    from test_lite import CHECK, every_part  # noqa: PLC0415
    from test_report import _module  # noqa: PLC0415

    from ctx.report.builder import ResultTable, render, shown  # noqa: PLC0415

    root = all_module_pages(tmp_path, "study")
    blank = _module(results=(ResultTable(CHECK, ("Run", "Reason"), (("r1", "low"),)),))
    setup = root / "outputs" / "00_setup" / "2026-10-22_setup_environment_v1.html"
    monkeypatch.setenv(editions.EDITION_ENV, "core")
    setup.write_text(render(blank))
    kept = tmp_path / "kept.html"
    kept.write_text(render(every_part()))
    assert not shown(setup) and shown(kept)
    assert _hub_modules(root) == list(CORE)
    assert [r.number for r in course.page(tmp_path, [(root, True)]).rows] == CORE_NUMBERS[1:]
    monkeypatch.delenv(editions.EDITION_ENV)
    setup.write_text(render(blank))
    assert shown(setup) and _hub_modules(root) == ["setup", *EVERY]
