"""scripts/python/build_kit.py, check_kit.py and ctx.report.course: what a student kit holds."""

import hashlib
import importlib.util
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

from ctx.report import course, text

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "python" / "build_kit.py"
CHECK = ROOT / "scripts" / "python" / "check_kit.py"
OLD, NEW = "2026-09-29", "2026-09-30"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def kit_script():
    return load(SCRIPT, "build_kit")


@pytest.fixture(scope="module")
def checker():
    return load(CHECK, "check_kit")


def page(path: Path, title: str, *hrefs: str, extra: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "".join(f'<a href="{h}">x</a>' for h in hrefs) + extra
    path.write_text(f"<html><head><title>{title}</title></head><body>{body}</body></html>")


def a_root(root: Path, name: str, slug: str) -> Path:
    """A ctx root: an old and a new M1 page, a linked and an unlinked undated file, a hub."""
    out = root / "outputs"
    page(out / "00_setup" / f"{NEW}_setup_environment_v1.html", "M0 Setup")
    page(out / "01_gather" / f"{OLD}_gather_{slug}_overview_v1.html", f"M1 Gather: {name}")
    page(out / "01_gather" / f"{NEW}_gather_{slug}_overview_v1.html", f"M1 Gather: {name}",
         f"{NEW}_gather_{slug}_counts_v1.csv", "../02_qc/multiqc/report.html", "#top")  # fmt: skip
    (out / "01_gather" / f"{NEW}_gather_{slug}_counts_v1.csv").write_text("gene,n\nA,1\n")
    page(out / "02_qc" / "multiqc" / "report.html", "MultiQC")
    (out / "02_qc" / "multiqc" / "run.log").write_text("not linked\n")
    page(out / "09_report" / f"{NEW}_report_{slug}_hub_v1.html", f"Course hub: {name}",
         f"../01_gather/{NEW}_gather_{slug}_overview_v1.html")  # fmt: skip
    return root


def build(
    tmp: Path, *extra: str, real: Path | None = None
) -> tuple[Path, subprocess.CompletedProcess]:
    real = real or a_root(tmp / "real", "Study", "study")
    made = a_root(tmp / "made", "Made-up", "made")
    out = tmp / "kit"
    done = subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(out), f"--real=step={real}",
         f"--made-up=preview={made}", *extra],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    return out, done


def nothing_written(kit: Path) -> bool:
    return not kit.exists() and not kit.with_name(kit.name + ".partial").exists()


def test_kit_carries_the_newest_outputs_and_what_they_link(tmp_path):
    kit, done = build(tmp_path, "--zip")
    assert done.returncode == 0, done.stderr
    shipped = {p.relative_to(kit).as_posix() for p in kit.rglob("*") if p.is_file()}
    assert f"step/outputs/01_gather/{NEW}_gather_study_overview_v1.html" in shipped
    assert f"step/outputs/01_gather/{OLD}_gather_study_overview_v1.html" not in shipped
    assert "step/outputs/02_qc/multiqc/report.html" in shipped  # undated, but linked
    assert "step/outputs/02_qc/multiqc/run.log" not in shipped  # undated and unlinked
    record = json.loads((kit / "MANIFEST.json").read_text())
    assert {e["path"] for e in record["entries"]} == shipped - {"MANIFEST.json"}
    for e in record["entries"]:
        assert hashlib.sha256((kit / e["path"]).read_bytes()).hexdigest() == e["sha256"]
    readme = (kit / "README.md").read_text()
    assert "{" not in readme and "`step/`: Study; real data." in readme
    assert not kit.with_name("kit.partial").exists()
    with zipfile.ZipFile(kit.with_suffix(".zip")) as archive:
        modes = {i.filename: (i.external_attr >> 16) & 0o777 for i in archive.infolist()}
    assert {n.split("/")[0] for n in modes} == {"kit"} and "kit/START.html" in modes


def test_home_page_links_every_module_inside_the_kit(tmp_path, checker):
    kit, done = build(tmp_path)
    assert done.returncode == 0, done.stderr
    assert checker.check(kit) == []
    home = course.page(kit, [(kit / "step", True), (kit / "preview", False)])
    assert [r.number for r in home.rows] == [f"M{n}" for n in range(11)]
    by_module = {r.number: r.links for r in home.rows}
    # M0 is the computer's check, the same in every root: listed once, from the real one.
    assert [link.real for link in by_module["M0"]] == [True]
    assert [(link.label, link.real) for link in by_module["M1"]] == [
        ("Study", True),
        ("Made-up", False),
    ]
    assert by_module["M5"] == ()
    assert [link.label for link in home.hubs] == ["Study", "Made-up"]
    # M0 has no data: of the ten data modules only M1 is real here.
    assert [k.value for k in home.key_numbers] == ["2 of 11", "1 of 10"]
    html = (kit / course.INDEX).read_text()
    assert 'data-kind="made-up"' in html and "https://" not in html and "unzip" in html


def test_newest_page_per_dataset(tmp_path):
    folder = a_root(tmp_path / "r", "Study", "study") / "outputs" / "01_gather"
    page(folder / f"{NEW}_gather_study_overview_v2.html", "M1 Gather: Study")
    assert [p.name for p in course.newest(folder, "overview")] == [
        f"{NEW}_gather_study_overview_v2.html"
    ]


@pytest.mark.parametrize(
    ("name", "content", "reason"),
    [
        (f"{NEW}_gather_study_r1_v1.fastq.gz", b"@r\nA\n+\nI\n", "reads or alignments never ship"),
        (f"{NEW}_gather_study_big_v1.csv", None, "over the 64,000,000 cap"),
    ],
)
def test_kit_refuses_reads_and_big_files_and_writes_nothing(tmp_path, name, content, reason):
    path = a_root(tmp_path / "real", "Study", "study") / "outputs" / "01_gather" / name
    if content is None:
        with path.open(
            "wb"
        ) as handle:  # sparse: Rytis's 64 MB cap (2026-10-09) plus a byte, nothing on disk
            handle.truncate(64_000_001)
    else:
        path.write_bytes(content)
    kit, done = build(tmp_path, real=tmp_path / "real")
    assert done.returncode == 1 and reason in done.stderr and nothing_written(kit)


@pytest.mark.parametrize(
    ("extra", "reason"),
    [
        ('<a href="missing.csv">x</a>', "does not resolve"),
        ('<link rel="stylesheet" href="https://fonts.example/a.css">', "loads https://fonts"),
        ('<script src="//cdn.example/x.js"></script>', "loads //cdn.example"),
    ],
)
def test_a_bad_link_in_a_source_page_writes_nothing(tmp_path, extra, reason):
    root = a_root(tmp_path / "real", "Study", "study")
    page(
        root / "outputs" / "05_de" / f"{NEW}_de_study_overview_v1.html", "M5 DE: Study", extra=extra
    )
    kit, done = build(tmp_path, real=root)
    assert done.returncode == 1 and reason in done.stderr and nothing_written(kit)


def test_a_symlinked_outputs_folder_builds(tmp_path):
    root = a_root(tmp_path / "store", "Study", "study")
    (tmp_path / "linked").mkdir()
    (tmp_path / "linked" / "outputs").symlink_to(root / "outputs")
    kit, done = build(tmp_path, real=tmp_path / "linked")
    assert done.returncode == 0, done.stderr


def test_a_failed_check_leaves_the_partial_kit_and_no_kit(tmp_path):
    root = a_root(tmp_path / "real", "Study", "study")
    long_name = f"{NEW}_gather_study_{'x' * 180}_v1.csv"  # too long once unpacked on Windows
    (root / "outputs" / "01_gather" / long_name).write_text("a\n")
    kit, done = build(tmp_path, real=root)
    assert done.returncode == 1 and "characters after extraction" in done.stderr
    assert not kit.exists() and kit.with_name("kit.partial").is_dir()


def test_the_page_is_the_only_way_in(tmp_path):
    """No START script beside START.html: a double-clicked script runs through the reader's shell
    and the system's prompts for downloaded programs, and the page alone opens everywhere."""
    kit, done = build(tmp_path)
    assert done.returncode == 0, done.stderr
    top = sorted(p.name for p in kit.iterdir() if p.is_file())
    assert top == ["MANIFEST.json", "README.md", course.INDEX, "check_kit.ps1", "check_kit.py",
                   "kit.js"]  # fmt: skip
    assert not [p for p in kit.rglob("*") if p.is_file() and p.stat().st_mode & 0o111]


def test_start_html_warns_when_copied_out_of_the_zip(tmp_path):
    """kit.js beside START.html loads; copied out alone, its failed load shows the warning."""
    kit, done = build(tmp_path)
    assert done.returncode == 0, done.stderr
    html = (kit / "START.html").read_text()
    assert '<div class="zipped" id="ctx-zipped" role="alert" hidden>' in html
    assert (
        '<script src="kit.js" onerror="document.getElementById(\'ctx-zipped\').hidden = false">'
        in html
    )
    assert (kit / "kit.js").read_text().startswith("window.ctxKit")


def test_check_kit_finds_a_changed_file(tmp_path):
    kit, done = build(tmp_path)
    assert done.returncode == 0, done.stderr
    ok = subprocess.run([sys.executable, str(kit / "check_kit.py")], capture_output=True,
                        text=True, check=False)  # fmt: skip
    assert ok.returncode == 0, ok.stdout
    target = kit / "step" / "outputs" / "01_gather" / f"{NEW}_gather_study_counts_v1.csv"
    target.write_text("gene,n\nA,2\n")
    bad = subprocess.run([sys.executable, str(kit / "check_kit.py"), str(kit)],
                         capture_output=True, text=True, check=False)  # fmt: skip
    assert bad.returncode == 1 and "not the file the course built" in bad.stdout


def test_what_fetches_from_the_web(checker, tmp_path):
    html = tmp_path / "p.html"
    html.write_text(
        '<a href="https://www.genecards.org/">gene</a><a href="mailto:x@y.org">m</a>'
        '<img srcset="data:image/png;base64,AA 1x, https://img.example/a.png 2x">'
        '<a href="file:///Users/someone/page.html">f</a>'
        '<div style="background:url(https://bg.example/a.png)"></div>'
        "<style>@import url('https://css.example/a.css');</style>"
    )
    found = checker.link_problems(tmp_path, [html])
    assert [f.split("loads ")[1][:12] for f in found] == [
        "https://img.", "file:///User", "https://bg.e", "https://css.",
    ]  # fmt: skip


def test_windows_names(checker):
    names = ["a/ok.csv", "a/b:c.csv", "a/CON.txt", "a/name.", "a/Same.csv", "a/same.csv",
             "a/" + "x" * 220]  # fmt: skip
    found = checker.windows_problems(names, "kit")
    assert len(found) == 5 and not any(f.startswith("a/ok.csv") for f in found)


def test_the_windows_checker_reads_its_kit_and_is_plain_ascii(tmp_path):
    """check_kit.ps1 runs under Windows PowerShell 5.1, which reads a file without a BOM as ANSI,
    and must check any kit: so ASCII only, and no value of one kit (a hash, a count, a name)."""
    kit, done = build(tmp_path)
    assert done.returncode == 0, done.stderr
    text = (kit / "check_kit.ps1").read_bytes()
    assert text == (ROOT / "scripts" / "powershell" / "check_kit.ps1").read_bytes()
    assert text.isascii()
    script = text.decode()
    assert "ConvertFrom-Json" in script and "MANIFEST.json" in script
    assert not re.search(r"[0-9a-fA-F]{64}|course_kit|\b10[56][56]\b", script)


def two_datasets(root: Path) -> Path:
    """A root holding a dataset and a second one whose slug starts with the first's."""
    a_root(root, "Study", "study")
    return a_root(root, "Study subset", "study_sub")


def shipped_files(kit: Path) -> set[str]:
    return {p.relative_to(kit).as_posix() for p in kit.rglob("*") if p.is_file()}


def test_a_left_out_dataset_ships_nothing_and_a_longer_slug_ships_whole(tmp_path):
    """--leave-out drops one dataset's outputs and home links, never a longer slug's (D-050).

    Mutation caught: ownership by prefix alone, with which leaving out "study" drops "study_sub".
    """
    real = two_datasets(tmp_path / "real")
    kit, done = build(tmp_path / "a", "--leave-out", "study_sub", real=real)
    assert done.returncode == 0, done.stderr
    files = shipped_files(kit)
    assert f"step/outputs/01_gather/{NEW}_gather_study_overview_v1.html" in files
    assert not [f for f in files if "study_sub" in f]
    assert "study_sub" not in (kit / course.INDEX).read_text()
    assert json.loads((kit / "MANIFEST.json").read_text())["left_out"] == ["study_sub"]
    kit, done = build(tmp_path / "b", "--leave-out", "study", real=real)
    assert done.returncode == 0, done.stderr
    files = shipped_files(kit)
    assert f"step/outputs/01_gather/{NEW}_gather_study_sub_overview_v1.html" in files
    assert f"step/outputs/09_report/{NEW}_report_study_sub_hub_v1.html" in files
    assert f"step/outputs/01_gather/{NEW}_gather_study_overview_v1.html" not in files


def test_a_shipped_page_linking_a_left_out_dataset_writes_nothing(tmp_path):
    """Mutation caught: what the shipped pages link left unchecked, which carries it back in."""
    real = two_datasets(tmp_path / "real")
    page(real / "outputs" / "02_qc" / f"{NEW}_qc_study_overview_v1.html", "M2 Read QC: Study",
         f"../01_gather/{NEW}_gather_study_sub_overview_v1.html")  # fmt: skip
    kit, done = build(tmp_path, "--leave-out", "study_sub", real=real)
    assert done.returncode != 0 and "its dataset is left out" in done.stderr
    assert nothing_written(kit)


def test_leaving_out_a_dataset_no_root_holds_writes_nothing(tmp_path):
    """Mutation caught: an unknown slug accepted, so a typo builds a kit that leaves nothing out."""
    kit, done = build(tmp_path, "--leave-out", "stdy")
    assert done.returncode != 0 and "no root holds that dataset" in done.stderr
    assert nothing_written(kit)


def test_a_kit_folder_builds_a_kit(tmp_path):
    """U.7: a dataset folder of a kit is a root the builder takes again, with the guide's page,
    so a student's kit rebuilds from the kit alone. Catches that folder refused, and the guide
    left behind on the way."""
    real = a_root(tmp_path / "real", "Study", "study")
    page(real / "outputs" / "09_report" / f"{NEW}_report_guide_v1.html", "Your own data: Study")
    kit, done = build(tmp_path, real=real)
    assert done.returncode == 0, done.stderr
    again = tmp_path / "again"
    done = subprocess.run([sys.executable, str(SCRIPT), "--out", str(again),
                           f"--real=step={kit / 'step'}"], capture_output=True, text=True,
                          check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr
    pages = {p.relative_to(kit / "step") for p in (kit / "step").rglob("*.html")}
    assert pages == {p.relative_to(again / "step") for p in (again / "step").rglob("*.html")}
    assert course.page(again, [(again / "step", True)]).guide == (
        f"step/outputs/09_report/{NEW}_report_guide_v1.html"
    )


def release(tmp: Path, real: Path, *extra: str) -> tuple[Path, subprocess.CompletedProcess]:
    """A --release build of one real root, with any other arguments."""
    out = tmp / "kit"
    done = subprocess.run([sys.executable, str(SCRIPT), "--out", str(out), "--release",
                           f"--real=step={real}", *extra], capture_output=True, text=True,
                          check=False)  # fmt: skip
    return out, done


def test_a_release_kit_refuses_a_made_up_root_and_writes_nothing(tmp_path):
    """K.6: --release with --made-up stops before writing. Mutation caught: the release check
    dropped from build(), which ships the preview root marked Made-up."""
    real, made = a_root(tmp_path / "real", "Study", "study"), a_root(tmp_path / "made", "M", "m")
    kit, done = release(tmp_path, real, f"--made-up=preview={made}")
    assert done.returncode != 0 and "a release kit holds no made-up root" in done.stderr
    assert nothing_written(kit)


def test_a_release_kit_refuses_a_made_up_page_in_a_real_root(tmp_path):
    """K.6: a preview page put in a real root is found by the preview's own note. Mutation
    caught: made_up_pages returning nothing, so the page ships labelled real."""
    real = a_root(tmp_path / "real", "Study", "study")
    note = text.word("ui", "preview.note")
    page(real / "outputs" / "12_clinical" / f"{NEW}_clinical_study_overview_v1.html",
         "M10 Clinical and survival: Made-up study", extra=f"<p>{note}</p>")  # fmt: skip
    kit, done = release(tmp_path, real)
    assert done.returncode != 0 and "a made-up page; a release kit holds none" in done.stderr
    assert nothing_written(kit)


def test_a_release_kit_counts_zero_made_up_links_and_says_so(tmp_path):
    """K.6: the release's MANIFEST.json and README.md state it, with 0 Made-up links counted;
    a development kit counts its own. Mutation caught: made_up_links returning 0 whatever the
    page holds, which a development kit with a preview root shows."""
    real = a_root(tmp_path / "real", "Study", "study")
    kit, done = release(tmp_path, real)
    assert done.returncode == 0, done.stderr
    record = json.loads((kit / "MANIFEST.json").read_text())
    assert (record["release"], record["made_up_links"]) == (True, 0)
    assert text.word("course", "readme.release") in (kit / "README.md").read_text()
    assert "release kit: 0 Made-up links" in done.stdout
    dev, done = build(tmp_path / "dev")
    assert done.returncode == 0, done.stderr
    record = json.loads((dev / "MANIFEST.json").read_text())
    assert record["release"] is False and record["made_up_links"] > 0
    assert text.word("course", "readme.release") not in (dev / "README.md").read_text()


def ctx_readme(table: Path, rows: dict[str, str]) -> None:
    """A table README as ctx.io.readme_companion writes it: the rows, then the file table."""
    lines = ["| Field | Value |", "| --- | --- |", *(f"| {k} | {v} |" for k, v in rows.items())]
    files = [f"| `{p.name}` | {p.stat().st_size} | `{'a' * 64}` |"
             for p in sorted(table.parent.glob(table.name.split(".")[0] + ".*"))]  # fmt: skip
    body = "\n".join([f"# {table.stem}", "", *lines, "", "## Files", "", *files, ""])
    table.with_name(table.name.split(".")[0] + "_README.md").write_text(body)


def gdc_root(tmp: Path, dataset_file: tuple[str, str] | None = None) -> tuple[Path, Path]:
    """Step-1's root (its counts with a README naming no GDC folder) plus a GDC dataset: its
    gather counts, CSV and Parquet, read from a GDC folder that holds a manifest and a fetch
    script, linked from its M1 page; and a DE table whose README names the folder too.
    ``dataset_file``, (path, SHA-256), adds a ``--dataset-file`` and its recorded hash row."""
    root = a_root(tmp / "real", "Study", "study")
    out = root / "outputs"
    ctx_readme(out / "01_gather" / f"{NEW}_gather_study_counts_v1.csv",
               {"Command": "`ctx gather --dataset study`"})  # fmt: skip
    folder = tmp / "raw" / "tcga_x"
    folder.mkdir(parents=True)
    (folder / "gdc_manifest_star_counts.txt").write_text("id\tfilename\n")
    (folder / "_fetch_tcga_x.sh").write_text("echo fetch\n")
    gather = out / "01_gather"
    counts = gather / f"{NEW}_gather_tcgax_counts_v1.csv"
    counts.write_text("gene_id,A\nG1,1\n")
    counts.with_suffix(".parquet").write_bytes(b"PAR1")
    command = f"`/x/src/ctx/__main__.py --root {root} gather --dataset tcgax --gdc-dir {folder}`"
    rows = {"Command": command, "GDC folder": str(folder)}
    if dataset_file is not None:
        path, sha = dataset_file
        rows["Command"] = command[:-1] + f" --dataset-file {path}`"
        rows[f"sha256:{Path(path).name}"] = sha
    ctx_readme(counts, rows)
    page(gather / f"{NEW}_gather_tcgax_overview_v1.html", "M1 Gather: TCGA x", counts.name)
    de = out / "05_de" / f"{NEW}_de_tcgax_results_v1.csv"
    page(out / "05_de" / f"{NEW}_de_tcgax_overview_v1.html", "M5 DE: TCGA x", de.name)
    de.write_text("gene_id,padj\nG1,0.1\n")
    ctx_readme(de, {"Command": "`ctx de --dataset tcgax`", "GDC folder": str(folder)})
    return root, folder


def test_a_gdc_count_matrix_stays_out_with_its_rebuild_commands(tmp_path, checker):
    """Rytis 2026-10-09: a raw GDC count matrix never ships; its manifest and fetch script do,
    and LEFT_OUT.json, README.md and MANIFEST.json name it with its bytes and the commands that
    rebuild it; a page linking it still passes the kit's own check. Mutations caught: the matrix
    shipped; a derived table left out; Step-1's counts left out; the gather command missing
    from the README."""
    real, folder = gdc_root(tmp_path)
    kit, done = build(tmp_path / "k", real=real)
    assert done.returncode == 0, done.stderr
    files = shipped_files(kit)
    gathered = "step/outputs/01_gather"
    assert f"{gathered}/{NEW}_gather_tcgax_counts_v1.csv" not in files
    assert f"{gathered}/{NEW}_gather_tcgax_counts_v1.parquet" not in files
    assert f"{gathered}/{NEW}_gather_tcgax_counts_v1_README.md" in files
    assert f"step/outputs/05_de/{NEW}_de_tcgax_results_v1.csv" in files
    assert f"{gathered}/{NEW}_gather_study_counts_v1.csv" in files
    assert {
        "step/gdc/tcga_x/gdc_manifest_star_counts.txt",
        "step/gdc/tcga_x/_fetch_tcga_x.sh",
    } <= files
    left = json.loads((kit / "LEFT_OUT.json").read_text())["files"]
    assert sorted(e["path"] for e in left) == [
        f"{gathered}/{NEW}_gather_tcgax_counts_v1.csv",
        f"{gathered}/{NEW}_gather_tcgax_counts_v1.parquet",
    ]
    gather_command = f"ctx --root {real} gather --dataset tcgax --gdc-dir {folder}"
    for e in left:
        assert (
            e["bytes"]
            == (real / "outputs" / e["path"].removeprefix("step/outputs/")).stat().st_size
        )
        assert e["rebuild"] == [
            f"OUT={folder} bash step/gdc/tcga_x/_fetch_tcga_x.sh",
            gather_command,
        ]
    readme = (kit / "README.md").read_text()
    assert gather_command in readme and all(e["path"] in readme for e in left)
    assert json.loads((kit / "MANIFEST.json").read_text())["left_out_files"] == left
    assert checker.check(kit) == []


def test_a_source_file_never_ships_whatever_its_name(tmp_path, kit_script, monkeypatch):
    """M7 plan, decision 3: a file equal to a declared resource stops the build by its bytes.

    Mutations caught: the check removed from select; build passing no sources to select; the
    sources read from no dataset file (the CPTAC freeze and Rieckmann's table missing)."""
    root = a_root(tmp_path / "real", "Study", "study")
    copy = root / "outputs" / "07_proteo" / f"{NEW}_proteo_study_table_v1.csv"
    copy.parent.mkdir(parents=True)
    copy.write_text("a freeze, renamed\n")
    page(root / "outputs" / "07_proteo" / f"{NEW}_proteo_study_overview_v1.html", "M7: Study",
         copy.name)  # fmt: skip
    planted = {hashlib.sha256(copy.read_bytes()).hexdigest(): "fake.toml freeze"}
    _, _, problems = kit_script.select(root, sources=planted)
    assert problems == [f"{copy}: the source file fake.toml freeze; derived tables ship, "
                        "sources never"]  # fmt: skip
    monkeypatch.setattr(kit_script, "declared_sources", lambda: planted)
    with pytest.raises(SystemExit, match="the source file fake.toml freeze"):
        kit_script.build(tmp_path / "kit", [("step", root, True)])
    assert nothing_written(tmp_path / "kit")


def test_the_declared_sources_are_every_dataset_files_resources(kit_script):
    """Mutation caught: the sources read from no dataset file, or from one alone."""
    declared = kit_script.declared_sources()
    cptac = "c27a0226cba4f14c96b905f86c30007b62dd2a67a176e608c96f642566862903"
    assert declared[cptac] == "cptac_brca.toml cptac_freeze"
    assert any(name.startswith("tier1.toml ") for name in declared.values())


def committed(name: str, at: str = "HEAD") -> str:
    """The SHA-256 of a dataset file as one commit of this repository holds it."""
    shown = subprocess.run(["git", "-C", str(ROOT), "show", f"{at}:conf/datasets/{name}"],
                           capture_output=True, check=True).stdout  # fmt: skip
    return hashlib.sha256(shown).hexdigest()


@pytest.mark.parametrize("name", ["tier1.toml", "tcga_kirc.toml"])
def test_every_committed_version_of_a_dataset_file_is_the_codebases(kit_script, name):
    """K.9: a gather run reads whichever version its commit held, so every version counts,
    a merge's included: ~/ctx_kidney's TCGA gather read the version merge 789f8ad made.
    Mutations caught: only the newest commit read (git log -n 1); merges skipped (no -m)."""
    codebase = kit_script.codebase_datasets()
    commits = subprocess.run(["git", "-C", str(ROOT), "log", "--format=%H", "--",
                              f"conf/datasets/{name}"], capture_output=True, text=True,
                             check=True).stdout.split()  # fmt: skip
    assert len(commits) > 1
    for at in commits:
        assert codebase[committed(name, at)].path == f"conf/datasets/{name}", at


@pytest.mark.parametrize("ours", [True, False])
def test_a_rebuild_command_names_the_codebases_dataset_file(tmp_path, ours):
    """K.9: kit v8 showed the build clone's folder as the reader's dataset file. A
    --dataset-file whose recorded bytes a commit holds, from any clone, is shown as
    conf/datasets/<file>; one no commit holds keeps its path, and README.md names it the
    builder's own; README.md names the builder's ctx root too; the table's own README, the
    provenance, ships byte for byte. Mutations caught: no mapping; a mapping by path, which
    takes a reader's file saved in conf/datasets/ of a clone for the codebase's; the root
    sentence dropped."""
    clone = tmp_path / "some_clone" / "conf" / "datasets" / "tier1.toml"
    sha = committed("tier1.toml") if ours else "f" * 64
    real, folder = gdc_root(tmp_path, (str(clone), sha))
    kit, done = build(tmp_path / "k", real=real)
    assert done.returncode == 0, done.stderr
    left = json.loads((kit / "LEFT_OUT.json").read_text())["files"]
    shown = "conf/datasets/tier1.toml" if ours else str(clone)
    gather = f"ctx --root {real} gather --dataset tcgax --gdc-dir {folder} --dataset-file {shown}"
    assert [e["rebuild"][-1] for e in left] == [gather, gather]
    assert {e["builder_dataset_file"] for e in left} == {None if ours else str(clone)}
    readme = (kit / "README.md").read_text()
    root_line = text.fill("course", "readme.left_out_root", root=str(real)).strip()
    own_line = text.fill("course", "readme.left_out_own", path=str(clone)).strip()
    assert readme.count(root_line) == 2 and (own_line in readme) is not ours
    name = f"{NEW}_gather_tcgax_counts_v1_README.md"
    shipped = kit / "step" / "outputs" / "01_gather" / name
    assert shipped.read_bytes() == (real / "outputs" / "01_gather" / name).read_bytes()


def test_a_root_note_is_written_by_the_build(tmp_path):
    """K.10: kit v8's M9 note was appended to README.md by hand and MANIFEST.json rewritten.
    Mutations caught: the note left out of README.md; left out of MANIFEST.json."""
    note = "M9 on this study shows the clinical genes only."
    kit, done = build(tmp_path, f"--note=step={note}")
    assert done.returncode == 0, done.stderr
    readme = (kit / "README.md").read_text()
    assert f"`step/`: Study; real data. {note}\n" in readme
    assert "`preview/`: Made-up; made-up data.\n" in readme
    roots = json.loads((kit / "MANIFEST.json").read_text())["roots"]
    assert [(r["folder"], r["note"]) for r in roots] == [("step", note), ("preview", None)]


@pytest.mark.parametrize(
    ("note", "reason"),
    [
        ("nowhere=M9 shows the clinical genes only.", "no --real or --made-up folder"),
        ("step=M9 is very short.", "banned word 'very'"),
        ("step=Is M9 short?", "--note step"),
        ("step=", "one note per folder, and not empty"),
    ],
)
def test_a_bad_note_writes_nothing(tmp_path, note, reason):
    """K.10: a note naming no root, failing the prose lint, or empty stops the build first.
    Mutations caught: the folder unchecked; the prose lint skipped."""
    kit, done = build(tmp_path, f"--note={note}")
    assert done.returncode != 0 and reason in done.stderr and nothing_written(kit)


def test_a_two_project_count_matrix_stays_out(tmp_path, checker):
    """A sheet joining two GDC projects names each folder after its project (D-076); the
    matrix it read stays out all the same. Kit v9's tcga_lung counts, 192,637,178 B, were
    taken for a file to ship and stopped the build at the cap (2026-10-16)."""
    real, folder = gdc_root(tmp_path)
    readme = real / "outputs" / "01_gather" / f"{NEW}_gather_tcgax_counts_v1_README.md"
    readme.write_text(readme.read_text().replace("| GDC folder |", "| PROJ-X GDC folder |"))
    kit, done = build(tmp_path / "k", real=real)
    assert done.returncode == 0, done.stderr
    files = shipped_files(kit)
    assert f"step/outputs/01_gather/{NEW}_gather_tcgax_counts_v1.csv" not in files
    assert "step/gdc/tcga_x/gdc_manifest_star_counts.txt" in files


def test_a_dataset_ships_its_newest_run_not_each_names_newest_copy(tmp_path):
    """A run writes its files on one day, so a dataset's newest day in a folder is its run.
    Kit v9 kept the newest copy per file name and shipped 114 files of an older M9 run whose
    genes the newer run did not draw (2026-10-16). Mutation caught: the per-name rule back."""
    real = a_root(tmp_path / "real", "Study", "study")
    cards = real / "outputs" / "11_genecard"
    for gene in ("a", "b"):
        page(cards / f"{OLD}_genecard_study_card_{gene}_v1.html", f"M9 card {gene}")
    old_cards = [f"{OLD}_genecard_study_card_{gene}_v1.html" for gene in ("a", "b")]
    page(cards / f"{OLD}_genecard_study_overview_v1.html", "M9 Gene cards: Study", *old_cards)
    page(cards / f"{NEW}_genecard_study_card_c_v1.html", "M9 card c")
    page(cards / f"{NEW}_genecard_study_overview_v1.html", "M9 Gene cards: Study",
         f"{NEW}_genecard_study_card_c_v1.html")  # fmt: skip
    page(cards / "index.html", "M9", f"{NEW}_genecard_study_overview_v1.html")
    kit, done = build(tmp_path / "k", real=real)
    assert done.returncode == 0, done.stderr
    shipped = {Path(f).name for f in shipped_files(kit) if "/11_genecard/" in f}
    assert f"{NEW}_genecard_study_card_c_v1.html" in shipped
    assert not {n for n in shipped if n.startswith(OLD)}, shipped


def test_m0_shows_the_computers_check_once_across_real_roots(tmp_path):
    """Every root's M0 page describes the one computer that built the kit, so the home page
    lists one, the newest (Rytis 2026-10-20); on a date tie, the first real root in the kit's
    order (K.17: kit v10 linked pleura's, the last of eight). Mutations caught: one link per
    real root, the last root on a tie, and the first root whatever its date."""
    kit = tmp_path / "kit"
    for folder, day in (("a", OLD), ("b", NEW), ("c", NEW), ("d", NEW)):
        a_root(kit / folder, folder.upper(), folder)
        setup = kit / folder / "outputs" / "00_setup"
        for old in setup.glob("*.html"):
            old.unlink()
        page(setup / f"{day}_setup_environment_v1.html", "M0 Setup")
    roots = [(kit / "a", True), (kit / "d", False), (kit / "b", True), (kit / "c", True)]
    (m0,) = [r.links for r in course.page(kit, roots).rows if r.number == "M0"]
    assert [link.href for link in m0] == [f"b/outputs/00_setup/{NEW}_setup_environment_v1.html"]


def stale_guide_root(root: Path) -> Path:
    """K.15's ~/ctx_tcga_hnsc: module pages from NEW, a guide from OLD linking OLD's glossary."""
    a_root(root, "Study", "study")
    report = root / "outputs" / "09_report"
    page(report / f"{OLD}_report_glossary_v1.html", "Glossary")
    page(
        report / f"{OLD}_report_guide_v1.html",
        "Your own data: Study",
        f"{OLD}_report_glossary_v1.html",
        extra='<code id="ctx-command">pixi run ctx --root ~/r guide --dataset-file x.toml</code>',
    )
    page(report / f"{NEW}_report_glossary_v1.html", "Glossary")
    hub = report / f"{NEW}_report_study_hub_v1.html"
    page(hub, "Course hub: Study", f"../01_gather/{NEW}_gather_study_overview_v1.html",
         f"{NEW}_report_glossary_v1.html")  # fmt: skip
    return root


def test_kit_refuses_a_guide_older_than_its_roots_module_pages(tmp_path):
    """K.15: kit v10 shipped ~/ctx_tcga_hnsc's 2026-10-18 guide, and the glossary it links,
    beside that root's 2026-10-20 pages, because `ctx guide` is no module and no run wrote it
    again. The refusal names both files and the guide's own command. Mutation caught: no check."""
    kit, done = build(tmp_path, real=stale_guide_root(tmp_path / "real"))
    assert done.returncode != 0 and nothing_written(kit)
    stale = [line for line in done.stderr.splitlines() if "older than" in line]
    assert [Path(line.split(":")[0]).name for line in stale] == [
        f"{OLD}_report_glossary_v1.html", f"{OLD}_report_guide_v1.html"]  # fmt: skip
    assert "`pixi run ctx --root ~/r guide --dataset-file x.toml`" in stale[1]
    assert f"{OLD}_report_guide_v1.html links it" in stale[0]


def test_kit_ships_an_older_glossary_a_module_page_links(tmp_path):
    """A module not run again keeps its day's glossary, its own copy (D-048): breast M9 rebuilt
    a day after breast's other pages must still ship (K.13b). Mutation caught: every guide or
    glossary older than the root's newest module page refused, linked or not."""
    real = stale_guide_root(tmp_path / "real")
    report = real / "outputs" / "09_report"
    (report / f"{OLD}_report_guide_v1.html").unlink()
    page(report / f"{NEW}_report_guide_v1.html", "Your own data: Study",
         f"{NEW}_report_glossary_v1.html")  # fmt: skip
    old_qc = real / "outputs" / "02_qc" / f"{OLD}_qc_study_overview_v1.html"
    page(old_qc, "M2 QC: Study", f"../09_report/{OLD}_report_glossary_v1.html")
    kit, done = build(tmp_path, real=real)
    assert done.returncode == 0, done.stderr
    assert (kit / "step" / "outputs" / "09_report" / f"{OLD}_report_glossary_v1.html").is_file()
