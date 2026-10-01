"""scripts/python/build_kit.py, check_kit.py and ctx.report.course: what a student kit holds."""

import hashlib
import importlib.util
import json
import os
import re
import shutil
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
    assert modes["kit/START.command"] == modes["kit/START.sh"] == 0o755


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
        (f"{NEW}_gather_study_big_v1.csv", None, "over the 50,000,000 cap"),
    ],
)
def test_kit_refuses_reads_and_big_files_and_writes_nothing(tmp_path, name, content, reason):
    path = a_root(tmp_path / "real", "Study", "study") / "outputs" / "01_gather" / name
    if content is None:
        with path.open("wb") as handle:  # sparse: 50 MB on paper, nothing on disk
            handle.truncate(50_000_001)
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


def test_start_scripts_open_the_page_beside_them(tmp_path):
    kit, done = build(tmp_path)
    assert done.returncode == 0, done.stderr
    top = sorted(p.name for p in kit.iterdir() if p.is_file())
    assert top == ["MANIFEST.json", "README.md", "START.cmd", "START.command", course.INDEX,
                   "START.sh", "check_kit.ps1", "check_kit.py", "kit.js"]  # fmt: skip
    runnable = {p.name for p in kit.rglob("*") if p.is_file() and p.stat().st_mode & 0o111}
    assert runnable == {"START.command", "START.sh"}
    # Run each script from elsewhere, as a double-click does, with open and xdg-open stubbed.
    stubs = tmp_path / "bin"
    stubs.mkdir()
    for opener in ("open", "xdg-open"):
        (stubs / opener).write_text('#!/bin/sh\nprintf "%s/%s" "$PWD" "$1" > "$OPENED"\n')
        (stubs / opener).chmod(0o755)
    zipped = text.word("course", "start.zipped")
    for script in ("START.command", "START.sh"):
        opened = tmp_path / f"{script}.opened"
        env = {**os.environ, "PATH": f"{stubs}:{os.environ['PATH']}", "OPENED": str(opened)}
        subprocess.run([str(kit / script)], cwd=tmp_path, env=env, check=True)
        assert Path(opened.read_text()).resolve() == (kit / "START.html").resolve()
        # Copied out of the zip alone, as Windows does with a file opened inside one: no page
        # beside it, so it says why and opens nothing.
        alone = tmp_path / f"alone_{script}"
        alone.mkdir()
        shutil.copy2(kit / script, alone / script)
        opened.unlink()
        run = subprocess.run([str(alone / script)], cwd=tmp_path, env=env, capture_output=True,
                             text=True, check=False)  # fmt: skip
        assert run.returncode == 1 and zipped in run.stdout and not opened.exists()
    cmd = (kit / "START.cmd").read_bytes()
    assert cmd.count(b"\r\n") == cmd.count(b"\n")  # cmd.exe needs CRLF to find a goto label
    lines = cmd.decode().split("\r\n")
    assert lines[1] == 'if exist "%~dp0START.html" goto open'
    assert lines[2] == f"echo {zipped}" and "pause" in lines and ":open" in lines
    assert lines[lines.index(":open") + 1] == 'start "" "%~dp0START.html"'


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


def test_a_start_message_cmd_or_sh_would_act_on_is_refused(tmp_path):
    content = tmp_path / "content"
    shutil.copytree(ROOT / "src" / "ctx" / "report" / "content", content)
    course_toml = content / "course.toml"
    zipped = text.word("course", "start.zipped")
    course_toml.write_text(course_toml.read_text().replace(zipped, zipped + " & del *"))
    real, made = a_root(tmp_path / "real", "Study", "study"), a_root(tmp_path / "made", "M", "made")
    done = subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(tmp_path / "kit"), f"--real=step={real}",
         f"--made-up=preview={made}"],
        capture_output=True, text=True, check=False,
        env={**os.environ, "CTX_CONTENT": str(content)},
    )  # fmt: skip
    assert done.returncode == 1 and "a character cmd or sh would act on" in done.stderr


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
