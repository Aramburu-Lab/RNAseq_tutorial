"""The ctx CLI end to end: `ctx setup`, parameter cards, profiles, and code kept apart from data."""

import hashlib
import os
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import polars as pl
import pytest
import typer
from test_packaging import build_wheel, run_from_wheel

from ctx.cli import app, parameter_cards
from ctx.io.naming import MODULES, lint
from ctx.runtime import BUDGETS

ROOT = Path(__file__).resolve().parents[1]
CTX = Path(sys.executable).parent / "ctx"
NO_SLURM = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}


def _ctx(
    *args: str, cwd: Path, env: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(CTX), *args], cwd=cwd, env={**NO_SLURM, **(env or {})},
        capture_output=True, text=True, check=False,
    )  # fmt: skip


def _report(root: Path) -> str:
    (page,) = (root / "outputs" / "00_setup").glob("*_setup_environment_v1.html")
    return page.read_text()


@pytest.fixture(scope="module")
def setup_run(tmp_path_factory):
    root = tmp_path_factory.mktemp("root")
    done = _ctx("--root", str(root), "setup", cwd=root)
    assert done.returncode == 0, done.stderr
    return root, _report(root)


FIXTURE = ROOT / "tests" / "fixtures" / "fetchngs" / "download"


@pytest.fixture(scope="module")
def gather_run(tmp_path_factory):
    root = tmp_path_factory.mktemp("gather_root")
    args = (
        "--root",
        str(root),
        "gather",
        "--dataset",
        "nfcore_test",
        "--fetchngs-dir",
        str(FIXTURE),
    )
    done = _ctx(*args, cwd=root)
    assert done.returncode == 0, done.stderr
    (page,) = (root / "outputs" / "01_gather").glob("*_gather_nfcore_test_overview_v1.html")
    return root, page.read_text()


def test_every_subcommand_has_a_card_test():
    # A new subcommand without a case below fails here, not silently.
    assert set(typer.main.get_command(app).commands) == {"setup", "gather"}


@pytest.mark.parametrize("command", ["setup", "gather"])
def test_parameter_cards_match_cli(request, command):
    _, html = request.getfixturevalue(f"{command}_run")
    group = typer.main.get_command(app)
    options = {
        p.opts[0]: p
        for cmd in (group, group.commands[command])
        for p in cmd.params
        if p.opts[0] not in {"--help", "--version"}
    }
    assert options
    shown = dict(
        re.findall(r'<tr id="card-(--[\w-]+)">.*?<code class="default">(.*?)</code>', html)
    )
    assert set(shown) == set(options)
    for name, param in options.items():
        expected = (
            "required"
            if param.required
            else ("unset" if param.default is None else str(param.default))
        )
        assert shown[name] == expected, name
    if command != "setup":
        return
    extra = typer.Typer()

    @extra.command()
    def uncarded(flag: int = 3) -> None:
        """A command whose option has no card."""

    with pytest.raises(ValueError, match="--flag"):
        parameter_cards(typer.main.get_command(extra))


def test_subcommands_in_modules():
    tokens = {token for token, _ in MODULES}
    names = list(typer.main.get_command(app).commands)
    assert names
    assert all(name.replace("-", "_") in tokens for name in names)


def test_threads_set_before_import(tmp_path):
    threads = 2 if os.cpu_count() != 2 else 3
    slurm = {"CTX_PROFILE": "hpc", "SLURM_CPUS_PER_TASK": str(threads), "SLURM_MEM_PER_CPU": "500"}
    done = _ctx("--root", str(tmp_path), "setup", cwd=tmp_path, env=slurm)
    assert done.returncode == 0, done.stderr
    assert f"polars thread pool: {threads}; BLAS thread pools: [{threads}]" in _report(tmp_path)


def test_help_without_slurm(tmp_path):
    hpc = {"CTX_PROFILE": "hpc"}
    assert _ctx("setup", "--help", cwd=tmp_path, env=hpc).returncode == 0
    no_root = _ctx("setup", cwd=tmp_path, env=hpc)
    assert no_root.returncode != 0
    assert "--root or set CTX_ROOT" in no_root.stderr
    no_job = _ctx("--root", str(tmp_path), "setup", cwd=tmp_path, env=hpc)
    assert no_job.returncode != 0
    assert "SLURM_CPUS_PER_TASK" in no_job.stderr


def _snapshot(directory: Path) -> dict[str, float]:
    return {str(p): p.stat().st_mtime for p in directory.rglob("*") if p.is_file()}


def test_code_is_not_data(tmp_path):
    root, elsewhere = tmp_path / "root", tmp_path / "elsewhere"
    elsewhere.mkdir()
    before = _snapshot(ROOT / "src")
    done = _ctx("setup", cwd=elsewhere, env={"CTX_ROOT": str(root)})
    assert done.returncode == 0, done.stderr
    written = [Path(line) for line in done.stdout.split()]
    assert written
    assert all(path.is_relative_to(root.resolve()) for path in written)
    assert not list(elsewhere.iterdir())
    assert _snapshot(ROOT / "src") == before
    (readme,) = (root / "outputs" / "00_setup").glob("*_README.md")
    assert f"| Directory | {(root / 'outputs' / '00_setup').resolve()} |" in readme.read_text()


def test_repro_block_shows_measured(tmp_path):
    stubbed = textwrap.dedent(f"""
        import ctx.runtime as runtime
        runtime.measure = lambda started: runtime.Usage(12.345, 678 * 2**20, 0)
        from ctx.cli import app
        app(["--root", {str(tmp_path)!r}, "setup"], standalone_mode=False)
    """)
    subprocess.run([sys.executable, "-c", stubbed], env=NO_SLURM, check=True, capture_output=True)
    html = _report(tmp_path)
    assert '<td class="wall">12.345 s</td>' in html
    assert f'<td class="rss-self">678.0 MiB ({678 * 2**20} bytes)</td>' in html


def test_setup_end_to_end(setup_run):
    root, html = setup_run
    out = root / "outputs" / "00_setup"
    scanned, violations = lint(out)
    assert scanned >= 3
    assert violations == []
    for data in [*out.glob("*.csv"), *out.glob("*.parquet")]:
        assert (out / f"{data.stem}_README.md").is_file()
    assert not re.search(r"<script[^>]+\bsrc\s*=|<link[^>]+\bhref\s*=", html, re.IGNORECASE)
    wall = float(re.search(r'<td class="wall">([\d.]+) s</td>', html).group(1))
    rss = int(re.search(r'<td class="rss-self">[\d.]+ MiB \((\d+) bytes\)</td>', html).group(1))
    assert 0 < wall < BUDGETS["B"].minutes * 60
    assert 0 < rss < BUDGETS["B"].ram_gb * 10**9


def test_setup_runs_from_wheel(tmp_path):
    site = build_wheel(tmp_path)
    run_from_wheel(
        site,
        f"""
        from ctx.cli import app
        app(["--root", {str(tmp_path / "root")!r}, "setup"], standalone_mode=False)
        """,
        cwd=tmp_path,
    )
    assert "installed wheel" in next((tmp_path / "root").rglob("*_README.md")).read_text()


def _snapshot_hashes(directory):
    return {
        str(p): (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
        for p in directory.rglob("*")
        if p.is_file()
    }


def test_gather_end_to_end(gather_run):
    root, html = gather_run
    out = root / "outputs" / "01_gather"
    scanned, violations = lint(out)
    assert scanned == 4
    assert violations == []
    (csv_path,) = out.glob("*_samplesheet_v1.csv")
    (parquet,) = out.glob("*_samplesheet_v1.parquet")
    assert (out / f"{csv_path.stem}_README.md").is_file()
    assert b"\r" not in csv_path.read_bytes()
    parquet_frame = pl.read_parquet(parquet)
    csv_frame = pl.read_csv(csv_path, schema=parquet_frame.schema)
    assert csv_frame.equals(parquet_frame)
    assert parquet_frame.height == 4
    assert sorted(parquet_frame["status_1"]) == ["no_md5", "no_md5", "verified", "verified"]
    wall = float(re.search(r'<td class="wall">([\d.]+) s</td>', html).group(1))
    rss = int(re.search(r'<td class="rss-self">[\d.]+ MiB \((\d+) bytes\)</td>', html).group(1))
    assert 0 < wall < BUDGETS["B"].minutes * 60
    assert 0 < rss < BUDGETS["B"].ram_gb * 10**9


def test_upstream_untouched(tmp_path):
    before = _snapshot_hashes(FIXTURE)
    args = ("--root", str(tmp_path), "gather", "--dataset", "x", "--fetchngs-dir", str(FIXTURE))
    assert _ctx(*args, cwd=tmp_path).returncode == 0
    assert _snapshot_hashes(FIXTURE) == before


def test_failure_writes_no_sheet(tmp_path):
    broken = tmp_path / "broken"
    shutil.copytree(FIXTURE, broken)
    (broken / "fastq" / "ERX1234253_ERR1160846.fastq.gz").unlink()
    root = tmp_path / "root"
    done = _ctx(
        "--root",
        str(root),
        "gather",
        "--dataset",
        "bad",
        "--fetchngs-dir",
        str(broken),
        cwd=tmp_path,
    )
    assert done.returncode != 0
    assert "ERX1234253_ERR1160846.fastq.gz is missing" in done.stderr
    assert not (root / "outputs").exists()


def test_gather_runs_from_wheel(tmp_path):
    site = build_wheel(tmp_path)
    run_from_wheel(
        site,
        f"""
        from ctx.cli import app
        app(["--root", {str(tmp_path / "root")!r}, "gather", "--dataset", "w",
             "--fetchngs-dir", {str(FIXTURE)!r}], standalone_mode=False)
        """,
        cwd=tmp_path,
    )
    assert list((tmp_path / "root" / "outputs" / "01_gather").glob("*_samplesheet_v1.csv"))
