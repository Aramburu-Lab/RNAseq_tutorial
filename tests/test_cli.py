"""The ctx CLI end to end: `ctx setup`, parameter cards, profiles, and code kept apart from data."""

import hashlib
import os
import re
import shutil
import subprocess
import sys
import textwrap
from datetime import date
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


DE_SPEC = ROOT / "tests" / "fixtures" / "datasets" / "de_synthetic.toml"


@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    """de, function and report on one root, from gather outputs written by gather's writer."""
    from synth import (  # noqa: PLC0415
        DE_GROUPS,
        design_counts,
        write_gather_outputs,
        write_gene_sets,
    )

    from ctx.dataset import load  # noqa: PLC0415

    root = tmp_path_factory.mktemp("pipeline")
    write_gather_outputs(root, "desynth", design_counts(DE_GROUPS, batch_lfc=2.0), date.today())
    write_gene_sets(root / load(DE_SPEC).resource("sets").dest)
    gathered = {p.name: p.read_bytes() for p in (root / "outputs" / "01_gather").iterdir()}
    runs = {}
    with_spec = ["--dataset-file", str(DE_SPEC)]
    for command, extra in (("de", with_spec), ("function", with_spec), ("report", with_spec)):
        done = _ctx("--root", str(root), command, "--dataset", "desynth", *extra, cwd=root)
        assert done.returncode == 0, (command, done.stderr[-3000:])
        runs[command] = done
    return root, gathered, runs


def _page(root: Path, number_module: str, what: str) -> str:
    (page,) = (root / "outputs" / number_module).glob(f"*_{what}_v1.html")
    return page.read_text()


@pytest.fixture(scope="module")
def de_run(pipeline):
    return pipeline[0], _page(pipeline[0], "05_de", "desynth_overview")


@pytest.fixture(scope="module")
def function_run(pipeline):
    return pipeline[0], _page(pipeline[0], "06_function", "desynth_overview")


@pytest.fixture(scope="module")
def report_run(pipeline):
    return pipeline[0], _page(pipeline[0], "09_report", "desynth_hub")


COMMANDS = ["setup", "gather", "de", "function", "report"]


def test_every_subcommand_has_a_card_test():
    # A new subcommand without a case below fails here, not silently.
    assert set(typer.main.get_command(app).commands) == set(COMMANDS)


@pytest.mark.parametrize("command", COMMANDS)
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
        if param.required:
            expected = "required"
        elif param.secondary_opts:  # an on/off switch reads "on" or "off", never "True"
            expected = "on" if param.default else "off"
        else:
            expected = "unset" if param.default is None else str(param.default)
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
    # The key number counts files really in fastq/, not rows that list one.
    on_disk = len(list((FIXTURE / "fastq").glob("*.fastq.gz")))
    shown = re.search(r"<dt>FASTQ files on disk</dt><dd class=\"key-value\">([\d,]+)</dd>", html)
    assert shown and int(shown.group(1).replace(",", "")) == on_disk == 6
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


def _offline(html: str) -> None:
    """P1.8: plotly.js inlined at most once, no external script, link or url()."""
    import plotly.offline  # noqa: PLC0415

    js = plotly.offline.get_plotlyjs()
    assert html.count(js) <= 1
    outside = html.replace(js, "")
    assert not re.search(r"<script[^>]+\bsrc\s*=", outside, re.IGNORECASE)
    assert not re.search(r"<link[^>]+\bhref\s*=", outside, re.IGNORECASE)
    assert not re.search(r"url\(\s*['\"]?(https?:)?//", outside, re.IGNORECASE)
    assert len(html.encode()) <= 15 * 10**6


@pytest.mark.parametrize("folder", ["05_de", "06_function", "09_report"])
def test_module_outputs_clean(pipeline, folder):
    root, _, _ = pipeline
    out = root / "outputs" / folder
    scanned, violations = lint(out)
    assert scanned > 0 and violations == []
    for data in [*out.glob("*.csv"), *out.glob("*.parquet")]:
        stem = data.name.rsplit(".", 1)[0]
        assert (out / f"{stem}_README.md").is_file(), data.name
    for page in out.glob("*.html"):
        _offline(page.read_text())


def test_de_end_to_end(de_run):
    root, html = de_run
    out = root / "outputs" / "05_de"
    results = pl.read_parquet(next(out.glob("*_desynth_a_case_results_v1.parquet")))
    planted = [f"g{i + 1:05d}" for i in range(100)]
    called = set(results.filter(pl.col("call") != "ns")["gene_id"])
    assert len(called & set(planted)) >= 80
    metrics = dict(
        pl.read_parquet(next(out.glob("*_desynth_metrics_v1.parquet")))
        .select("metric", "value")
        .rows()
    )
    # The registered S1 rule, not zero: a 3-vs-3 null split called 1 to 4 genes on half the
    # seeds swept, as a Wald test on estimated dispersions does (D-017).
    assert metrics["calls.a_null.significant"] <= 0.01 * metrics["calls.a_null.tested"]
    # Recounted from the table. That significance counts any fold change, unlike up and down,
    # is tested on a constructed gene in test_de.py::test_significant_counts_any_fold_change.
    significant = int((results["padj"] < 0.05).sum())
    up, down = int((results["call"] == "up").sum()), int((results["call"] == "down").sum())
    assert metrics["calls.a_case.significant"] == significant
    assert (metrics["calls.a_case.up"], metrics["calls.a_case.down"]) == (up, down)
    assert metrics["calls.a_case.tested"] == results["padj"].is_not_null().sum()
    assert metrics["sex.disagree"] == 0 and metrics["sex.compared"] == 33
    wall = float(re.search(r'<td class="wall">([\d.]+) s</td>', html).group(1))
    assert 0 < wall < BUDGETS["B"].minutes * 60


def test_function_end_to_end(function_run):
    root, _ = function_run
    out = root / "outputs" / "06_function"
    gsea = dict(
        (s, (nes, fdr))
        for s, nes, fdr in pl.read_parquet(next(out.glob("*_desynth_a_case_gsea_v1.parquet")))
        .select("gene_set", "nes", "fdr")
        .rows()
    )
    assert gsea["UP_SET"][0] > 0 and gsea["UP_SET"][1] < 0.05
    assert gsea["DOWN_SET"][0] < 0 and gsea["DOWN_SET"][1] < 0.05
    ora = pl.read_parquet(next(out.glob("*_desynth_a_case_ora_v1.parquet")))
    up = ora.filter((pl.col("direction") == "up") & (pl.col("gene_set") == "UP_SET"))
    assert up["fdr"][0] < 0.05


def test_inputs_untouched(pipeline):
    root, gathered, _ = pipeline
    now = {p.name: p.read_bytes() for p in (root / "outputs" / "01_gather").iterdir()}
    assert now == gathered


def test_hub_links_resolve_after_move(pipeline, tmp_path):
    root, _, _ = pipeline
    moved = tmp_path / "copied_outputs"
    shutil.copytree(root / "outputs", moved)
    (hub,) = (moved / "09_report").glob("*_desynth_hub_v1.html")
    # Links off the page; the skip link stays on it.
    hrefs = [h for h in re.findall(r'href="([^"]+)"', hub.read_text()) if not h.startswith("#")]
    # Only de and function wrote pages here; setup and gather are listed as not run.
    assert [h.split("/")[1] for h in hrefs] == ["05_de", "06_function"]
    assert "Not run yet: M0 Setup, M1 Gather" in hub.read_text()
    for href in hrefs:
        assert not href.startswith(("http", "/")), href
        assert (hub.parent / href).resolve().is_file(), href


def test_newest_by_date_then_version(tmp_path):
    from ctx.de.run import newest  # noqa: PLC0415

    for day, dataset, version in (("01", "x", 3), ("02", "x", 1), ("02", "x", 2), ("02", "y", 9)):
        (tmp_path / f"2026-09-{day}_gather_{dataset}_counts_v{version}.parquet").write_bytes(b"")
    assert newest(tmp_path, "x", "counts").name == "2026-09-02_gather_x_counts_v2.parquet"


def test_min_count_default_single_source():
    from ctx.de.filters import MIN_COUNT  # noqa: PLC0415

    command = typer.main.get_command(app).commands["de"]
    (option,) = [p for p in command.params if p.opts[0] == "--min-count"]
    assert option.default == MIN_COUNT


def test_lfc_shrink_default_single_source():
    from ctx.de.deseq import LFC_SHRINK  # noqa: PLC0415

    command = typer.main.get_command(app).commands["de"]
    (option,) = [p for p in command.params if p.opts[0] == "--lfc-shrink"]
    assert option.default is LFC_SHRINK
    assert option.secondary_opts == ["--no-lfc-shrink"]


def test_lfc_shrink_reaches_the_fit(monkeypatch):
    # In process, with the profile and the module stubbed: the flag is all that is tested.
    from typer.testing import CliRunner  # noqa: PLC0415

    import ctx.de.run  # noqa: PLC0415
    from ctx import runtime  # noqa: PLC0415

    seen: dict = {}
    monkeypatch.setattr(runtime, "start", lambda profile, root: None)
    monkeypatch.setattr(ctx.de.run, "run", lambda *a, **k: seen.update(k) or [])
    for flag, expected in (([], True), (["--no-lfc-shrink"], False), (["--lfc-shrink"], True)):
        done = CliRunner().invoke(app, ["de", "--dataset", "d", "--dataset-file", "d.toml", *flag])
        assert done.exit_code == 0, done.output
        assert seen["lfc_shrink"] is expected, flag


def test_de_runs_from_wheel(pipeline, tmp_path):
    root, _, _ = pipeline
    fresh = tmp_path / "root"
    shutil.copytree(root / "outputs" / "01_gather", fresh / "outputs" / "01_gather")
    site = build_wheel(tmp_path)
    run_from_wheel(
        site,
        f"""
        import ctx
        from ctx.cli import app
        assert ctx.__file__.startswith({str(site)!r}), ctx.__file__
        app(["--root", {str(fresh)!r}, "de", "--dataset", "desynth",
             "--dataset-file", {str(DE_SPEC)!r}], standalone_mode=False)
        """,
        cwd=tmp_path,
    )
    assert list((fresh / "outputs" / "05_de").glob("*_desynth_overview_v1.html"))


def test_report_decides_sanity(pipeline):
    root, _, _ = pipeline
    out = root / "outputs" / "09_report"
    verdicts = pl.read_parquet(next(out.glob("*_desynth_sanity_v1.parquet")))
    rules = dict(verdicts.select("rule", "rule_verdict").unique().rows())
    assert rules == {"S1": "PASS", "S2": "PASS", "S3": "FAIL"}
    hub = next(out.glob("*_desynth_hub_v1.html")).read_text()
    assert hub.count('class="verdict"') == 4  # S1, both S2 conditions, S3
