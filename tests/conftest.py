"""Session-wide guards, and the needs_star rule (R6.2, D-036)."""

import os
import platform
import shutil
import subprocess
import sys
from collections.abc import Callable
from datetime import date
from pathlib import Path

import pytest

# The needs_star hook is tested by running a one-test file through it (test_process_align).
pytest_plugins = ("pytester",)


def tested_tree_problem(package_file: str, checkout: Path) -> str | None:
    """Why a test run would test another checkout's code, or None when it tests this one."""
    src = (checkout / "src").resolve()
    # A copy of this file outside a checkout (pytester's runs) has no src to compare with.
    if not (src / "ctx").is_dir() or Path(package_file).resolve().is_relative_to(src):
        return None
    return (f"ctx imports from {Path(package_file).parent}, not from {src}: the shared editable "
            "install was re-pointed by `pixi run` or `pixi install` in another checkout. Run "
            "`pixi install` in this checkout, or set PYTHONPATH to its src.")  # fmt: skip


def pytest_configure(config: pytest.Config) -> None:
    # One environment serves several worktrees; on 2026-10-01 it pointed at a worktree, and
    # the main tree's tests ran that worktree's code. A run of the wrong code stops here.
    import ctx  # noqa: PLC0415 - after pytest's own path setup

    problem = tested_tree_problem(ctx.__file__, Path(__file__).resolve().parents[1])
    if problem:
        raise pytest.UsageError(problem)
    if not hasattr(config, "workerinput"):
        config.stash[REPO_BEFORE] = repo_state(Path(__file__).resolve().parents[1])


REPO_BEFORE = pytest.StashKey[str | None]()


def repo_state(checkout: Path) -> str | None:
    """``git status`` of the checkout, untracked files listed one by one; None outside git.

    Every worker shares the checkout: a test that writes into it is caught by comparing this
    before and after the session.
    """
    done = subprocess.run(["git", "-C", str(checkout), "status", "--porcelain",
                           "--untracked-files=all"], capture_output=True, text=True,
                          check=False)  # fmt: skip
    return done.stdout if done.returncode == 0 else None


# Optional: a run without xdist (-p no:xdist) has no such hook to validate it against.
@pytest.hookimpl(optionalhook=True)
def pytest_xdist_auto_num_workers(config: pytest.Config) -> int:
    """Leave two cores free under `-n auto`: with every core busy, Defender's scan of the
    browser tests' profile files made a Chromium launch take 44 s against 0.7 s (K.13; Rytis
    2026-10-17, K.15)."""
    return max(1, (os.cpu_count() or 1) - 2)


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    # The template's `test` task maps pytest exit 5 (nothing collected) to success, so a
    # wrong path or an over-eager -k would pass CI silently. That one code fails here instead.
    if exitstatus == pytest.ExitCode.NO_TESTS_COLLECTED:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
    # Under xdist the controller, which ran no test, checks what its workers built.
    if hasattr(session.config, "workerinput"):
        return
    if repeated := built_twice(session.config._tmp_path_factory.getbasetemp()):  # noqa: SLF001
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
        print(f"\nsession runs built more than once: {repeated}", file=sys.stderr)
    before = session.config.stash.get(REPO_BEFORE, None)
    after = repo_state(Path(__file__).resolve().parents[1])
    if before is not None and after != before:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
        changed = set((after or "").splitlines()) ^ set(before.splitlines())
        print(f"\nthe test run changed the checkout: {sorted(changed)}", file=sys.stderr)


BUILT_LOG = "session_runs.log"


def built_twice(shared: Path) -> dict[str, int]:
    """Session runs logged by ``once`` more than once in a run's shared folder."""
    log = shared / BUILT_LOG
    names = log.read_text().split() if log.is_file() else []
    return {name: names.count(name) for name in sorted(set(names)) if names.count(name) > 1}


def _shared(factory: pytest.TempPathFactory) -> Path:
    """The folder every worker of this run sees: the controller's basetemp under xdist, whose
    workers each take a ``popen-gwN`` folder inside it; this process's basetemp otherwise."""
    base = factory.getbasetemp()
    return base.parent if os.environ.get("PYTEST_XDIST_WORKER") else base


def _seal(folder: Path) -> None:
    """Make a session run read-only: workers read it at once, so a test writing into it would
    change what another test is reading. Such a write now fails where it happens."""
    for path in [*folder.rglob("*"), folder]:
        if not path.is_symlink():
            path.chmod(path.stat().st_mode & ~0o222)


def writable_copy(source: Path, target: Path) -> Path:
    """A copy of a sealed session run that a test may change."""
    shutil.copytree(source, target, symlinks=True)
    for path in [*target.rglob("*"), target]:
        if not path.is_symlink():
            path.chmod(path.stat().st_mode | 0o200)
    return target


def once[T](factory: pytest.TempPathFactory, name: str, build: Callable[[Path], T]) -> T:
    """``build(folder)`` run once per test run, however many xdist workers ask for it.

    The first worker to take the lock builds into a folder of the run's shared directory,
    seals it and pickles the value; the others wait on the lock and unpickle it. Every build
    is logged, and the session fails if a name is logged twice (``built_twice``).
    """
    import fcntl  # noqa: PLC0415 - POSIX; the suite runs on macOS and Linux
    import pickle  # noqa: PLC0415
    import tempfile  # noqa: PLC0415

    shared = _shared(factory)
    value_file = shared / f"{name}.pickle"
    with (shared / f"{name}.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not value_file.is_file():
            folder = Path(tempfile.mkdtemp(prefix=f"{name}_", dir=shared))
            value = build(folder)
            _seal(folder)
            with (shared / BUILT_LOG).open("a") as log:
                log.write(f"{name}\n")
            value_file.write_bytes(pickle.dumps(value))
    return pickle.loads(value_file.read_bytes())  # noqa: S301 - written by this run


def star_may_skip() -> bool:
    """R6.2 (Rytis, 2026-10-01): a STAR test may skip only on darwin-arm64 with CI unset.

    The pinned osx-arm64 STAR cannot map (D-036); on Linux, in CI or under Rosetta a failing
    probe is a fault to see, not a platform to forgive.
    """
    return (platform.system() == "Darwin" and platform.machine() == "arm64"
            and "CI" not in os.environ)  # fmt: skip


def star_outcome(failure: str | None) -> tuple[str, str] | None:
    """What a needs_star test does: run (None), or ("skip" | "fail", why).

    Args:
        failure: The STAR probe's message, None when it mapped its reads.
    """
    if failure is None:
        return None
    if star_may_skip():
        return "skip", f"STAR cannot map on darwin-arm64 with CI unset (R6.2, D-036): {failure}"
    return "fail", f"STAR probe failed; only darwin-arm64 with CI unset may skip (R6.2): {failure}"


def pytest_runtest_setup(item: pytest.Item) -> None:
    if item.get_closest_marker("needs_star") is None:
        return
    from ctx import tools  # noqa: PLC0415 - only a STAR test pays for the analysis imports
    from ctx.process import align  # noqa: PLC0415

    try:
        align.star_probe()
        failure = None
    except tools.ToolError as error:
        failure = str(error)
    outcome = star_outcome(failure)
    if outcome is None:
        return
    kind, why = outcome
    if kind == "skip":
        pytest.skip(why)
    pytest.fail(why, pytrace=False)


def threads_unlike_default() -> int:
    """A thread count no library would pick by itself and OpenBLAS can honour.

    One below the cores this process may run on: OpenBLAS clamps OPENBLAS_NUM_THREADS to those
    cores, so asking for more read back as the core count (3 asked, 2 seen on a 2-core runner,
    run 37807649888), and polars sizes its pool to them by default, so asking for them proves
    nothing. Skips on one core, where no count is both.
    """
    cores = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    if not cores or cores < 2:
        pytest.skip("one core: no thread count differs from the default and is honoured")
    return cores - 1


CTX = Path(sys.executable).parent / "ctx"


def ctx_run(*args: str, cwd: Path, env: dict[str, str] | None = None
            ) -> subprocess.CompletedProcess[str]:  # fmt: skip
    """``ctx`` in a subprocess, SLURM_ and CTX_ variables cleared unless ``env`` sets them."""
    clean = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    return subprocess.run([str(CTX), *args], cwd=cwd, env={**clean, **(env or {})},
                          capture_output=True, text=True, check=False)  # fmt: skip


# The session run's cap: no fixture run holds more reads, so every read is kept and the
# planted truths hold, while the cap's note and its realised counts are on the page.
CAP = "1000"
# The Salmon twin's off-contig counts for se_e: Salmon puts 0.4 of its reads on the slice, below
# the share assigned here, as on a real partial reference (D-039).
SE_E_ELSEWHERE = 500


@pytest.fixture(scope="session")
def process_root(tmp_path_factory):
    """`ctx process --aligner hisat2` on the chr22_se fixture, every setting and a Salmon twin.

    One run for the whole session: tests/test_process_run.py reads it, and test_cli reruns a
    copy of it (``process_run``), so the alignment is paid for once (M3.14).
    """

    def build(root: Path) -> tuple[Path, Path, subprocess.CompletedProcess[str]]:
        from process_root import (  # noqa: PLC0415
            COMPARE,
            DATASET,
            SLICE,
            dataset_file,
            plant,
            plant_compare,
        )

        plant(root)
        plant_compare(root, off={"se_e": SE_E_ELSEWHERE})
        spec = dataset_file(root, root)
        args = ("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                str(spec), "--genome", str(SLICE), "--aligner", "hisat2", "--settings",
                "cut,untrimmed", "--compare-dataset", COMPARE, "--reads-per-run", CAP)  # fmt: skip
        done = ctx_run(*args, cwd=root)
        assert done.returncode == 0, done.stderr[-3000:]
        return root, spec, done

    return once(tmp_path_factory, "process", build)


@pytest.fixture(scope="session")
def process_run(process_rerun):
    """(root, overview html) of ``process_rerun``: test_cli's card test reads its page."""
    root, html, _ = process_rerun
    return root, html


@pytest.fixture(scope="session")
def process_rerun(process_root, tmp_path_factory):
    """A copy of ``process_root`` rerun at the defaults but the aligner and the cap: (root, html).

    Its dataset file now declares itself a subset (D-050). The tool folders' modification
    times before the rerun ride along as the third value, so a test can see what was reused.
    """

    def build(folder: Path) -> tuple[Path, str, dict[str, int]]:
        from process_root import (  # noqa: PLC0415
            DATASET,
            PARENT,
            PARENT_LABEL,
            SLICE,
            dataset_file,
        )

        root = writable_copy(process_root[0], folder / "root")
        tools = root / "outputs" / "03_process" / "tools"
        before = {str(p.relative_to(tools)): p.stat().st_mtime_ns for p in tools.rglob("*")
                  if p.is_file()}  # fmt: skip
        dataset_file(root, root, PARENT, label=PARENT_LABEL)
        copied = dataset_file(root, root, subset_of=PARENT)
        done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                       str(copied), "--genome", str(SLICE), "--aligner", "hisat2",
                       "--reads-per-run", CAP, cwd=root)  # fmt: skip
        assert done.returncode == 0, done.stderr[-3000:]
        (page,) = (root / "outputs" / "03_process").glob(f"*_process_{DATASET}_overview_v1.html")
        return root, page.read_text(), before

    return once(tmp_path_factory, "process_again", build)


@pytest.fixture(scope="session")
def quant_built(tmp_path_factory):
    """`ctx quant` on the quant fixture, once per session: (root, dataset file, files written)."""

    def build(root: Path) -> tuple[Path, Path, list[Path]]:
        from quant_root import DATASET, RNASEQ, SALMON, quant_root  # noqa: PLC0415

        spec = quant_root(root)
        done = ctx_run("--root", str(root), "quant", "--dataset", DATASET, "--dataset-file",
                       str(spec), "--rnaseq-dir", str(RNASEQ), "--salmon-dir", str(SALMON),
                       cwd=root)  # fmt: skip
        assert done.returncode == 0, done.stderr[-3000:]
        return root, spec, [Path(p) for p in done.stdout.split()]

    return once(tmp_path_factory, "quant", build)


@pytest.fixture(scope="session")
def quant_run(quant_built):
    """(root, overview html) of ``quant_built``: test_cli's card test reads its page."""
    root, _, written = quant_built
    return root, written[0].read_text()


# Every module's end-to-end run, once per session: test_cli reads their pages for the cards,
# test_outputs_contract every file they wrote (M3.14, V4).
ROOT = Path(__file__).resolve().parents[1]
FETCHNGS_DOWNLOAD = ROOT / "tests" / "fixtures" / "fetchngs" / "download"
DE_SPEC = ROOT / "tests" / "fixtures" / "datasets" / "de_synthetic.toml"


def setup_page(root: Path) -> str:
    """The environment page ``ctx setup`` wrote under ``root``."""
    (page,) = (root / "outputs" / "00_setup").glob("*_setup_environment_v1.html")
    return page.read_text()


def module_page(root: Path, number_module: str, what: str) -> str:
    """The one ``*_<what>_v1.html`` page in ``outputs/<number_module>``."""
    (page,) = (root / "outputs" / number_module).glob(f"*_{what}_v1.html")
    return page.read_text()


def _built_wheel(folder: Path) -> Path:
    from test_packaging import build_wheel  # noqa: PLC0415

    return build_wheel(folder)


@pytest.fixture(scope="session")
def wheel_site(tmp_path_factory):
    """The wheel built and unzipped once; every test that runs ctx from it only reads it.

    test_io builds its own, inside an unrelated git repository, which is what it tests.
    """
    return once(tmp_path_factory, "wheel", _built_wheel)


def _ran(root: Path, *args: str) -> None:
    done = ctx_run("--root", str(root), *args, cwd=root)
    assert done.returncode == 0, (args[0], done.stderr[-3000:])


def _setup(root: Path) -> tuple[Path, str]:
    _ran(root, "setup")
    return root, setup_page(root)


@pytest.fixture(scope="session")
def setup_run(tmp_path_factory):
    """(root, page) of ``ctx setup``."""
    return once(tmp_path_factory, "root", _setup)


def _gather(root: Path) -> tuple[Path, str]:
    _ran(root, "gather", "--dataset", "nfcore_test", "--fetchngs-dir", str(FETCHNGS_DOWNLOAD))
    return root, module_page(root, "01_gather", "gather_nfcore_test_overview")


@pytest.fixture(scope="session")
def gather_run(tmp_path_factory):
    """(root, overview page) of ``ctx gather`` on the fetchngs download fixture."""
    return once(tmp_path_factory, "gather_root", _gather)


@pytest.fixture(scope="session")
def pipeline(tmp_path_factory, wheel_site):
    """de, function and report on one root, from gather outputs written by gather's writer:
    (root, the gather files' bytes before, the ctx package de imported).

    de runs from the wheel, so its one fit also proves the wheel runs it; it was fitted twice,
    80 s of CI run 37898263890's 1,778 s for the second (M3.14).
    """

    def build(root: Path) -> tuple[Path, dict[str, bytes], Path]:
        from synth import (  # noqa: PLC0415
            DE_GROUPS,
            design_counts,
            write_gather_outputs,
            write_gene_sets,
        )
        from test_packaging import run_from_wheel  # noqa: PLC0415

        from ctx.dataset import load  # noqa: PLC0415

        counts = design_counts(DE_GROUPS, batch_lfc=2.0)
        write_gather_outputs(root, "desynth", counts, date.today())
        write_gene_sets(root / load(DE_SPEC).resource("sets").dest, synthetic=counts)
        gathered = {p.name: p.read_bytes() for p in (root / "outputs" / "01_gather").iterdir()}
        args = ["--root", str(root), "de", "--dataset", "desynth", "--dataset-file", str(DE_SPEC)]
        imported = run_from_wheel(wheel_site, f"""
            import ctx
            from ctx.cli import app
            print(ctx.__file__)
            app({args!r}, standalone_mode=False)
            """, cwd=root).splitlines()[0]  # fmt: skip
        for command in ("function", "report"):
            _ran(root, command, "--dataset", "desynth", "--dataset-file", str(DE_SPEC))
        return root, gathered, Path(imported)

    return once(tmp_path_factory, "pipeline", build)


@pytest.fixture(scope="session")
def de_run(pipeline):
    """(root, overview page) of the pipeline's ``ctx de``."""
    return pipeline[0], module_page(pipeline[0], "05_de", "desynth_overview")


@pytest.fixture(scope="session")
def function_run(pipeline):
    """(root, overview page) of the pipeline's ``ctx function``."""
    return pipeline[0], module_page(pipeline[0], "06_function", "desynth_overview")


@pytest.fixture(scope="session")
def report_run(pipeline):
    """(root, hub page) of the pipeline's ``ctx report``."""
    return pipeline[0], module_page(pipeline[0], "09_report", "desynth_hub")


def _qc(root: Path) -> tuple[Path, str]:
    from test_qc_run import plant  # noqa: PLC0415

    plant(root)
    _ran(root, "qc", "--dataset", "desynth", "--dataset-file", str(DE_SPEC))
    return root, module_page(root, "02_qc", "desynth_overview")


@pytest.fixture(scope="session")
def qc_run(tmp_path_factory):
    """(root, overview page) of ``ctx qc`` on the planted gather outputs of the qc tests."""
    return once(tmp_path_factory, "qc_cli", _qc)


def _clinical(folder: Path) -> tuple[Path, str]:
    from clinical_root import make  # noqa: PLC0415

    root, spec = make(folder)
    _ran(root, "clinical", "--dataset", "fake_brca", "--dataset-file", str(spec))
    return root, module_page(root, "12_clinical", "fake_brca_overview")


@pytest.fixture(scope="session")
def clinical_run(tmp_path_factory):
    """(root, overview page) of ``ctx clinical`` on the made-up GDC root of the M10 tests."""
    return once(tmp_path_factory, "clinical_cli", _clinical)


def _proteo(folder: Path) -> tuple[Path, str]:
    from proteo_root import cross  # noqa: PLC0415

    root = cross(folder)
    return root, module_page(root, "07_proteo", "fake_cells_overview")


@pytest.fixture(scope="session")
def proteo_run(tmp_path_factory):
    """(root, overview page) of the made-up cross-dataset root of the M7 tests."""
    return once(tmp_path_factory, "proteo_cli", _proteo)


def _genecard(root: Path) -> tuple[Path, str]:
    from genecard_root import SLUG, make  # noqa: PLC0415

    spec = make(root)
    _ran(root, "genecard", "--dataset", SLUG, "--dataset-file", str(spec), "--genes", "NPTX2",
         "--offline")  # fmt: skip
    return root, module_page(root, "11_genecard", f"{SLUG}_overview")


@pytest.fixture(scope="session")
def genecard_run(tmp_path_factory):
    """(root, overview page) of ``ctx genecard --offline`` on the made-up root of the M9 tests."""
    return once(tmp_path_factory, "genecard_cli", _genecard)


def _ncrna(folder: Path) -> tuple[Path, str]:
    from ncrna_root import DATASET, make  # noqa: PLC0415

    root, spec = make(folder, date.today())
    _ran(root, "ncrna", "--dataset", DATASET, "--dataset-file", str(spec))
    return root, module_page(root, "08_ncrna", f"{DATASET}_overview")


@pytest.fixture(scope="session")
def ncrna_run(tmp_path_factory):
    """(root, overview page) of ``ctx ncrna`` on the made-up root of the M8 tests."""
    return once(tmp_path_factory, "ncrna_cli", _ncrna)


def _guide(root: Path) -> tuple[Path, str]:
    _ran(root, "guide")
    return root, module_page(root, "09_report", "report_guide")


@pytest.fixture(scope="session")
def guide_run(tmp_path_factory):
    """(root, page) of ``ctx guide``."""
    return once(tmp_path_factory, "guide_cli", _guide)


@pytest.fixture(scope="session")
def run_run(guide_run):
    """`ctx run` writes no page of its own: its cards are on the guide, which teaches it."""
    return guide_run
