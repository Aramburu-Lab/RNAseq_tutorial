"""Session-wide guards, and the needs_star rule (R6.2, D-036)."""

import os
import platform
import shutil
import subprocess
import sys
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


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    # The template's `test` task maps pytest exit 5 (nothing collected) to success, so a
    # wrong path or an over-eager -k would pass CI silently. That one code fails here instead.
    if exitstatus == pytest.ExitCode.NO_TESTS_COLLECTED:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED


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
    from process_root import (  # noqa: PLC0415
        COMPARE,
        DATASET,
        SLICE,
        dataset_file,
        plant,
        plant_compare,
    )

    root = tmp_path_factory.mktemp("process")
    plant(root)
    plant_compare(root, off={"se_e": SE_E_ELSEWHERE})
    spec = dataset_file(root, root)
    args = ("--root", str(root), "process", "--dataset", DATASET, "--dataset-file", str(spec),
            "--genome", str(SLICE), "--aligner", "hisat2", "--settings", "cut,untrimmed",
            "--compare-dataset", COMPARE, "--reads-per-run", CAP)  # fmt: skip
    done = ctx_run(*args, cwd=root)
    assert done.returncode == 0, done.stderr[-3000:]
    return root, spec, done


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
    from process_root import DATASET, PARENT, PARENT_LABEL, SLICE, dataset_file  # noqa: PLC0415

    source, _, _ = process_root
    root = tmp_path_factory.mktemp("process_again") / "root"
    shutil.copytree(source, root)
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


@pytest.fixture(scope="session")
def quant_built(tmp_path_factory):
    """`ctx quant` on the quant fixture, once per session: (root, dataset file, files written)."""
    from quant_root import DATASET, RNASEQ, SALMON, quant_root  # noqa: PLC0415

    root = tmp_path_factory.mktemp("quant")
    spec = quant_root(root)
    done = ctx_run("--root", str(root), "quant", "--dataset", DATASET, "--dataset-file",
                   str(spec), "--rnaseq-dir", str(RNASEQ), "--salmon-dir", str(SALMON),
                   cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return root, spec, [Path(p) for p in done.stdout.split()]


@pytest.fixture(scope="session")
def quant_run(quant_built):
    """(root, overview html) of ``quant_built``: test_cli's card test reads its page."""
    root, _, written = quant_built
    return root, written[0].read_text()
