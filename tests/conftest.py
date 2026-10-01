"""Session-wide guards, and the needs_star rule (R6.2, D-036)."""

import os
import platform

import pytest

# The needs_star hook is tested by running a one-test file through it (test_process_align).
pytest_plugins = ("pytester",)


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
