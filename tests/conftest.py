"""Session-wide guards."""

import pytest


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    # The template's `test` task maps pytest exit 5 (nothing collected) to success, so a
    # wrong path or an over-eager -k would pass CI silently. That one code fails here instead.
    if exitstatus == pytest.ExitCode.NO_TESTS_COLLECTED:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
