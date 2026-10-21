"""conftest.once: a session run is built once across xdist workers and is read-only."""

import os

import pytest
from conftest import BUILT_LOG, built_twice, module_group, writable_copy


@pytest.mark.skipif(os.name == "posix" and os.geteuid() == 0, reason="root ignores modes")
def test_a_session_run_refuses_a_write(setup_run, tmp_path):
    """Workers read a session run at once, so a test that writes into one fails where it
    writes; a writable copy is the way to change one.

    Mutation caught: ``_seal`` dropped from ``once`` (the write succeeds)."""
    root, _ = setup_run
    with pytest.raises(PermissionError):
        (root / "outputs" / "planted.txt").write_text("x")
    copy = writable_copy(root, tmp_path / "copy")
    (copy / "outputs" / "planted.txt").write_text("x")


def test_a_run_built_twice_is_named(tmp_path):
    """The log ``once`` appends to, read as the controller reads it at the end of the session.

    Mutation caught: a name logged twice not reported (the count compared with 2, not 1)."""
    assert built_twice(tmp_path) == {}
    (tmp_path / BUILT_LOG).write_text("process\nquant\nprocess\n")
    assert built_twice(tmp_path) == {"process": 2}


@pytest.fixture(scope="module")
def module_scoped():
    """A module fixture: built once per worker that runs a test asking for it."""
    return object()


def test_a_test_on_a_module_fixture_is_scheduled_with_its_module(request, module_scoped):
    """Under ``--dist loadgroup`` (the test task) such a test is in its module's xdist group,
    which xdist writes into its node ID.

    Mutation caught: the group added only where the dist option says loadgroup, which a
    worker's never does (the hand-off module ran on 9 workers, a Chromium run on each)."""
    assert module_group(request.node) == __name__
    if getattr(request.config.option, "loadgroup", False):
        assert request.node.nodeid.endswith(f"@{__name__}"), request.node.nodeid


def test_a_test_on_session_and_function_fixtures_goes_to_any_worker(request, tmp_path):
    """Mutation caught: every test grouped by its module, which is loadscope again."""
    assert module_group(request.node) is None
    assert "@" not in request.node.nodeid
