"""conftest.once: a session run is built once across xdist workers and is read-only."""

import os

import pytest
from conftest import BUILT_LOG, built_twice, writable_copy


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
