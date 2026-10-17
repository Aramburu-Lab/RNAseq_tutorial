"""K.13: a browser driver past its time limit takes its browser with it.

Killing only the node driver left its Chromium running: four headless Chromiums of 2026-10-15
and 2026-10-16 were still up on 2026-10-17. A stand-in browser that never answers is started
through the driver; past the limit, both are gone. The mutation caught: the driver killed
alone (``proc.kill()`` in place of the session's ``killpg``), which leaves the stand-in up.
"""

import os
import time

import pytest
from test_export_browser import _browse


def test_a_driver_past_its_limit_leaves_no_browser_running(tmp_path, monkeypatch):
    pidfile = tmp_path / "browser.pid"
    fake = tmp_path / "chromium"
    fake.write_text('#!/bin/sh\necho $$ > "$CTX_FAKE_PID"\nexec sleep 300\n')
    fake.chmod(0o755)
    monkeypatch.setenv("CTX_FAKE_PID", str(pidfile))
    with pytest.raises(AssertionError, match="browser driver past 5 s"):
        _browse([{"name": "never", "url": "about:blank"}], tmp_path, chromium=str(fake), timeout=5)
    pid = int(pidfile.read_text())
    for _ in range(50):  # the kill is sent; give the kernel a moment to reap
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return
        time.sleep(0.1)
    os.kill(pid, 9)
    pytest.fail(f"the stand-in browser {pid} outlived its driver")
