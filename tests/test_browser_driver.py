"""K.13: a browser driver past its time limit takes its browser with it.

Killing only the node driver left its Chromium running: four headless Chromiums of 2026-10-15
and 2026-10-16 were still up on 2026-10-17. A stand-in browser that never answers is started
through the driver; past the limit, both are gone. The mutation caught: the driver killed
alone (``proc.kill()`` in place of the session's ``killpg``), which leaves the stand-in up.
"""

import os
import re
import time

import pytest
from test_export_browser import CHROMIUM, _browse


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


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_a_helper_that_outlives_the_browser_does_not_hold_the_driver(tmp_path, monkeypatch):
    """H.6: Chromium's helpers and its detached crashpad handlers share its stderr, and on a
    loaded Mac outlived the browser by minutes, while the driver read that stderr and so stayed
    alive after "done", past its 900 s limit. A stand-in helper, a sleep holding the stderr of
    the real browser, outlives it here by 40 s; the driver is done within 20 s anyway. The
    mutation caught: Chromium's stderr left open at exit (no ``proc.stderr.destroy()``)."""
    pidfile = tmp_path / "helper.pid"
    wrapper = tmp_path / "chromium"
    wrapper.write_text(
        f'#!/bin/sh\nsleep 40 &\necho $! > "$CTX_HELPER_PID"\nexec "{CHROMIUM}" "$@"\n'
    )
    wrapper.chmod(0o755)
    monkeypatch.setenv("CTX_HELPER_PID", str(pidfile))
    try:
        got = _browse([{"name": "one", "url": "about:blank", "act": "1 + 1"}], tmp_path,
                      chromium=str(wrapper), timeout=20)  # fmt: skip
        assert got["one"]["acted"] == 2  # noqa: PLR2004 - the act's own answer
    finally:
        os.kill(int(pidfile.read_text()), 9)


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_a_browser_that_dies_mid_task_is_named_at_once(tmp_path):
    """H.6: Chromium gone mid-act, its command is never answered. The driver says the browser
    closed its connection, and in which step, well within its limit. The mutation caught: no
    close handler, which leaves the act waiting with no cause named."""
    wrapper = tmp_path / "chromium"
    wrapper.write_text(f'#!/bin/sh\n(sleep 6; kill -9 $$) &\nexec "{CHROMIUM}" "$@"\n')
    wrapper.chmod(0o755)
    with pytest.raises(AssertionError) as raised:
        _browse([{"name": "wait", "url": "about:blank", "act": "new Promise(() => {})"}],
                tmp_path, chromium=str(wrapper), timeout=40)  # fmt: skip
    told = str(raised.value)
    assert "past 40 s" not in told, told
    assert "the browser closed its DevTools connection" in told and "in wait: act" in told, told


SPIN = """<!doctype html><title>spin</title><script>
function spin() { return Promise.resolve().then(spin); }
addEventListener('load', () => setTimeout(spin, 0));
</script>"""


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
def test_a_page_that_never_yields_is_reported_with_its_frame(tmp_path):
    """H.6: a page in a microtask loop answers no evaluation, and the driver, past its limit,
    names the function, script and line it is in. The mutation caught: no SIGTERM handler, or
    a diagnosis that reads no call frames, which leaves only the trace's step."""
    page = tmp_path / "spin.html"
    page.write_text(SPIN)
    act = "new Promise(r => setTimeout(r, 100))"
    with pytest.raises(AssertionError) as raised:
        _browse([{"name": "spin", "url": page.as_uri(), "act": act}], tmp_path, timeout=15)
    told = str(raised.value)
    assert "hang diagnosis: in spin: act" in told, told
    assert re.search(rf"at spin {re.escape(page.as_uri())}:2:\d+", told), told
