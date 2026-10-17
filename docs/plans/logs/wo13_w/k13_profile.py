"""Package W, K.13: the files one headless Chromium launch writes into its profile.

Starts Chromium with tests/browser.mjs's flags (and, with ``extra``, more), waits for its
DevTools line, lets it idle a few seconds, kills it, and counts the files and bytes left in
the profile, so a file operation's cost can be multiplied by how many a launch makes.
Run: python -I k13_profile.py <scratch folder> <label> [extra flags ...]
"""

import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

FLAGS = ["--headless=new", "--remote-debugging-port=0", "--no-first-run",
         "--no-default-browser-check", "--disable-gpu", "--disable-background-networking",
         "--disable-component-update", "--disable-sync"]  # fmt: skip
CHROMIUM = os.environ.get("CTX_CHROMIUM", "/opt/homebrew/bin/chromium")


def main() -> None:
    """One launch: the seconds to its DevTools line, then the profile's file count and bytes."""
    profile = Path(tempfile.mkdtemp(prefix="w_k13_profile_", dir=sys.argv[1]))
    start = time.perf_counter()
    proc = subprocess.Popen([CHROMIUM, *FLAGS, *sys.argv[3:], f"--user-data-dir={profile}",
                             "about:blank"], stderr=subprocess.PIPE, text=True)  # fmt: skip
    assert proc.stderr is not None
    for line in proc.stderr:
        if "DevTools listening" in line:
            break
    ready = time.perf_counter() - start
    time.sleep(5)
    proc.kill()
    proc.wait()
    files = [p for p in profile.rglob("*") if p.is_file()]
    size = sum(p.stat().st_size for p in files)
    print(f"{sys.argv[2]}\tready {ready:.1f} s\t{len(files)} files\t{size:,} B", flush=True)


if __name__ == "__main__":
    main()
