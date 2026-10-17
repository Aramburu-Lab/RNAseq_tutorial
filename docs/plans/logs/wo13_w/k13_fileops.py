"""Package W, K.13: how long the file operations a browser launch is made of take, and CPU work.

Times, in one process: CPU work (a fixed loop), writing then reading 500 new small files,
deleting them, and 30 starts of /usr/bin/true. Run once quiet and once with every core busy:
CPU work slows by the share of the cores the process gets, so a file operation that slows far
more than that waits on something other than this process's own CPU.
Run: python -I k13_fileops.py <scratch folder> <label>
"""

import subprocess
import sys
import tempfile
import time
from pathlib import Path


def timed(what, work):
    """Seconds ``work`` takes, printed with ``what``."""
    start = time.perf_counter()
    work()
    took = time.perf_counter() - start
    print(f"{sys.argv[2]}\t{what}\t{took:.3f} s", flush=True)
    return took


def cpu():
    """A fixed amount of arithmetic."""
    total = 0
    for i in range(3_000_000):
        total += i * i % 7
    return total


def main():
    """Each measurement once, in a fresh folder."""
    folder = Path(tempfile.mkdtemp(prefix="w_k13_fileops_", dir=sys.argv[1]))
    files = [folder / f"f{i}.txt" for i in range(500)]
    timed("cpu loop", cpu)
    timed("write 500 new files", lambda: [f.write_text("x" * 64) for f in files])
    timed("read 500 files", lambda: [f.read_text() for f in files])
    timed("delete 500 files", lambda: [f.unlink() for f in files])
    timed("30 starts of /usr/bin/true", lambda: [subprocess.run(["/usr/bin/true"], check=True)
                                                for _ in range(30)])  # fmt: skip
    folder.rmdir()


if __name__ == "__main__":
    main()
