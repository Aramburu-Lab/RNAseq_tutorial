"""NC's planted mutations of src/ctx/genecard/cache.py: each must fail at least one test.

Counts failed tests from pytest's summary line, never errors (lessons, PK 2026-10-20).
Run through wo13_nc_py.sh; restores the file after every mutation.
"""

import re
import subprocess
import sys
from pathlib import Path

SOURCE = Path("src/ctx/genecard/cache.py")
TESTS = ["tests/test_genecard_absent.py", "tests/test_genecard.py"]
MUTATIONS = {
    "5xx cached": ("if error.code != ABSENT:", "if error.code < 400:"),
    "404 not cached": ("if error.code != ABSENT:", "if True:"),
    "timeout cached": (
        "        except (urllib.error.URLError, OSError) as error:\n"
        "            return Response(url, None, error=str(error))",
        '        except (urllib.error.URLError, OSError):\n            body, status = b"", ABSENT',
    ),
    "offline wording differs": (
        'return Response(url, None, entry["retrieved"], absent(ABSENT, entry["retrieved"]))',
        'return Response(url, None, entry["retrieved"], absent(ABSENT, ""))',
    ),
    "recorded absence asked again before --refresh": (
        'if entry["status"] == ABSENT:\n',
        'if entry["status"] == ABSENT and self.offline:\n',
    ),
    "--refresh deletes the absence": ("shutil.move(folder, target)", "shutil.rmtree(folder)"),
    "release not recorded": ('"release": release, ', ""),
    "status not recorded": ('"status": status, ', '"status": 200, '),
}


def main() -> int:
    """Plant each mutation, run the tests, restore; print one line per mutation."""
    original = SOURCE.read_text()
    survived = 0
    try:
        for name, (old, new) in MUTATIONS.items():
            assert original.count(old) == 1, name
            SOURCE.write_text(original.replace(old, new))
            command = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-n", "auto"]
            done = subprocess.run([*command, *TESTS], capture_output=True, text=True, check=False)
            summary = done.stdout.strip().splitlines()[-1]
            failed = int(m.group(1)) if (m := re.search(r"(\d+) failed", summary)) else 0
            survived += failed == 0
            print(f"{name}: {failed} failed tests | {summary}", flush=True)
    finally:
        SOURCE.write_text(original)
    print(f"caught {len(MUTATIONS) - survived} of {len(MUTATIONS)}")
    return 1 if survived else 0


if __name__ == "__main__":
    sys.exit(main())
