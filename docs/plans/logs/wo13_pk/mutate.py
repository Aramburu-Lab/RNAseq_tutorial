"""PK: plant each mutation, run the tests that should catch it, restore the file.

Usage: bash py.sh docs/plans/logs/wo13_pk/mutate.py <basetemp>
A mutation is caught when a test fails, not when one errors: a fixture that cannot build says
nothing about the mutation. Each file is restored from the bytes read first.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[4]
T = "src/ctx/report/templates/"
HAND = "tests/test_handoff_browser.py"
MUTATIONS = [
    ("M6 reveal left out", T + "_page_controls.js.j2", "      ctxLive.reveal(form);\n", "", HAND),
    ("M6 names left out", T + "_page_controls.js.j2",
     '(every.length ? " " + ctxPick.fill(words.named, { names: names }) : "")', '""', HAND),
    ("M6 cap ignored", T + "_page_controls.js.j2", "every.slice(0, cfg.named)", "every.slice(0)",
     HAND),
    ("replaceState unguarded", T + "_page_pick.js.j2",
     "      } catch (err) { /* the address alone is left behind */ }\n",
     "      } finally { /* unguarded */ }\n", HAND),
    ("entry page drops the fragment", T + "entry.html.j2", ' + location.hash); });', "); });",
     HAND),
    ("M9 reveal left out", T + "_page_m9.js.j2", "        ctxLive.reveal(panel);\n", "",
     "tests/test_genecard_browser.py"),
    ("M10 reveal left out", T + "_live.js.j2", "        ctxLive.reveal(form);\n", "",
     "tests/test_clinical_live.py"),
]  # fmt: skip


def main() -> None:
    """Run every mutation and print whether its tests caught it."""
    base = Path(sys.argv[1])
    base.mkdir(parents=True, exist_ok=True)
    for i, (name, rel, old, new, tests) in enumerate(MUTATIONS):
        path = W / rel
        before = path.read_bytes()
        text = before.decode()
        assert text.count(old) == 1, (name, text.count(old))
        path.write_text(text.replace(old, new))
        try:
            done = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                                   f"--basetemp={base / f'm{i}'}", tests], cwd=W,
                                  capture_output=True, text=True, check=False)  # fmt: skip
        finally:
            path.write_bytes(before)
        last = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else done.stderr[-300:]
        caught = " failed" in last and " error" not in last
        print(f"{'caught' if caught else 'NOT CAUGHT'}: {name}: {last}", flush=True)


if __name__ == "__main__":
    main()
