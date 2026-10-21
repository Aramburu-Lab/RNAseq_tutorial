"""PK.5: plant each mutation, run the tests that should catch it, restore the file.

Usage: py.sh python docs/plans/logs/wo13_pk5/mutate.py <basetemp>. A mutation is caught when
a test fails, not when one errors: a fixture that cannot build says nothing about the
mutation. Each file is restored from the bytes read first.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[4]
T = "src/ctx/report/templates/"
ORA, PAGE = T + "_ora.js.j2", T + "_page_controls.js.j2"
F = "src/ctx/function/"
LIVE, RUN = F + "live.py", F + "run.py"
HAND, NODE = "tests/test_handoff_browser.py", "tests/test_ora_page.py"
STATIC = "tests/test_function_words.py"
LINE = '        hline=(-math.log10(ALPHA), figure("ora", "line", alpha=ALPHA)),\n'
MUTATIONS = [
    ("JS twin length -log10 p", ORA, "length: -Math.log10(r.fdr)", "length: -Math.log10(r.p)",
     NODE),
    ("page redraw length -log10 p", PAGE, "return d.length;", "return -Math.log10(d.row.p);",
     HAND),
    ("page redraw drops the line", PAGE, "          { yaxis: Object.assign(div.layout.yaxis",
     "          { shapes: [], yaxis: Object.assign(div.layout.yaxis", HAND),
    ("Python length -log10 p", LIVE, '**{LENGTH: -pl.col("fdr").log10()}',
     '**{LENGTH: -pl.col("p").log10()}', f"{NODE} {STATIC}"),
    ("static length -log10 p", RUN, "    top = _labelled(live.bars(frame))",
     '    top = _labelled(live.bars(frame).with_columns(**{live.LENGTH: -pl.col("p").log10()}))',
     STATIC),
    ("static line left out", RUN, LINE, "", STATIC),
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
                                   f"--basetemp={base / f'm{i}'}", *tests.split()], cwd=W,
                                  capture_output=True, text=True, check=False)  # fmt: skip
        finally:
            path.write_bytes(before)
        last = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else done.stderr[-300:]
        caught = " failed" in last and " error" not in last
        print(f"{'caught' if caught else 'NOT CAUGHT'}: {name}: {last}", flush=True)


if __name__ == "__main__":
    main()
