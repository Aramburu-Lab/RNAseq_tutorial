"""PK3: plant each mutation, run the tests that should catch it, restore the file.

Usage: python docs/plans/logs/wo13_pk3/mutate.py <basetemp>, with PATH holding main's
environment (node), PYTHONPATH this worktree's src and CTX_CHROMIUM set. A mutation is caught
when a test fails, not when one errors: a fixture that cannot build says nothing about the
mutation. Each file is restored from the bytes read first.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[4]
T = "src/ctx/report/templates/"
ORA, PAGE = T + "_ora.js.j2", T + "_page_controls.js.j2"
F = "src/ctx/function/"
LIVE, RUN, WORDS = F + "live.py", F + "run.py", F + "page.py"
HAND, NODE = "tests/test_handoff_browser.py", "tests/test_ora_page.py"
STATIC = "tests/test_function_words.py"
FILL = ("    return passing.concat(held.filter(function (d) { return !d.passes; })\n"
        "      .slice(0, top - passing.length));")  # fmt: skip
MUTATIONS = [
    ("JS fill dropped", ORA, FILL, "    return passing;", NODE),
    ("JS fill past top", ORA, ".slice(0, top - passing.length)", ".slice(0, top)", NODE),
    ("JS set with no hit drawn", ORA, "return r.hits > 0;", "return r.hits >= 0;", NODE),
    ("JS tie by name reversed", ORA, "(a.name < b.name ? -1 : a.name > b.name ? 1 : 0)",
     "(a.name < b.name ? 1 : a.name > b.name ? -1 : 0)", NODE),
    ("JS families unsorted", ORA,
     "      .sort(function (a, b) { return a[1][0].row.p - b[1][0].row.p; });", ";", NODE),
    ("page fill dropped", ORA, FILL, "    return passing;", HAND),
    ("page grey in the family colour", PAGE, "passes ? colours[f[0]] : colours.failing",
     "colours[f[0]]", HAND),
    ("page length -log10 FDR", PAGE, "return -Math.log10(r.p);", "return -Math.log10(r.fdr);",
     HAND),
    ("page grey sentence left out", PAGE, '(failing ? " " + ctxPick.fill(words.grey',
     '(false ? " " + ctxPick.fill(words.grey', HAND),
    ("Python fill dropped", LIVE,
     "    return pl.concat([passing, held.filter(~pl.col(\"passes\")).head(top - passing.height)])",
     "    return passing", NODE + " " + STATIC),
    ("Python set with no hit drawn", LIVE, 'pl.col("hits") > 0', 'pl.col("hits") >= 0', NODE),
    ("Python families unsorted", LIVE,
     'families = sorted((f for f in families if f.height), key=lambda f: f["p"][0])',
     "families = [f for f in families if f.height]", NODE),
    ("static grey in the up colour", RUN, '"up_failing": None,',
     '"up_failing": ROLES["comparison"],', STATIC),
    ("static length -log10 FDR", RUN, 'neg_log10_p=-pl.col("p").log10()',
     'neg_log10_p=-pl.col("fdr").log10()', STATIC),
    ("static legend says passing", RUN,
     'names[f"{d}_failing"] = fill(PAGE, "live.legend.failing"',
     'names[f"{d}_failing"] = fill(PAGE, "live.legend.passing"', STATIC),
    ("static reading line kept", WORDS, "if failing\n", "if False\n", STATIC),
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
