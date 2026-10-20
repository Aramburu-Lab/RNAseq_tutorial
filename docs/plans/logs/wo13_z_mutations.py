"""wo13-z: plant each mutation, run the tests named to catch it, restore, and log the verdict.

Usage: wo13_z_py.sh python docs/plans/logs/wo13_z_mutations.py <item> [<item> ...]
Each mutation replaces one exact text, found exactly once, and is caught when its tests fail.
Run with PYTHONDONTWRITEBYTECODE=1 (wo13_z_py.sh sets it); __pycache__ is removed after each.
"""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LOG = Path(__file__).with_suffix(".log")
BR = "tests/test_reproduce_browser.py"
EACH = f"{BR}::test_each_tab_copies_and_downloads_its_whole_script"
CAP = "tests/test_page_cap.py"
JS = "src/ctx/report/templates/_page.js.j2"
REPRO = "src/ctx/report/templates/_reproduce.html.j2"

# item: [(name, file, old, new, tests)]
MUTATIONS: dict[str, list[tuple[str, str, str, str, list[str]]]] = {
    "panel": [
        ("the page's join without the blank line", JS,
         'code.textContent = code.textContent + "\\n" + steps.textContent;',
         "code.textContent = code.textContent + steps.textContent;", [EACH]),
        ("the download without the script head", JS,
         "new TextEncoder().encode(head + code.textContent)",
         "new TextEncoder().encode(code.textContent)", [EACH]),
        ("the shared block left shown", JS, "box.hidden = true;", "box.hidden = false;", [EACH]),
        ("every tab's download holding the first tab's block", JS,
         "panel.querySelector('a[data-script=\"' + code.id + '\"]').href =",
         "panel.querySelector('a[data-script]').href =", [EACH]),
        ("body without the blank line after step 1", "src/ctx/figures/panel.py",
         'return f"{head(spec, tab)}\\n{shared(spec)}"',
         'return f"{head(spec, tab)}{shared(spec)}"',
         [f"{BR}::test_body_and_script_are_main_s_text"]),
        ("each tab printing the shared steps again", REPRO,
         ">{{ t.head }}</code>", ">{{ t.head }}{{ r.shared }}</code>",
         [f"{BR}::test_the_shared_steps_are_printed_once"]),
        ("the download link shown without JavaScript", REPRO,
         '<p class="js-only"><a download=', "<p><a download=",
         [f"{BR}::test_without_javascript_every_step_reads_once"]),
        ("the image proof's join without the blank line", "containers/prove_figures.py",
         'pasted = f"{html.unescape(step_1)}\\n{html.unescape(steps)}"',
         'pasted = f"{html.unescape(step_1)}{html.unescape(steps)}"',
         ["tests/test_figure_panel.py::test_the_image_proof_reads_the_docker_tab_the_page_prints"]),
    ],
    "cap": [
        ("M2's view found by the stem suffix as before", "src/ctx/qc/flag_map.py",
         "BAND_MAP.search(item.stem)", 'item.stem.endswith("_band_map")',
         [f"{CAP}::test_m2_drops_its_flag_view_under_a_cap_one_byte_below_its_page"]),
        ("M9's cards fitted once, before their panels", "src/ctx/genecard/page.py",
         "if (fits := cap.kept(pages, kept, least=1)) == kept:", "if (fits := kept) == kept:",
         [f"{CAP}::test_m9_decides_its_cards_on_the_page_with_its_panels"]),
        ("M9's refit keeping the files of a card that left", "src/ctx/genecard/page.py",
         "            path.unlink(missing_ok=True)\n", "            pass\n",
         [f"{CAP}::test_m9_decides_its_cards_on_the_page_with_its_panels"]),
        ("builder.write holding its own 15 MB", "src/ctx/report/builder.py",
         "    if cap.over(html) > 0:", "    if len(html.encode()) > 15 * 10**6:",
         [f"{CAP}::test_one_cap_for_every_reader"]),
        ("an unmeasured part freeing nothing", "src/ctx/report/cap.py",
         "if count - least < len(sizes) else excess", "if count - least < len(sizes) else 0",
         [f"{CAP}::test_parts_without_sizes_leave_one_at_a_time"]),
    ],
}  # fmt: skip


def _run(tests: list[str]) -> int:
    """The tests' exit code, run serially in this worktree."""
    argv = [sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", *tests]
    done = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, check=False)
    return done.returncode


def _clean() -> None:
    """Remove bytecode a mutated run may have left (lessons, 2026-09-30)."""
    for cache in [*ROOT.joinpath("src").rglob("__pycache__"),
                  *ROOT.joinpath("tests").rglob("__pycache__")]:  # fmt: skip
        shutil.rmtree(cache, ignore_errors=True)


def main(items: list[str]) -> int:
    """Plant every mutation of ``items``; return how many survived."""
    survived = 0
    with LOG.open("a") as log:
        for item in items:
            chosen = MUTATIONS[item]
            for name, file, old, new, tests in chosen:
                path = ROOT / file
                source = path.read_text()
                assert source.count(old) == 1, (name, source.count(old))
                path.write_text(source.replace(old, new))
                try:
                    code = _run(tests)
                finally:
                    path.write_text(source)
                    _clean()
                verdict = "caught" if code != 0 else "SURVIVED"
                survived += code == 0
                line = f"item {item}: {verdict} (exit {code}): {name} [{' '.join(tests)}]"
                print(line, flush=True)
                log.write(line + "\n")
            clean = _run(sorted({t for m in chosen for t in m[4]}))
            line = f"item {item}: restored tree, the same tests exit {clean}"
            print(line, flush=True)
            log.write(line + "\n")
    return survived


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
