"""wo13-y: plant each mutation, run the tests named to catch it, restore, and log the verdict.

Usage: wo13_y_py.sh python docs/plans/logs/wo13_y_mutations.py <item> [<item> ...]
Each mutation replaces one exact text, found exactly once, and is caught when its tests fail.
Run with PYTHONDONTWRITEBYTECODE=1 (wo13_y_py.sh sets it); __pycache__ is removed after each.
"""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LOG = Path(__file__).with_suffix(".log")
M7 = "tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte"
M7_MADE, M7_ACROSS = f"{M7}[made]", f"{M7}[across_root]"
M9 = "tests/test_figures_m9.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte"

# item: [(name, file, old, new, tests)]
MUTATIONS: dict[str, list[tuple[str, str, str, str, list[str]]]] = {
    "1": [
        (
            "the paired box redrawn without its pairing",
            "src/ctx/figures/draw_m5.py",
            "sum(t.p_adj < ALPHA for t in tests), bool(data.pairs))",
            "sum(t.p_adj < ALPHA for t in tests))",
            ["tests/test_de_box_order.py"],
        ),
        (
            "the running kind naming the function V2 moved",
            "src/ctx/figures/kinds_m6.py",
            '"ctx.function.running:scatter"',
            '"ctx.function.run:_running"',
            ["tests/test_figures.py::test_every_kind_names_inputs_its_drawer_reads"],
        ),
    ],
    "2": [
        (
            "the normals counted as the tumours",
            "src/ctx/proteo/figure_tables.py",
            '(shown["tissue"] == "normal").sum()',
            '(shown["tissue"] == "tumour").sum()',
            [
                "tests/test_figures_m7.py::test_the_normals_class_figure_counts_the_normals_read_back"
            ],
        ),
        (
            "the redraw without the normals' class test",
            "src/ctx/figures/draw_m7.py",
            "        m7.normal_metrics(measured, spec.proteo, ctx.options.method)\n",
            "        pass\n",
            [M7_MADE],
        ),
        (
            "the normals' kind left out",
            "src/ctx/figures/kinds_m7.py",
            '_normals, lambda _: "rho_by_class_normal"',
            '(lambda _: []), lambda _: "rho_by_class_normal"',
            [M7_MADE],
        ),
        (
            "the normals' kind made for every design that correlates",
            "src/ctx/figures/kinds_m7.py",
            '_normals, lambda _: "rho_by_class_normal"',
            '_across, lambda _: "rho_by_class_normal"',
            [M7_ACROSS],
        ),
        (
            "the normals' drawer keeping the tumours' figure",
            "src/ctx/figures/draw_m7.py",
            'return _correlation(inputs, spec, "rho_by_class_normal")',
            'return _correlation(inputs, spec, "rho_by_class")',
            [M7_MADE],
        ),
    ],
    "3": [
        (
            "one figure's panel losing its pixi tab",
            "src/ctx/figures/panel.py",
            '"tabs": [_tab(spec, tab) for tab in TABS],',
            '"tabs": [_tab(spec, tab) for tab in TABS'
            ' if not (tab == "pixi" and "immune" in spec.figure_id)],',
            [M9],
        ),
        (
            "every panel losing its pixi tab",
            "src/ctx/figures/panel.py",
            '"tabs": [_tab(spec, tab) for tab in TABS],',
            '"tabs": [_tab(spec, tab) for tab in TABS[:2]],',
            [M9],
        ),
        (
            "the page's lines counted over every tab, as before",
            "tests/test_figures_m9.py",
            "    return each[TABS[0]]",
            '    return sorted(re.findall(r"ctx figure \\S+ --params", html.unescape(shown)))',
            [M9],
        ),
    ],
}


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
    """Plant every mutation of ``items``, or one as ``item:index``; return how many survived."""
    survived = 0
    with LOG.open("a") as log:
        for asked in items:
            item, _, index = asked.partition(":")
            chosen = MUTATIONS[item][int(index) : int(index) + 1] if index else MUTATIONS[item]
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
                print(line)
                log.write(line + "\n")
            clean = _run(sorted({t for m in chosen for t in m[4]}))
            line = f"item {item}: restored tree, the same tests exit {clean}"
            print(line)
            log.write(line + "\n")
    return survived


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
