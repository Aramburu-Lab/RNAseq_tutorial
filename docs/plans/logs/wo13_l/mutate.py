"""Package wo13-l: plant each mutation, run the tests that should catch it, restore the file.

Usage: <main env python> docs/plans/logs/wo13_l/env.py python docs/plans/logs/wo13_l/mutate.py
A mutation is caught when its tests fail; one whose tests pass, or skip, survives.
"""

import subprocess
import sys
from pathlib import Path

WT = Path(__file__).absolute().parents[4]
CSS = "src/ctx/report/templates/_page.css.j2"
BAR = "src/ctx/viz/bar.py"
LEGEND = "tests/test_legend_browser.py::test_long_bar_labels_wrap_inside_the_plot_on_a_phone"
WIDTH = "tests/test_page_width_browser.py::test_a_long_word_in_page_text_wraps_inside_its_box"
TN = "tests/test_proteo_tumour_normal.py"
MUTATIONS = [
    (
        "css break-word removed",
        CSS,
        "antialiased;overflow-wrap:break-word}",
        "antialiased}",
        [WIDTH],
    ),
    (
        "no narrow bar labels",
        BAR,
        "fixed = self.horizontal and self.live is None and",
        "fixed = False and",
        [LEGEND],
    ),
    (
        "bar labels wrapped to the whole phone plot",
        BAR,
        "    room = PHONE_PLOT_PX / 2",
        "    room = PHONE_PLOT_PX",
        [LEGEND],
    ),
    (
        "built labels left out of meta",
        BAR,
        '{"yticktext": categories} if labels',
        "{} if labels",
        [LEGEND],
    ),
    (
        "checks titled blind to the project",
        "src/ctx/gather/gdc_run.py",
        'if several else word("checks.title")',
        'if False else word("checks.title")',
        ["tests/test_gdc_projects.py"],
    ),
    (
        "M7 variant text not used",
        "src/ctx/proteo/page.py",
        '"variants.tumour_normal" if paired else None',
        "None",
        [f"{TN}::test_page_text_is_written_for_tumour_against_normal"],
    ),
    (
        "normals' figure drawn from tumour rho",
        "src/ctx/proteo/page.py",
        'rho="normal_rho", rho_padj="normal_rho_padj"',
        '"rho", "rho_padj"',
        [f"{TN}::test_normals_draw_their_own_class_figure"],
    ),
    (
        "normals' class medians over tumour-tested genes",
        "src/ctx/proteo/run.py",
        'tested = measured.genes.filter(pl.col("normal_tested"))\n    short = "normal_rho"',
        'tested = measured.genes.filter(pl.col("tested"))\n    short = "normal_rho"',
        [f"{TN}::test_normal_class_medians_read_the_genes_tested_in_normals"],
    ),
    (
        "variants left out of the prose lint",
        "src/ctx/report/prose.py",
        '        yield from check_page(f"{name}: variants.{variant}", sections)\n',
        "        pass\n",
        ["tests/test_prose.py::test_a_variants_sections_keep_the_budgets"],
    ),
    (
        "variants left out of the glossary sections",
        "src/ctx/report/glossary_lint.py",
        'if page and head == "variants" and len(path) > _VARIANT_DEPTH:',
        "if False:",
        ["tests/test_glossary.py::test_a_variants_sections_render_marks_as_their_own_section"],
    ),
    (
        "variant grouped with the page's own section",
        "src/ctx/report/glossary_lint.py",
        'return None if inner is None else f"{path[1]}: {inner}"',
        "return inner",
        ["tests/test_glossary.py::test_a_variants_sections_render_marks_as_their_own_section"],
    ),
]


def run(tests: list[str]) -> str:
    """The tests' summary line."""
    args = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-x", *tests]
    done = subprocess.run(args, cwd=WT, capture_output=True, text=True, check=False)
    return done.stdout.strip().splitlines()[-1] if done.stdout.strip() else done.stderr[-300:]


for name, file, old, new, tests in MUTATIONS:
    path = WT / file
    source = path.read_text()
    assert source.count(old) == 1, (name, source.count(old))
    path.write_text(source.replace(old, new))
    try:
        result = run(tests)
    finally:
        path.write_text(source)
    caught = "failed" in result or "error" in result
    print(f"{'CAUGHT' if caught else 'SURVIVED'}  {name}: {result}", flush=True)
