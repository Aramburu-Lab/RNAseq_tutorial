"""wo13-f1 M7: plant each mutation, run the tests that should catch it, restore the files.

Usage: wo13_f1_proteo_py.sh python docs/plans/logs/wo13_f1/proteo/mutations.py [name ...]
Each mutation is one or more exact one-match replacements; every file is restored even when
the run fails. Modelled on docs/plans/logs/wo13_f1/mutations.py.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[5]
T = "tests/test_figures_m7.py"
EVERY = f"{T}::test_every_figure_has_a_spec_and_redraws_byte_for_byte"
OPTIONS = f"{T}::test_a_run_with_other_options_redraws_with_them"
PAGE, TABLES = "src/ctx/proteo/page.py", "src/ctx/proteo/figure_tables.py"
DRAW, KINDS = "src/ctx/figures/draw_m7.py", "src/ctx/figures/kinds_m7.py"
EARLY = '    writer.table(m.genes, "genes", figures.columns(*m.genes.columns))\n'
MUTATIONS = [
    ("genes table written after the figures", [
        (PAGE, EARLY, ""),
        (PAGE, "    results.append(_sanity(writer, ctx))\n",
         EARLY + "    results.append(_sanity(writer, ctx))\n")], f"{EVERY}[across_root]"),
    ("figure tables not written", [
        (PAGE, "    figure_tables.write(writer, m, ctx.spec, ctx.options, PAIRED_TOP)\n", "")],
     f"{EVERY}[cross_root]"),
    ("heatmap rows one short", [
        (TABLES, "chosen = figures.heatmap_genes(measured.genes, named, top)",
         "chosen = figures.heatmap_genes(measured.genes, named, top - 1)")],
     f"{EVERY}[across_root] {T}::test_drawn_rows_hold_every_heatmap_gene_past_the_cap"),
    ("annotation rows reversed", [
        (TABLES, 'found = [(TABLES["paired_annotation"], extra["annotation"])]',
         'found = [(TABLES["paired_annotation"], extra["annotation"].reverse())]')],
     f"{EVERY}[across_root]"),
    ("abundance reversed", [
        (TABLES, '"abundance": np.asarray(measured.abundance, dtype=float)',
         '"abundance": np.asarray(measured.abundance, dtype=float)[::-1]')],
     f"{EVERY}[cross_root]"),
    ("presence read inverted", [
        (TABLES, "== 1, 0.0, np.nan)", "== 0, 0.0, np.nan)")], f"{EVERY}[cross_root]"),
    ("method from the default", [
        (TABLES, 'return Options(dataset_file, row["method"],',
         'return Options(dataset_file, "spearman",')], OPTIONS),
    ("min_quantified from the default", [
        (TABLES, 'row["min_quantified"], row["lfc_min"],', '0.7, row["lfc_min"],')], OPTIONS),
    ("gene class from the default", [
        (TABLES, '                   row["shown_class"])', '                   "all")')],
     OPTIONS),
    ("contrast rho read from its low end", [
        (TABLES, '(row["fold_rho"], row["fold_rho_low"]',
         '(row["fold_rho_low"], row["fold_rho_low"]')],
     f"{EVERY}[across_root]"),
    ("cross-study n swapped", [
        (DRAW, '"n_tested": measured.groups.count(declared.tested)',
         '"n_tested": measured.groups.count(declared.reference)')],
     f"{T}::test_a_cross_study_n_counts_each_level_of_the_contrast"),
    ("worked samples not an input of the paired design", [
        (KINDS, '        names.append("worked_samples")\n', "        pass\n")], f"{EVERY}[made]"),
    ("worked gene not checked", [
        (KINDS, "_correlation, _gene_check),", "_correlation),")],
     f"{T}::test_a_worked_gene_the_dataset_file_does_not_name_is_refused"),
    ("class test not run again", [
        (DRAW, "    m7.correlation_metrics(measured, spec.proteo, ctx.options.method)\n",
         '    measured.extra["class_test"] = m7.measure.class_test({"x": __import__("numpy")'
         '.array([0.0, 1.0]), "y": __import__("numpy").array([0.5, 2.0])})\n')],
     f"{EVERY}[across_root]"),
]  # fmt: skip


def run(names: list[str]) -> int:
    """Plant, test, restore; print one line per mutation and return the number survived."""
    survived = 0
    for name, edits, tests in MUTATIONS:
        if names and name not in names:
            continue
        originals = {path: (W / path).read_text() for path, _, _ in edits}
        planted = dict(originals)
        for path, old, new in edits:
            assert planted[path].count(old) == 1, (name, old, planted[path].count(old))
            planted[path] = planted[path].replace(old, new)
        try:
            for path, body in planted.items():
                (W / path).write_text(body)
            done = subprocess.run(["pytest", "-q", "-x", "-n", "2", "-p", "no:cacheprovider",
                                   *tests.split()],
                                  cwd=W, capture_output=True, text=True, check=False)  # fmt: skip
        finally:
            for path, body in originals.items():
                (W / path).write_text(body)
        caught = done.returncode != 0
        survived += not caught
        lines = done.stdout.strip().splitlines()
        last = lines[-1] if lines else ""
        why = next((x.strip() for x in lines if x.startswith("E ")), "")[:160]
        print(f"{'CAUGHT ' if caught else 'SURVIVED'} {name}: {last} | {why}", flush=True)
    return survived


if __name__ == "__main__":
    sys.exit(run(sys.argv[1:]))
