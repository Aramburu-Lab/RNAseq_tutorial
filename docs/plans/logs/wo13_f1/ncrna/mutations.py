"""wo13-f1 M8: plant each mutation, run the tests that should catch it, restore the file.

Usage: wo13_f1_ncrna_py.sh python docs/plans/logs/wo13_f1/ncrna/mutations.py [name ...]
Each mutation is one exact replacement; the file is restored even when the run fails.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[5]
T = "tests/test_figures_m8.py"
STEP1 = f"{T}::test_step1_figures_redraw_byte_for_byte"
OPTIONS = f"{T}::test_the_run_s_options_are_the_ones_redrawn"
GDC = f"{T}::test_gdc_figures_redraw_byte_for_byte"
DOOR = f"{T}::test_a_drawer_takes_the_sheet_through_the_gdc_door"
MIRNA_DOOR = f"{T}::test_the_mirna_sheet_passes_the_gdc_door"
MOVED = f"{T}::test_a_moved_target_file_is_refused"
F5 = f"{T}::test_library_types_reads_every_file_of_the_sheet"
DRAW = "src/ctx/figures/draw_m8.py"
MUTATIONS = [
    ("records and options written after the first figure", "src/ctx/ncrna/page.py",
     "    figure_tables.write(writer, ctx)\n    shown, metrics = library(writer, ctx)\n",
     "    shown, metrics = library(writer, ctx)\n    figure_tables.write(writer, ctx)\n", STEP1),
    ("gene records out of the counts' order", "src/ctx/ncrna/figure_tables.py",
     '(ctx.inputs.genes, text.raw(PAGE, "tables.gene_biotypes"))',
     '(ctx.inputs.genes.sort("gene_name", descending=True), '
     'text.raw(PAGE, "tables.gene_biotypes"))', STEP1),
    ("default groups drawn", DRAW, "chosen = figure_tables.biotypes_of(options)",
     "chosen = biotypes.TESTED", OPTIONS),
    ("default count drawn", DRAW,
     'min_count, gtf = int(options["min_count"][0]), Path(options["gtf"][0])',
     'min_count, gtf = MIN_COUNT, Path(options["gtf"][0])', OPTIONS),
    ("groups read back reversed", "src/ctx/ncrna/figure_tables.py",
     'tuple(options["biotypes"][0].split(SEPARATOR))',
     'tuple(reversed(options["biotypes"][0].split(SEPARATOR)))', OPTIONS),
    ("sheet past the door", DRAW,
     'return gdc.chosen(sheet.unique("sample", keep="first", maintain_order=True))',
     'return (sheet.unique("sample", keep="first", maintain_order=True))', DOOR),
    ("miRNA volcano's sheet past the door", DRAW,
     'select_samples(gdc.chosen(_table(inputs, "mirna_samplesheet")), contrast, {})',
     'select_samples(_table(inputs, "mirna_samplesheet"), contrast, {})', MIRNA_DOOR),
    ("pairs' miRNA sheet past the door", DRAW, "gdc.chosen(mirnas) if mirnas.height else mirnas",
     "mirnas", MIRNA_DOOR),
    ("target files unread", DRAW, "targets.read(named) if named else NO_TARGETS", "NO_TARGETS",
     GDC),
    ("target files out of the shift's spec", "src/ctx/figures/kinds_m8.py",
     '_with("gene_biotypes", "results", "mirna_results", targets=True)',
     '_with("gene_biotypes", "results", "mirna_results")', MOVED),
    ("F5 handed the chosen files", DRAW,
     '_table(inputs, "samplesheet"), _table(inputs, "mirna_mature_counts")',
     '_sheet(inputs), _table(inputs, "mirna_mature_counts")', F5),
]  # fmt: skip


def run(names: list[str]) -> int:
    """Plant, test, restore; print one line per mutation and return the number survived."""
    survived = 0
    for name, path, old, new, tests in MUTATIONS:
        if names and name not in names:
            continue
        file = W / path
        original = file.read_text()
        assert original.count(old) == 1, (name, original.count(old))
        file.write_text(original.replace(old, new))
        try:
            done = subprocess.run(["pytest", "-q", "-x", "-p", "no:cacheprovider", *tests.split()],
                                  cwd=W, capture_output=True, text=True, check=False)  # fmt: skip
        finally:
            file.write_text(original)
        caught = done.returncode != 0
        survived += not caught
        lines = done.stdout.strip().splitlines()
        last = lines[-1] if lines else ""
        why = next((x.strip() for x in lines if x.startswith("E ")), "")[:160]
        print(f"{'CAUGHT ' if caught else 'SURVIVED'} {name}: {last} | {why}", flush=True)
    return survived


if __name__ == "__main__":
    sys.exit(run(sys.argv[1:]))
