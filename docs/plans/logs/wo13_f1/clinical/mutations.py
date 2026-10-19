"""wo13-f1 M10: plant each mutation, run the tests that should catch it, restore the file.

Usage: wo13_f1_clinical_py.sh python docs/plans/logs/wo13_f1/clinical/mutations.py [name ...]
Each mutation is one exact replacement; the file is restored even when the run fails.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[5]
T = "tests/test_figures_m10.py"
EVERY = f"{T}::test_every_figure_redraws_byte_for_byte"
OPTIONS = f"{T}::test_the_runs_options_are_the_ones_redrawn"
DRAW = "src/ctx/figures/draw_m10.py"
TABLES = "src/ctx/clinical/figure_tables.py"
KINDS = "src/ctx/figures/kinds_m10.py"
REFUSED = f"{T}::test_an_endpoint_or_box_the_file_does_not_draw_is_refused"
# The two tables written once every figure is drawn: their specs find them only afterwards.
PAGE = (W / "src/ctx/clinical/page.py").read_text()
WRITE = "    figure_tables.write(writer, ctx)\n"
START = PAGE.index(WRITE)
SPAN = PAGE[START : PAGE.index("    ]  # fmt: skip\n", START) + len("    ]  # fmt: skip\n")]
AFTER = (SPAN, SPAN.replace(WRITE, "") + WRITE)
MUTATIONS = [
    ("panel values not written", TABLES,
     '    writer.table(frame, "panel_values", {**figures.columns("columns.panel_values"), '
     "**samples})\n", "", EVERY),
    ("build genes not marked", TABLES, '"build": [r in build for r in rows]',
     '"build": [r == rows[0] for r in rows]', EVERY),
    ("values matched by position", TABLES,
     "if (found := frame.columns[len(GENE_COLUMNS) :]) != samples:",
     "if len(found := frame.columns[len(GENE_COLUMNS) :]) != len(samples):",
     f"{T}::test_values_of_other_samples_are_refused"),
    ("values read as floats", TABLES, "values = frame.select(samples).to_numpy().astype(np.int16)",
     "values = frame.select(samples).to_numpy() // 2", EVERY),
    ("live block of the build genes alone", DRAW,
     "carried = panel.Panel(np.arange(genes.height), panel.encode(values), genes.height,",
     "carried = panel.Panel(np.arange(genes.height), panel.encode(values[chosen]), genes.height,",
     EVERY),
    ("screen left out of the live panel", DRAW, "page.screen_table(ctx, result))",
     "page.screen_table(ctx, result).head(0))", EVERY),
    ("split taken as the default", DRAW, 'split=row["split"]', 'split="median"', OPTIONS),
    ("horizon taken as the default", DRAW, 'horizon=row["horizon"]', "horizon=60.0", OPTIONS),
    ("covariates taken as the default", DRAW, 'covariates=row["covariates"]',
     'covariates=("age", "stage")', OPTIONS),
    ("strata taken as none", DRAW, 'strata=row["strata"])', "strata=())", OPTIONS),
    ("endpoint of the dataset file", DRAW, 'spec=spec, options=options, endpoint=row["endpoint"],',
     "spec=spec, options=options, endpoint=spec.clinical.endpoint,", OPTIONS),
    ("source of the dataset file", DRAW, 'source=row["source"]',
     'source=text.word("clinical", "live.build")', OPTIONS),
    ("endpoint not checked", KINDS,
     'lambda p: f"er_{p.endpoint.lower()}", _read, _endpoint_check),',
     'lambda p: f"er_{p.endpoint.lower()}", _read),', REFUSED),
    ("box not checked", KINDS, 'lambda p: f"expression_{p.box}", _scored, _box_check),',
     'lambda p: f"expression_{p.box}", _scored),', REFUSED),
    ("cramers v in hashed order", "src/ctx/clinical/analyses.py",
     "both.group_by(a, b, maintain_order=True)", "both.group_by(a, b)",
     f"{T}::test_cramers_v_is_the_same_every_call"),
    ("tables written after the figures", "src/ctx/clinical/page.py", AFTER[0], AFTER[1], EVERY),
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
        planted = original.replace(old, new)
        if "text.word" in new and "from ctx.report import text" not in planted:
            head = "from ctx.report.builder import"
            planted = planted.replace(head, f"from ctx.report import text\n{head}", 1)
        file.write_text(planted)
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
