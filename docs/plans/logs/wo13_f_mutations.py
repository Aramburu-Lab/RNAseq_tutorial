"""wo13-f: plant each mutation, run the tests that should catch it, restore the file.

Usage: wo13_f_py.sh python docs/plans/logs/wo13_f_mutations.py <piece> [name ...]
Each mutation is one exact replacement; the file is restored even when the run fails.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[3]
FIG = "tests/test_figures.py"
MUTATIONS = {
    "f1": [
        (
            "fit rows reversed",
            "src/ctx/figures/draw_m5.py",
            "    table = pl.read_parquet(results)\n",
            "    table = pl.read_parquet(results).reverse()\n",
            f"{FIG}::test_every_spec_redraws_its_data_byte_for_byte",
        ),
        (
            "registry copy ignored",
            "src/ctx/figures/draw_m1.py",
            '    with _registry(inputs["qc_registry"]):\n',
            "    with _registry(registry.resolve(None)):\n",
            f"{FIG}::test_gather_figures_redraw_from_the_sheet",
        ),
        (
            "loader stops at the first problem",
            "src/ctx/figures/params.py",
            "    values = _fields(raw, kind, problems)\n",
            "    values = _fields(raw, kind, problems)\n"
            "    if problems:\n        raise ParamsError(problems[0])\n",
            f"{FIG}::test_a_bad_params_file_lists_every_problem",
        ),
        (
            "hash not checked",
            "src/ctx/figures/params.py",
            "        if name in hashes and hashes[name] != digest:\n",
            "        if False:\n",
            f"{FIG}::test_a_bad_params_file_lists_every_problem",
        ),
        (
            "a root not recognised",
            "src/ctx/figures/command.py",
            "    for place in (folder.resolve(), *folder.resolve().parents):\n",
            "    for place in ():\n",
            f"{FIG}::test_figure_refuses_a_folder_inside_a_ctx_root",
        ),
        (
            "setting not passed",
            "src/ctx/figures/draw_m5.py",
            "    data = m5_figures.calls_bars(p.lesson, lesson, fits, sizes, setting)\n",
            "    data = m5_figures.calls_bars(p.lesson, lesson, fits, sizes)\n",
            f"{FIG}::test_a_lesson_redrawn_at_another_setting_counts_as_the_grid",
        ),
        (
            "README loses the ID",
            "src/ctx/viz/save.py",
            '            readme.write("\\n" + record.readme(spec))\n',
            '            readme.write("")\n',
            f"{FIG}::test_every_spec_redraws_its_data_byte_for_byte",
        ),
        (
            "a kind matches nothing",
            "src/ctx/figures/kinds.py",
            'lambda p: f"lfc_{p.x}_vs_{p.y}"',
            'lambda p: f"lfc_{p.x}_{p.y}"',
            f"{FIG}::test_every_spec_redraws_its_data_byte_for_byte",
        ),
        (
            "module number from the folder",
            "src/ctx/figures/spec.py",
            '"ncrna": 8, "genecard": 9, "clinical": 10,',
            '"ncrna": 8, "genecard": 11, "clinical": 12,',
            f"{FIG}::test_ids_round_trip",
        ),
    ],
}


def run(piece: str, names: list[str]) -> int:
    """Plant, test, restore; print one line per mutation and return the number survived."""
    survived = 0
    for name, path, old, new, tests in MUTATIONS[piece]:
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
        last = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else ""
        print(f"{'CAUGHT ' if caught else 'SURVIVED'} {name}: {last}", flush=True)
    return survived


if __name__ == "__main__":
    sys.exit(run(sys.argv[1], sys.argv[2:]))
