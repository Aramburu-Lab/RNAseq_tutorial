"""wo13-f: plant each mutation, run the tests that should catch it, restore the file.

Usage: wo13_f_py.sh python docs/plans/logs/wo13_f_mutations.py <piece> [name ...]
Each mutation is one exact replacement; the file is restored even when the run fails.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[3]
FIG = "tests/test_figures.py"
PANEL = "tests/test_figure_panel.py"
BYTES = "tests/test_figure_bytes.py"
COHORTS = "tests/test_cohorts.py"
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
    "f2": [
        (
            "pasted block stops the shell",
            "src/ctx/figures/panel.py",
            '        *_environment(spec), "",\n',
            '        "set -e", *_environment(spec), "",\n',
            f"{PANEL}::test_pasted_block_never_stops_the_shell_and_the_script_does",
        ),
        (
            "script without pipefail",
            "src/ctx/figures/panel.py",
            'SCRIPT_HEAD = "#!/usr/bin/env bash\\nset -euo pipefail\\n"',
            'SCRIPT_HEAD = "#!/usr/bin/env bash\\n"',
            f"{PANEL}::test_pasted_block_never_stops_the_shell_and_the_script_does",
        ),
        (
            "data hash of another file",
            "src/ctx/figures/record.py",
            "data_sha256=sha256(data_csv)",
            "data_sha256=sha256(paths['dataset_file'])",
            f"{PANEL}::test_the_downloaded_script_redraws_the_data_byte_for_byte",
        ),
        (
            "copy line under the wrong base",
            "src/ctx/figures/panel.py",
            "        return f'cp \"$ROOT/{item.source}\" .'\n",
            "        return f'cp \"$ROOT/outputs/{item.source}\" .'\n",
            f"{PANEL}::test_the_downloaded_script_redraws_the_data_byte_for_byte",
        ),
        (
            "an input left out",
            "src/ctx/figures/panel.py",
            "        *(_copy(i) for i in spec.inputs),\n",
            "        *(_copy(i) for i in spec.inputs[1:]),\n",
            f"{PANEL}::test_the_downloaded_script_redraws_the_data_byte_for_byte",
        ),
        (
            "module page without the panel",
            "src/ctx/viz/save.py",
            "        reproduce=_PANELS.get(figure_stem),\n",
            "        reproduce=None,\n",
            f"{PANEL}::test_every_spec_has_its_panel_on_the_module_page",
        ),
        (
            "commit typed",
            "src/ctx/figures/panel.py",
            'f"git -C {CLONE} checkout {commit}",',
            'f"git -C {CLONE} checkout main",',
            f"{PANEL}::test_environment_follows_where_ctx_came_from",
        ),
    ],
    "f3": [
        (
            "PDF keeps its CreationDate",
            "src/ctx/viz/save.py",
            'NO_DATE = {"pdf": {"CreationDate": None}, "svg": {"Date": None}}',
            'NO_DATE = {"pdf": {}, "svg": {"Date": None}}',
            BYTES,
        ),
        (
            "SVG keeps its dc:date",
            "src/ctx/viz/save.py",
            'NO_DATE = {"pdf": {"CreationDate": None}, "svg": {"Date": None}}',
            'NO_DATE = {"pdf": {"CreationDate": None}, "svg": {}}',
            BYTES,
        ),
        (
            "SVG ids salted at random",
            "src/ctx/viz/save.py",
            '    with mpl.rc_context({"svg.hashsalt": SVG_SALT}):\n',
            "    with mpl.rc_context({}):\n",
            BYTES,
        ),
    ],
    "f4": [
        (
            "bundle loader stops at the first problem",
            "src/ctx/cohorts.py",
            "    slugs = [e.slug for e in entries if e.slug]\n",
            "    if problems:\n        raise CohortError(problems[0])\n"
            "    slugs = [e.slug for e in entries if e.slug]\n",
            f"{COHORTS}::test_a_bad_bundle_lists_every_problem",
        ),
        (
            "a doubled module passes",
            "src/ctx/cohorts.py",
            "    if len(set(value)) != len(value):\n",
            "    if False:\n",
            f"{COHORTS}::test_a_bad_bundle_lists_every_problem",
        ),
        (
            "a zero cap passes",
            "src/ctx/cohorts.py",
            "isinstance(value, int) and not isinstance(value, bool) and value > 0",
            "isinstance(value, int) and not isinstance(value, bool)",
            f"{COHORTS}::test_a_bad_bundle_lists_every_problem",
        ),
        (
            "the listing preferred over the choice",
            "src/ctx/cohorts.py",
            "    if chosen := active():\n",
            "    if (chosen := active()) and False:\n",
            f"{COHORTS}::test_cohort_of_takes_the_choice_then_the_listing_then_own",
        ),
        (
            "the option not passed to the figures",
            "src/ctx/cohorts.py",
            "    os.environ[COHORT_ENV] = name\n",
            "    os.environ.get(COHORT_ENV)\n",
            f"{COHORTS}::test_cohort_option_names_figures_and_refuses_an_unknown_cohort",
        ),
    ],
    "f5": [
        (
            "ctx run gives de another option",
            "src/ctx/pipeline.py",
            '        "de": base if spec.contrasts and spec.de and counts else None,\n',
            '        "de": (*base, "--min-count", "20") if spec.contrasts and spec.de and counts'
            " else None,\n",
            "tests/test_pipeline.py::test_run_matches_the_module_commands_through_de_and_function",
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
        lines = done.stdout.strip().splitlines()
        last = lines[-1] if lines else ""
        why = next((x.strip() for x in lines if x.startswith("E ")), "")[:160]
        print(f"{'CAUGHT ' if caught else 'SURVIVED'} {name}: {last} | {why}", flush=True)
    return survived


if __name__ == "__main__":
    sys.exit(run(sys.argv[1], sys.argv[2:]))
