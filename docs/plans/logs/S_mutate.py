# ruff: noqa: T201, S603, S607, E501 - a proof record: one planted mutation per new rule, run once
"""Package S: plant each mutation, run the test that should catch it, restore, report.

Run by S_py.sh from the worktree; pytest runs with PYTHONDONTWRITEBYTECODE=1 and every
__pycache__ under src is removed after each restore (lessons.md, mutation harness).
"""

import os
import shutil
import subprocess
from pathlib import Path

WT = Path(__file__).resolve().parents[3]
T_TN, T_GP = "tests/test_proteo_tumour_normal.py", "tests/test_gdc_projects.py"
MUTATIONS = [
    ("exclusion first", "src/ctx/gather/gdc.py",
     'pl.when(pl.col("case").is_in(list(excluded))).then(pl.lit("excluded"))',
     'pl.when(pl.lit(value=False)).then(pl.lit("excluded"))', T_GP + "::test_an_excluded_case_drops_first_and_m5_never_sees_it"),
    ("exclude_reason required", "src/ctx/dataset.py",
     "if cases and (not isinstance(reason, str) or not reason):", "if False:",
     T_GP + "::test_an_excluded_case_drops_first_and_m5_never_sees_it"),
    ("accessions renumbered", "src/ctx/gather/gdc_run.py",
     "sheet = pl.concat([r.sheet for r in found]).with_columns(run_accession=gdc.run_accession())",
     "sheet = pl.concat([r.sheet for r in found])", T_GP + "::test_two_projects_make_one_sheet_with_a_project_column"),
    ("file's project order", "src/ctx/gather/gdc_run.py",
     "    for project in spec.gdc.projects:\n        one = spec.gdc.one(project)",
     "    for project in by:\n        one = spec.gdc.one(project)", T_GP + "::test_two_projects_make_one_sheet_with_a_project_column"),
    ("repeatable --gdc-dir", "src/ctx/cli.py", "gdc_dirs=gdc_dir,", "gdc_dirs=gdc_dir[-1:],",
     T_GP + "::test_two_project_gather_writes_the_sheet_m5_fits_on"),
    ("projects schema", "src/ctx/dataset.py", "not len(set(listed)) == len(listed) > 1",
     "not len(set(listed)) == len(listed) > 0", T_GP + "::test_schema_refuses_project_and_projects_or_one_listed"),
    ("biospecimens logged", "src/ctx/gather/gdc.py", "    if logged.get(name) != spec.release:\n",
     "    if False:\n", T_GP + "::test_biospecimens_read_and_checked"),
    ("missing not zero", "src/ctx/proteo/measure.py", "    differences = tested - reference\n",
     "    differences = np.nan_to_num(tested) - np.nan_to_num(reference)\n",
     T_TN + "::test_paired_calls_leave_a_missing_value_out_never_zero"),
    ("share over both groups", "src/ctx/proteo/folds.py", "np.hstack([tested, reference]), (lfc, pvalue)",
     "tested, (lfc, pvalue)", T_TN + "::test_paired_calls_leave_a_missing_value_out_never_zero"),
    ("positive rho only", "src/ctx/proteo/measure.py", "(rho > 0) & (padj < fdr)", "(padj < fdr)",
     T_TN + "::test_positive_share_counts_defined_genes_only"),
    ("below, not above", "src/ctx/proteo/run.py",
     'below, r = measure.shifted(other, inside["difference"].to_numpy())',
     'below, r = measure.shifted(inside["difference"].to_numpy(), other)',
     T_TN + "::test_planted_shifts_found_on_the_paired_cases"),
    ("excluded case's columns", "src/ctx/proteo/run.py",
     'chosen = dict(gdc.chosen(sheet, sample_type).select("case", "sample").rows())',
     'chosen = dict(sheet.filter(pl.col("sample_type") == sample_type)'
     '.select("case", "sample").rows())', T_TN + "::test_both_correlations_on_their_own_samples"),
    ("normals' own correlation", "src/ctx/proteo/run.py", 'normal = frames["normal"].drop("gene")',
     'normal = frames["tumour"].drop("gene")', T_TN + "::test_both_correlations_on_their_own_samples"),
    ("paired t, not Welch", "src/ctx/proteo/folds.py", "    lfc, pvalue = measure.paired(tested, reference)\n",
     "    lfc, pvalue = measure.welch(np.hstack([tested, reference]),"
     " np.arange(2 * tested.shape[1]) < tested.shape[1])\n",
     T_TN + "::test_paired_calls_leave_a_missing_value_out_never_zero"),
    ("proteome rows only", "src/ctx/proteo/tables.py",
     "if line[at] is None or not re.search(rows, str(line[at])):", "if line[at] is None:",
     T_TN + "::test_plex_map_keeps_the_proteome_rows_only"),
    ("TSV header checked", "src/ctx/proteo/tables.py",
     "if not rows or tuple(rows[0][: len(descriptors)]) != descriptors:", "if not rows:",
     T_TN + "::test_tsv_reader_refuses_another_header_and_reads_na"),
    ("Other Metadata asked", "src/ctx/proteo/pdc.py",
     'CATEGORIES = ("Publication Supplementary Material", "Other Metadata")',
     'CATEGORIES = ("Publication Supplementary Material",)', "tests/test_pdc.py::test_other_metadata_resolved"),
    ("one listing over all", "src/ctx/proteo/pdc.py", "    if len(files) != 1 or not",
     "    files = files[:1]\n    if len(files) != 1 or not", "tests/test_pdc.py::test_listed_in_two_categories_refused"),
    ("rules as registered", "conf/datasets/cptac_ccrcc.toml",
     'conditions = [{ metric = "pergene.tumour.positive_share", op = ">", value = 0.5 }]',
     'conditions = [{ metric = "pergene.tumour.positive_share", op = ">", value = 0.4 }]',
     T_TN + "::test_ccrcc_rules_in_force_equal_the_registered_ones"),
]  # fmt: skip


def clean() -> None:
    """Every bytecode cache under src, so no mutated .pyc outlives its source."""
    for cache in (WT / "src").rglob("__pycache__"):
        shutil.rmtree(cache)


caught = 0
for name, path, old, new, test in MUTATIONS:
    file = WT / path
    good = file.read_text()
    assert good.count(old) == 1, (name, "not found once")
    file.write_text(good.replace(old, new))
    try:
        done = subprocess.run(
            [
                "bash",
                str(WT / "docs/plans/logs/S_pytest.sh"),
                "-q",
                "-x",
                "-p",
                "no:cacheprovider",
                test,
            ],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        )
    finally:
        file.write_text(good)
        clean()
    failed = done.returncode != 0
    caught += failed
    print(f"{'CAUGHT' if failed else 'SURVIVED'}  {name}  ({test.split('::')[-1]})")
print(f"{caught} of {len(MUTATIONS)} mutations caught")
