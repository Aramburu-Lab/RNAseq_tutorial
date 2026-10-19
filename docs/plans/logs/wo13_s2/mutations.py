# ruff: noqa: E501 - each planted string is a source line verbatim, so it cannot wrap
"""wo13-s2: plant one mutation per new rule, run the tests that should catch it, restore.

Run from the worktree with main's environment's python and PYTHONPATH at this worktree's src:
    python docs/plans/logs/wo13_s2/mutations.py
Each mutation must make its tests fail; the log is mutations.log beside this file. Bytecode is
never written and __pycache__ under src is cleared after each restore (lessons.md).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
CHECK = "src/ctx/gather/gdc_check.py"
HARMONISE = "src/ctx/clinical/harmonise.py"
PAGE = "src/ctx/clinical/page.py"
COUNTS = "src/ctx/gather/counts.py"
RUN = "src/ctx/gather/counts_run.py"
GDC_TESTS = ("tests/test_dataset_template_gdc.py", "tests/test_clinical_page.py")
PAGE_TESTS = ("tests/test_clinical_page.py",)
HISTORY_TESTS = ("tests/test_gather_counts_history.py", "tests/test_guide_counts.py")
# fmt: off
MUTATIONS = [
    # R1.4: --check names a missing or empty [clinical.published.table3_use]
    ("check never reads table3_use", CHECK,
     "    problems = _missing(spec) + _uses(spec)\n", "    problems = _missing(spec)\n",
     GDC_TESTS, "names_each_edit"),
    ("one blank endpoint passes", HARMONISE,
     "    blank = tuple(e for e in ENDPOINTS if not marks.get(e))\n",
     "    blank = ENDPOINTS if not any(marks.get(e) for e in ENDPOINTS) else ()\n",
     GDC_TESTS, "table3_use_blank or endpoint_note"),
    ("a mark the page refuses passes the check", CHECK,
     '    return [Problem(f"[clinical.published.{TABLE3_USE}]: {fault}")] if fault else []\n',
     "    return []\n", GDC_TESTS, "table3_use_blank"),
    ("lines printed for the template's type", CHECK,
     '        lines = table(f"[clinical.published.{TABLE3_USE}]").replace("XXXX", kind)\n',
     '        lines = table(f"[clinical.published.{TABLE3_USE}]")\n', GDC_TESTS, "table3_use_blank"),
    ("a blank copy does not stop the run", CHECK,
     '                        "the dataset file holds no copy", lines)]  # fmt: skip\n',
     '                        "the dataset file holds no copy", lines, blocks=False)]  # fmt: skip\n',
     GDC_TESTS, "table3_use_blank"),
    ("a caution mark without its reason passes", HARMONISE,
     '    if caution and not marks.get("caution"):\n', "    if False:\n",
     GDC_TESTS, "table3_use_blank or endpoint_note"),
    ("a file without [clinical] asked for marks", CHECK,
     "    if spec.clinical is None or not spec.gdc.cdr_file:\n", "    if not spec.gdc.cdr_file:\n",
     GDC_TESTS, "table3_use_blank"),
    # R1.5: a variable with nothing missing says so
    ("nothing missing printed as not estimable", PAGE,
     '    return "; ".join(f"{r} {n:,}" for r, n in counted) or word("patients.missing_nothing")\n',
     '    return "; ".join(f"{r} {n:,}" for r, n in counted) or word("none")\n',
     PAGE_TESTS, "nothing_missing"),
    ("a gap with no reason passes", PAGE,
     "    if sum(n for _, n in counted) != patients[column].null_count():\n", "    if False:\n",
     PAGE_TESTS, "nothing_missing"),
    # R2.5: Ensembl IDs the annotation lacks, through Ensembl's ID history
    ("history never consulted", COUNTS,
     "    placed = _place(raw.ids, targets, entrez is not None, history)\n",
     "    placed = _place(raw.ids, targets, entrez is not None)\n", HISTORY_TESTS, "history"),
    ("walk stops after one event", COUNTS,
     "                frontier.append(after)\n", "                pass\n", HISTORY_TESTS, "only_one_current"),
    ("walk stops at a current gene", COUNTS,
     "                frontier.append(after)\n",
     "                frontier.extend(() if after in held else (after,))\n", HISTORY_TESTS, "only_one_current"),
    ("version steps read as successors", COUNTS,
     "                          & (old != new))  # fmt: skip\n", "                          )  # fmt: skip\n",
     HISTORY_TESTS, "only_one_current"),
    ("transcript events read", COUNTS,
     '    events = frame.filter((pl.col("type") == "gene") & old.is_not_null() & new.is_not_null()\n',
     "    events = frame.filter(old.is_not_null() & new.is_not_null()\n", HISTORY_TESTS, "only_one_current"),
    ("a history row unseats the matrix's own", COUNTS,
     "        own = [r for r in rows if placed[r][3] not in FALLBACKS]\n",
     "        own = [r for r in rows if placed[r][3] != THROUGH_HGNC]\n", HISTORY_TESTS, "only_one_current"),
    ("a split given to one side", COUNTS,
     '    if len(genes) == 1:\n        return genes, THROUGH_HISTORY, ""\n',
     '    if genes:\n        return genes[:1], THROUGH_HISTORY, ""\n', HISTORY_TESTS, "only_one_current"),
    ("fallback where the own map holds the key", COUNTS,
     "        if key not in targets[kind]:\n", "        if True:\n", HISTORY_TESTS, "only_one_current"),
    ("resource declared, never read", RUN,
     "                          _history(raw, root, spec, releases))  # fmt: skip\n",
     "                          None)  # fmt: skip\n", HISTORY_TESTS, "history_resource or geo_table"),
    ("rows the annotation lacks miscounted", COUNTS,
     '        1 for i in raw.ids if (m := ENSEMBL.fullmatch(i)) and m[1] not in targets["ensembl"]\n',
     "        1 for i in raw.ids if (m := ENSEMBL.fullmatch(i))\n", HISTORY_TESTS, "only_one_current"),
    ("a wider file read by position", COUNTS,
     "    if frame.width != len(HISTORY_COLUMNS):\n", "    if False:\n", HISTORY_TESTS, "another_shape"),
]
# fmt: on


def run(tests: tuple[str, ...], selector: str) -> subprocess.CompletedProcess[str]:
    """The tests that should catch a mutation, by keyword."""
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    return subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider",
                           *tests, "-k", selector], cwd=ROOT, env=env, capture_output=True,
                          text=True, check=False)  # fmt: skip


def main() -> int:
    """Each mutation applied, tested and restored; a surviving one fails the run."""
    log, survived = [], 0
    for name, path, old, new, tests, selector in MUTATIONS:
        file = ROOT / path
        source = file.read_text()
        if source.count(old) != 1:
            log.append(f"NOT APPLIED {name}: {old!r} found {source.count(old)} times")
            survived += 1
            continue
        file.write_text(source.replace(old, new))
        try:
            done = run(tests, selector)
        finally:
            file.write_text(source)
            for cache in ROOT.joinpath("src").rglob("__pycache__"):
                shutil.rmtree(cache)
        caught = done.returncode != 0 and "no tests ran" not in done.stdout
        survived += not caught
        tail = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else done.stderr[-200:]
        log.append(f"{'CAUGHT' if caught else 'SURVIVED'} {name}: {tail}")
    log.append(f"{len(MUTATIONS) - survived} of {len(MUTATIONS)} caught")
    (Path(__file__).with_suffix(".log")).write_text("\n".join(log) + "\n")
    print("\n".join(log))
    return 1 if survived else 0


if __name__ == "__main__":
    raise SystemExit(main())
