# ruff: noqa: E501 - each mutation is one line of the file it mutates, quoted whole
"""wo13-r1: plant one mutation per new rule, run the tests that should catch it, restore.

Each mutation must make its tests fail; the file is restored byte for byte after each run, and
no bytecode is written (lessons.md, 2026-09-30). Prints one line per mutation and exits 1 if any
survived.
"""

import os
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
PYTHON = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin/python"
S, C, G, P = (
    "scripts/bash/fetch_gdc.sbatch",
    "src/ctx/gather/gdc_check.py",
    "src/ctx/report/guide_gdc.py",
    "src/ctx/pipeline.py",
)
F, T, U = "tests/test_fetch_gdc.py", "tests/test_dataset_template_gdc.py", "tests/test_guide_gdc.py"
MUTATIONS = [
    (
        "md5 on arrival dropped",
        S,
        'if [[ "$(md5of "$dest.part")" == "$3" ]]; then',
        "if true; then",
        F,
    ),
    ("another query's folder taken", S, 'if [[ -s "$saved" ]] && ! cmp', "if false && ! cmp", F),
    (
        "missing files not counted",
        S,
        'must("data files md5 OK, bad, missing", mismatched == gone == 0',
        'must("data files md5 OK, bad, missing", mismatched == 0',
        F,
    ),
    (
        "size run fetches data",
        S,
        'if [[ "$DATA" == 0 ]]; then\n        echo "PASTE bytes',
        'if false; then\n        echo "PASTE bytes',
        F,
    ),
    ("CDR left unlogged", S, '&& log "$cdr_name" "GET $GDC_DATA/$CDR_ID"', "", F),
    (
        "read groups not fetched",
        S,
        "        post gdc_read_groups.json files",
        "        : post gdc_read_groups.json files",
        F,
    ),
    (
        "home check dropped",
        S,
        '"$home"/*) [[ "${SITE_HOME_ALLOWED:-0}" == 1 ]] ||',
        '"$home"/*) true ||',
        F,
    ),
    ("partial record passed", S, "hits == total and not warnings", "not warnings", F),
    (
        "release not compared",
        C,
        "    if found == spec.release:\n        return []",
        "    if True:\n        return []",
        T,
    ),
    (
        "Table 3 zeros pass",
        C,
        "    if all(n == ZEROS for n in given.values()):",
        "    if False:",
        T,
    ),
    (
        "Table 3 mismatch silent",
        C,
        "differ = [e for e in counted if given[e] != counted[e]]",
        "differ = []",
        T,
    ),
    (
        "CDR type not guessed",
        C,
        "                        f'cdr_type = \"{guess}\"')",
        "                        '')",
        T,
    ),
    ("miRNA n not checked", C, 'problems += _n(small, spec, [mirna], "[ncrna.mirna]")', "pass", T),
    (
        "rules time not checked",
        C,
        "    elif not REGISTERED.fullmatch(registered):",
        "    elif False:",
        T,
    ),
    (
        "misspelt gene passed",
        C,
        '"a gene symbol of the counts files")  # fmt: skip',
        '"a gene symbol of the counts files")[:0]  # fmt: skip',
        T,
    ),
    (
        "misspelt miRNA passed",
        C,
        "    known = set(targets.read_names(",
        "    known = {*spec.ncrna.mirnas} | set(targets.read_names(",
        T,
    ),
    ("slug left a placeholder", G, '.replace("<slug>", slug)', "", U),
    ("copy overwrites", G, 'f"cp -n {', 'f"cp {', U),
    ("isoform folder not named", G, 'KINDS = ("mirna", "isoform")', 'KINDS = ("mirna",)', U),
    (
        "GDC check before folder refusals",
        P,
        "and not any(p.blocks for p in problems)",
        "",
        "tests/test_pipeline.py",
    ),
]


def main() -> int:
    """Run every mutation; return 1 if one survived."""
    env = {**os.environ, "PYTHONPATH": str(WT / "src"), "PYTHONDONTWRITEBYTECODE": "1"}
    survived = 0
    for name, path, old, new, tests in MUTATIONS:
        file = WT / path
        original = file.read_bytes()
        text = original.decode()
        if text.count(old) != 1:
            print(f"NOT PLANTED {name}: {old!r} found {text.count(old)} times")
            survived += 1
            continue
        file.write_text(text.replace(old, new))
        try:
            done = subprocess.run([PYTHON, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider",
                                   tests], cwd=WT, env=env, capture_output=True, text=True,
                                  check=False)  # fmt: skip
        finally:
            file.write_bytes(original)
        caught = done.returncode != 0
        survived += not caught
        last = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else ""
        print(f"{'caught' if caught else 'SURVIVED'} {name}: {tests}: {last}", flush=True)
    print(f"{len(MUTATIONS) - survived} of {len(MUTATIONS)} caught")
    return 1 if survived else 0


if __name__ == "__main__":
    sys.exit(main())
