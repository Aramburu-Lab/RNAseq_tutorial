"""Package M9: plant each mutation of a set, run its tests, restore the file; print the verdicts.

    bash M9_py.sh M9_mutate.py <set>

Each mutation replaces exactly one occurrence. Bytecode is never written and every
``__pycache__`` is removed after the set, so a same-size mutation cannot outlive its restore.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[3]
PYTEST = WT / "docs/plans/logs/M9_pytest.sh"
EXONS = "tests/test_genecard_exons.py"
SETS: dict[str, list[tuple[str, str, str, list[str], str]]] = {
    "h3": [
        (
            "src/ctx/genecard/gtex.py",
            'pl.col("median") >= floor',
            'pl.col("median") > floor',
            [EXONS],
            "the 1 TPM floor taken as strictly above",
        ),
        (
            "src/ctx/genecard/annotation.py",
            'middle = held.group_by("exon")',
            'middle = shares.group_by("exon")',
            [EXONS],
            "reference over every tissue, restricted after",
        ),
        (
            "src/ctx/genecard/annotation.py",
            "        return math.inf\n",
            "        return math.nan\n",
            [EXONS],
            "share over a reference of 0 left blank",
        ),
        (
            "src/ctx/genecard/annotation.py",
            "return -math.inf if reference > 0 else math.nan",
            "return math.nan",
            [EXONS],
            "no share against a reference above 0 left blank",
        ),
        (
            "src/ctx/genecard/annotation.py",
            "return -math.inf if reference > 0 else math.nan",
            "return -math.inf",
            [EXONS],
            "0 over 0 given the bottom step",
        ),
        (
            "src/ctx/genecard/annotation.py",
            "if defined else None",
            'if defined else float("nan")',
            [EXONS],
            "an undefined CV as NaN",
        ),
        (
            "src/ctx/genecard/annotation.py",
            "values.size >= 2  #",
            "values.size >= 1  #",
            [EXONS],
            "a CV over one tissue",
        ),
        (
            "src/ctx/genecard/annotation.py",
            "if x is None or y is None or",
            "if",
            [EXONS],
            "the qPCR pair taking an exon with no CV",
        ),
        (
            "src/ctx/genecard/figures.py",
            "if len(use.tissues) < 2:",
            "if len(use.tissues) < 1:",
            [EXONS],
            "a one-tissue panel",
        ),
        (
            "src/ctx/genecard/figures.py",
            "row_stats=finite,",
            'row_stats={e: v if v is not None else float("nan") for e, v in use.cv.items()},',
            [EXONS],
            '"nan" printed beside a row',
        ),
    ],
}


def plant(path: Path, old: str, new: str) -> str:
    """Replace the one occurrence of ``old``; return the original text."""
    text = path.read_text()
    if text.count(old) != 1:
        raise SystemExit(f"{path.name}: {text.count(old)} matches for {old!r}")
    path.write_text(text.replace(old, new))
    return text


def main(name: str) -> int:
    """Run one set; the exit code counts the survivors."""
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    survived = 0
    for file, old, new, tests, label in SETS[name]:
        path = WT / file
        original = plant(path, old, new)
        try:
            done = subprocess.run(["bash", str(PYTEST), "-q", "-x", "-p", "no:cacheprovider",
                                   *tests], cwd=WT, env=env, capture_output=True,
                                  check=False)  # fmt: skip
        finally:
            path.write_text(original)
        caught = done.returncode != 0
        survived += not caught
        print(f"{'caught' if caught else 'SURVIVED'}: {label}", flush=True)
    for cache in (WT / "src").rglob("__pycache__"):
        shutil.rmtree(cache)
    print(f"{len(SETS[name]) - survived} of {len(SETS[name])} caught")
    return survived


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
