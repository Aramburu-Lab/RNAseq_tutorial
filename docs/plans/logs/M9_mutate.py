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
    "e1": [
        (
            "src/ctx/genecard/figures.py",
            'spread.filter(pl.col("sex").is_in(donors))',
            'spread.filter(pl.col("sex").is_in(SEXES))',
            ["tests/test_genecard_views.py"],
            "F2: a sex view drawing both sexes",
        ),
        (
            "src/ctx/genecard/figures.py",
            "own = every if donors == SEXES else",
            "own = every if True else",
            ["tests/test_genecard_views.py"],
            "F2: a view's lead from the all-donor medians",
        ),
        (
            "src/ctx/genecard/figures.py",
            "for d in DONOR_VIEWS if set(d) & present]",
            "for d in DONOR_VIEWS]",
            ["tests/test_genecard_views.py"],
            "F2: a view of a sex with no donor",
        ),
        (
            "src/ctx/genecard/figures.py",
            '.sort("slot", "sex").select(',
            '.sort("value", "sex").select(',
            ["tests/test_genecard_views.py"],
            "F2: a view in its own order",
        ),
        (
            "src/ctx/genecard/figures.py",
            "DONOR_VIEWS = (SEXES, (gtex.ALL,),",
            "DONOR_VIEWS = ((gtex.ALL,), SEXES,",
            ["tests/test_genecard_views.py"],
            "F2: the built figure not first",
        ),
        (
            "src/ctx/genecard/summary.py",
            "if log else raw",
            "if True else raw",
            ["tests/test_genecard_twins.py"],
            "F3: the scale ignored in Python",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "return log === false ? v : Math.log2(v + 1);",
            "return Math.log2(v + 1);",
            ["tests/test_genecard_twins.py"],
            "F3 twin: the scale ignored",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "if (top === null || r[2] > top[2]) top = r;",
            "if (top === null || r[2] >= top[2]) top = r;",
            ["tests/test_genecard_twins.py"],
            "F3 twin: the last of tied tops",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "return lineages.indexOf(r[1]) >= 0; }), top",
            "return lineages.indexOf(r[1]) >= -1; }), top",
            ["tests/test_genecard_twins.py"],
            "F3 twin: lineages not filtered",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "Math.max.apply(null, values) < ctxM9.DETECTED",
            "Math.max.apply(null, values) < 0",
            ["tests/test_genecard_twins.py"],
            "F3 twin: tau below the detection limit",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "var held = rows.filter(function (r) { return tissues.indexOf(r[1]) >= 0; });",
            "var held = rows;",
            ["tests/test_genecard_twins.py"],
            "F4 twin: tissues not restricted",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "if (!(mean > 0) || n < 2) return null;",
            "if (!(mean > 0) || n < 1) return null;",
            ["tests/test_genecard_twins.py"],
            "F4 twin: a CV over one tissue",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "    if (share > 0) return Infinity;\n",
            "    if (share > 0) return NaN;\n",
            ["tests/test_genecard_twins.py"],
            "F4 twin: the top step lost",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "while (i < bounds.length && bounds[i] <= step) i++;",
            "while (i < bounds.length && bounds[i] < step) i++;",
            ["tests/test_genecard_twins.py"],
            "F4 twin: bounds left-closed",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "if (best === null || worst < best.cv)",
            "if (best === null || worst <= best.cv)",
            ["tests/test_genecard_twins.py"],
            "F4 twin: the last of tied pairs",
        ),
        (
            "src/ctx/report/templates/_m9.js.j2",
            "Object.keys(cv).sort(function (a, b) { return Number(a.slice(1)) - Number(b.slice(1)); });",  # noqa: E501
            "Object.keys(cv).sort();",
            ["tests/test_genecard_twins.py"],
            "F4 twin: exons ordered as text",
        ),
        (
            "src/ctx/genecard/figures.py",
            "levels = [s if c else (none if p is not None else None)",
            "levels = [s if True else (none if p is not None else None)",
            ["tests/test_genecard_evidence.py"],
            "evidence: an uncalled cell in its step",
        ),
        (
            "src/ctx/genecard/figures.py",
            "(none if p is not None else None)",
            "(none if True else None)",
            ["tests/test_genecard_evidence.py"],
            "evidence: an untested cell grey",
        ),
        (
            "src/ctx/genecard/figures.py",
            "EVIDENCE_ROWS = 100",
            "EVIDENCE_ROWS = 1000",
            ["tests/test_genecard_evidence.py"],
            "evidence: no row limit",
        ),
        (
            "src/ctx/genecard/figures.py",
            "return heat, numbers, rest",
            "return heat, numbers, []",
            ["tests/test_genecard_evidence.py"],
            "evidence: the rest lost",
        ),
        (
            "src/ctx/genecard/page.py",
            "if len(ctx.spec.contrasts) < 2:",
            "if len(ctx.spec.contrasts) < 1:",
            ["tests/test_genecard_evidence.py"],
            "evidence: drawn for one contrast",
        ),
        (
            "src/ctx/genecard/page.py",
            "results: list[Result] = [] if evidence is None else [evidence.group]",
            "results: list[Result] = []",
            ["tests/test_genecard_evidence.py"],
            "evidence: left off the page",
        ),
        (
            "src/ctx/genecard/page.py",
            'c.symbol if seen[c.symbol] == 1 else f"{c.symbol} {c.gene_id}"',
            "c.symbol",
            ["tests/test_genecard_evidence.py"],
            "evidence: two genes of one symbol on one row",
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
