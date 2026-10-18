# ruff: noqa: E501 - each planted string is a source line verbatim, so it cannot wrap
"""wo13-r2: plant one mutation per new rule, run the tests that should catch it, restore.

Run from the worktree with main's environment's python and PYTHONPATH at this worktree's src:
    python docs/plans/logs/wo13_r2_mutations.py
Each mutation must make its tests fail; the log is wo13_r2_mutations.log beside this file.
Bytecode is never written and __pycache__ is cleared after each restore (lessons.md).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TESTS = ["tests/test_gather_counts.py", "tests/test_guide_counts.py"]
COUNTS = "src/ctx/gather/counts.py"
RUN = "src/ctx/gather/counts_run.py"
MUTATIONS = [
    ("fractional refused", COUNTS, "if fractional.any():", "if False:", "scale or listed"),
    ("negative refused", COUNTS, "(seen < 0).sum()", "(seen < -1e9).sum()", "scale"),
    ("scaled sums refused", COUNTS, "    elif _scaled(sums):", "    elif False:", "scale or edge"),
    (
        "tolerance edge <=",
        COUNTS,
        "<= SCALED_TOLERANCE * PER_MILLION",
        "< SCALED_TOLERANCE * PER_MILLION",
        "edge",
    ),
    ("missing refused", COUNTS, "if (missing := int((~finite).sum())) > 0:", "if False:", "scale"),
    (
        "duplicate gene by stable ID",
        COUNTS,
        "stable = [m[1] if (m := ENSEMBL.fullmatch(i)) else i for i in raw.ids]",
        "stable = list(raw.ids)",
        "listed",
    ),
    (
        "duplicate sample",
        COUNTS,
        '(("gene", stable), ("sample", list(raw.samples)))',
        '(("gene", stable),)',
        "listed",
    ),
    (
        "ensembl by stable part",
        COUNTS,
        '        by["ensembl"][stable[1] if stable else gene_id].add(gene_id)',
        '        by["ensembl"][gene_id].add(gene_id)',
        "end_to_end",
    ),
    (
        "ambiguous left out",
        COUNTS,
        "        elif len(genes) == 1:",
        "        elif len(genes) >= 1:",
        "by_kind",
    ),
    (
        "shared gene left out",
        COUNTS,
        "        if len(rows) > 1:",
        "        if len(rows) > 99:",
        "by_kind",
    ),
    (
        "entrez through GTF",
        COUNTS,
        "found[entrez.strip()].add(gene_of[transcript])",
        "found[entrez.strip()].add(transcript)",
        "through_the_gtf",
    ),
    ("entrez needs table", RUN, "    if n and path is None:", "    if False:", "need_the_table"),
    (
        "version tally",
        COUNTS,
        'tally["equal" if ours is not None and ours[2] == given[2] else "differ"] += 1',
        'tally["equal"] += 1',
        "end_to_end",
    ),
    (
        "featureCounts columns",
        COUNTS,
        "    if tuple(names[1:6]) == FEATURECOUNTS:",
        "    if False:",
        "featurecounts",
    ),
    (
        "summary rows removed",
        COUNTS,
        "summary = ids.str.starts_with(HTSEQ_PREFIX) | ids.is_in(list(STAR_ROWS))",
        "summary = ids.is_in(list(STAR_ROWS))",
        "featurecounts",
    ),
    (
        "csv by ending",
        COUNTS,
        'return "," if name.endswith(".csv") else "\\t"',
        'return "\\t"',
        "csv_by_ending",
    ),
    ("text cell named", COUNTS, "        if bad.any():", "        if False:", "csv_by_ending"),
    ("left-out counts kept", COUNTS, "float(totals[i])", "0.0", "end_to_end or by_kind"),
    (
        "matrix sheet marked",
        RUN,
        "**{counts.SOURCE_COLUMN: pl.lit(matrix.name)}",
        "**{'matrix_name': pl.lit(matrix.name)}",
        "end_to_end or refuse",
    ),
    (
        "process refuses",
        "src/ctx/process/run.py",
        " or refuse_matrix(sheet, path, MODULE)",
        "",
        "refuse",
    ),
    (
        "quant refuses",
        "src/ctx/quant/run.py",
        " or counts.refuse_sheet(sheet, path, MODULE)",
        "",
        "refuse",
    ),
    (
        "qc --reads refuses",
        "src/ctx/qc/run.py",
        "    if options.reads and (",
        "    if False and (",
        "refuse",
    ),
    (
        "ctx run carries --counts",
        "src/ctx/pipeline.py",
        '"gather": base + (gdc_args if gdc else matrix or fetchngs + rnaseq),',
        '"gather": base + (gdc_args if gdc else fetchngs + rnaseq),',
        "plans",
    ),
    (
        "ctx run counts for de",
        "src/ctx/pipeline.py",
        "counts = gdc or options.rnaseq_dir is not None or bool(matrix)",
        "counts = gdc or options.rnaseq_dir is not None",
        "plans",
    ),
    (
        "--check reads the matrix",
        "src/ctx/pipeline.py",
        "        return table.is_file() and options.counts.is_file() and spec.counts is not None",
        "        return False",
        "plans",
    ),
    (
        "guide slug placeholder",
        "src/ctx/report/guide.py",
        "spec.counts is not None and spec.slug != MINE",
        "spec.counts is not None",
        "placeholder",
    ),
    (
        "guide size apart",
        "src/ctx/report/guide.py",
        '            \'cd "$RAW" && curl -fL -o "$FILE" "$URL"\',',
        '            \'curl -fsSIL "$URL" | grep -i content-length; cd "$RAW" && curl -fL -o "$FILE" "$URL"\',',
        "size_before",
    ),
    (
        "geo donor pattern refused",
        "src/ctx/gather/geo_table.py",
        "        if found is None or not found.groups():",
        "        if found is None and False:",
        "pattern_without",
    ),
    (
        "geo key normalised",
        "src/ctx/gather/geo_table.py",
        'return re.sub(r"[^a-z0-9]+", "_", text.strip().lower()).strip("_") or "characteristic"',
        "return text.strip()",
        "writes_samples",
    ),
]


def run(selector: str) -> subprocess.CompletedProcess[str]:
    """The tests that should catch a mutation, by keyword."""
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    return subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider",
                           *TESTS, "-k", selector], cwd=ROOT, env=env, capture_output=True,
                          text=True, check=False)  # fmt: skip


def main() -> int:
    """Each mutation applied, tested and restored; a surviving one fails the run."""
    log, survived = [], 0
    for name, path, old, new, selector in MUTATIONS:
        file = ROOT / path
        source = file.read_text()
        if source.count(old) != 1:
            log.append(f"NOT APPLIED {name}: {old!r} found {source.count(old)} times")
            survived += 1
            continue
        file.write_text(source.replace(old, new))
        try:
            done = run(selector)
        finally:
            file.write_text(source)
            for cache in ROOT.joinpath("src").rglob("__pycache__"):
                shutil.rmtree(cache)
        caught = done.returncode != 0
        survived += not caught
        tail = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else done.stderr[-200:]
        log.append(f"{'CAUGHT' if caught else 'SURVIVED'} {name}: {tail}")
    log.append(f"{len(MUTATIONS) - survived} of {len(MUTATIONS)} caught")
    (Path(__file__).with_suffix(".log")).write_text("\n".join(log) + "\n")
    print("\n".join(log))
    return 1 if survived else 0


if __name__ == "__main__":
    raise SystemExit(main())
