# ruff: noqa: E501 - each planted string is a source line verbatim, so it cannot wrap
"""wo13-u: plant one mutation per new rule, run the tests that should catch it, restore.

Run from the worktree with main's environment's python and PYTHONPATH at this worktree's src:
    python docs/plans/logs/wo13_u_mutations.py
Each mutation must make its tests fail; the log is wo13_u_mutations.log beside this file.
Bytecode is never written and __pycache__ is cleared after each restore (lessons.md).
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TESTS = ["tests/test_gather_counts_rsem.py", "tests/test_gather_counts.py",
         "tests/test_guide_counts.py", "tests/test_dataset_load.py"]  # fmt: skip
COUNTS = "src/ctx/gather/counts.py"
RSEM = "src/ctx/gather/counts_rsem.py"
RUN = "src/ctx/gather/counts_run.py"
GEO = "src/ctx/gather/geo_table.py"
MUTATIONS = [
    # RSEM expected counts
    (
        "declared fractions pass",
        COUNTS,
        "    if fractional.any() and not expected:",
        "    if fractional.any():",
        "rounded_and_stated or folder_combined",
    ),
    (
        "declared log2 refused",
        COUNTS,
        "not _scaled(sums) and seen.max(initial=0) > LOG_CEILING",
        "not _scaled(sums)",
        "refuses_log",
    ),
    (
        "declared TPM refused",
        COUNTS,
        "kind == RSEM_EXPECTED and not _scaled(sums) and",
        "kind == RSEM_EXPECTED and",
        "refuses_log",
    ),
    (
        "raw kind names the declaration",
        COUNTS,
        '        hint = hint if kind == RAW else ""',
        '        hint = ""',
        "undeclared",
    ),
    (
        "round after the value check",
        RUN,
        "    raw = counts_rsem.read(matrix, spec.counts)\n",
        "    raw = counts_rsem.read(matrix, spec.counts)\n    raw = counts_rsem.rounded(raw)[0] if kind == counts.RSEM_EXPECTED else raw\n",
        "rounds_after",
    ),
    (
        "round, never truncate",
        RSEM,
        "    whole = np.rint(raw.values)",
        "    whole = np.trunc(raw.values)",
        "rounded_and_stated",
    ),
    (
        "halves to even",
        RSEM,
        "    whole = np.rint(raw.values)",
        "    whole = np.floor(raw.values + 0.5)",
        "rounded_and_stated",
    ),
    (
        "rounding in every README",
        RUN,
        '        releases["values"] = text.fill(PAGE, "rounding.readme", **vars(rounding))',
        "        pass",
        "rounded_and_stated",
    ),
    (
        "rounding on the page",
        RUN,
        '         else fill("checks.rounding_value", **vars(found.rounding))),',
        '         else word("checks.rounding_none")),',
        "rounded_and_stated",
    ),
    (
        "RSEM needs the declaration",
        RSEM,
        "    if raw.format == FORMAT and kind != counts.RSEM_EXPECTED:",
        "    if False:",
        "folder_problems",
    ),
    (
        "expected_count, never TPM",
        RSEM,
        "value=pl.col(EXPECTED).str.strip_chars())",
        'value=pl.col("TPM").str.strip_chars())',
        "folder_combined",
    ),
    (
        "joined on gene, not row",
        RSEM,
        "        columns.append(pl.Series(sample, [by_id.get(i) for i in ids], pl.String))",
        '        columns.append(pl.Series(sample, frame["value"].to_list(), pl.String))',
        "folder_combined",
    ),
    (
        "a gene missing is a problem",
        RSEM,
        "        if given != first:",
        "        if False:",
        "folder_problems",
    ),
    (
        "sample pattern used",
        RSEM,
        "    named = [(sample_name(p, spec.sample_pattern), p) for p in paths]",
        '    named = [(sample_name(p, ""), p) for p in paths]',
        "folder_combined",
    ),
    (
        "unmatched name a problem",
        RSEM,
        '    problems: list[Problem] = [("rsem", p.name, said) for s, p in named if s is None]',
        "    problems: list[Problem] = []",
        "folder_problems",
    ),
    (
        "kind schema",
        "src/ctx/dataset.py",
        "    if kind not in COUNT_KINDS:",
        "    if False:",
        "kind_and_pattern",
    ),
    (
        "pattern needs a group",
        "src/ctx/dataset.py",
        "            if re.compile(pattern).groups < 1:",
        "            if False:",
        "kind_and_pattern",
    ),
    (
        "ctx run takes a folder",
        "src/ctx/pipeline.py",
        "    if options.counts is not None and not options.counts.exists():",
        "    if options.counts is not None and not options.counts.is_file():",
        "takes_an_rsem_folder",
    ),
    (
        "gather takes a folder",
        RUN,
        "    if not matrix.is_file() and not matrix.is_dir():",
        "    if not matrix.is_file():",
        "takes_an_rsem_folder",
    ),
    (
        "folder digest",
        RUN,
        '    lines = "".join(f"{_sha256(f)}  {f.name}\\n" for f in sorted(path.iterdir()) if f.is_file())',
        '    lines = "".join(f"{f.name}\\n" for f in sorted(path.iterdir()) if f.is_file())',
        "takes_an_rsem_folder",
    ),
    # HGNC
    (
        "GENCODE first",
        COUNTS,
        '        through = kind == "entrez" and key not in targets[kind] and key in targets["hgnc"]',
        '        through = kind == "entrez" and key in targets["hgnc"]',
        "hgnc_places",
    ),
    (
        "HGNC one gene only",
        COUNTS,
        '            genes = sorted(targets["hgnc"][key])',
        '            genes = sorted(targets["hgnc"][key])[:1]',
        "hgnc_places",
    ),
    (
        "HGNC gene must be annotated",
        COUNTS,
        '    by["hgnc"] = {key: {g for e in stable for g in by["ensembl"].get(e, set())}',
        '    by["hgnc"] = {key: {g for e in stable for g in by["ensembl"].get(e, {e})}',
        "hgnc_places",
    ),
    (
        "HGNC never unseats GENCODE",
        COUNTS,
        "        leave = [r for r in rows if r not in own] if len(own) == 1 else rows",
        "        leave = rows",
        "hgnc_places",
    ),
    (
        "two HGNC rows both leave",
        COUNTS,
        "        leave = [r for r in rows if r not in own] if len(own) == 1 else rows",
        "        leave = ([r for r in rows if r not in own] if len(own) == 1 else rows) if own else rows[1:]",
        "two_hgnc",
    ),
    (
        "HGNC by column name",
        COUNTS,
        "        entrez, ensembl = (header.index(c) for c in HGNC_COLUMNS)",
        "        entrez, ensembl = 2, 4",
        "hgnc_places",
    ),
    (
        "HGNC missing column refused",
        COUNTS,
        "        if missing := [c for c in HGNC_COLUMNS if c not in header]:",
        "        if missing := []:",
        "missing_column",
    ),
    (
        "gather reads HGNC",
        RUN,
        "    mapped = counts.place(raw, gene_types(gtf), entrez, _hgnc(raw, root, spec, releases))",
        "    mapped = counts.place(raw, gene_types(gtf), entrez, None)",
        "reads_the_hgnc",
    ),
    (
        "largest first",
        RUN,
        '    top = left.sort("total", "source_id", descending=[True, False]).head(LARGEST)',
        '    top = left.sort("total", "source_id", descending=[False, False]).head(LARGEST)',
        "largest",
    ),
    # R2.3
    (
        "invalid TOML caught",
        "src/ctx/dataset.py",
        "    except tomllib.TOMLDecodeError as error:",
        "    except KeyError as error:",
        "invalid_toml or no_traceback",
    ),
    (
        "report catches DatasetError",
        "src/ctx/cli.py",
        '        written = hub.run(ctx.obj["settings"], cards, dataset, dataset_file)\n    except DatasetError as error:',
        '        written = hub.run(ctx.obj["settings"], cards, dataset, dataset_file)\n    except KeyError as error:',
        "no_traceback",
    ),
    # R2.2: GEO descriptions and the guide's notes
    (
        "geo descriptions read",
        GEO,
        "            elif row and row[0] in (CHARACTERISTICS, DESCRIPTIONS):",
        "            elif row and row[0] == CHARACTERISTICS:",
        "descriptions",
    ),
    (
        "keyless description",
        GEO,
        "            if not colon and key.strip():",
        "            if False:",
        "descriptions",
    ),
    (
        "guide authors note",
        "src/ctx/report/guide.py",
        '        word("part1c.table_note"), word("part1c.authors"),',
        '        word("part1c.table_note"),',
        "authors_join",
    ),
    (
        "guide RSEM note",
        "src/ctx/report/guide.py",
        '        word("part1c.other"), word("part1c.rsem"),',
        '        word("part1c.other"),',
        "rsem_note",
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
