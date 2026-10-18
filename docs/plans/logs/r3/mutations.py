"""R3: plant each mutation of a new rule in place, run the tests named for it, restore.

Run through env.sh: python docs/plans/logs/r3/mutations.py; the exit code counts survivors.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[4]
BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
ENV = {**os.environ, "PATH": f"{BIN}:{os.environ['PATH']}", "PYTHONPATH": str(W / "src"),
       "PYTHONDONTWRITEBYTECODE": "1"}  # fmt: skip
STUDY = "tests/test_proteo_study.py"
PREP = "tests/test_proteo_prepare.py"
Q, A, I, S = (f"src/ctx/proteo/{n}.py" for n in ("quantities", "adapters", "identifiers", "study"))
UNIT = f"{STUDY}::test_missing_never_zero_and_undeclared_zero_refused"
MAP = f"{STUDY}::test_mapping_reasons_in_order"
MQ = f"{STUDY}::test_maxquant_drops_each_flag_once_in_order"
E2E = f"{STUDY}::test_study_end_to_end"
MUTATIONS = [
    ("declared token read as a value", Q, "np.where(token.to_numpy(), np.nan, number",
     "np.where(False, np.nan, number", [UNIT]),
    ("undeclared zero accepted", Q, "np.sum(values <= 0)", "np.sum(values < 0)", [UNIT]),
    ("text cell accepted", Q, "    if text.len():\n", "    if False:\n", [UNIT]),
    ("log10 left as log10", Q, "10: math.log2(10)", "10: 1.0", [UNIT]),
    ("problems raised one at a time", Q,
     '        problems.append(f"no identifier column {id_col!r} in the protein table")\n',
     '        raise ProteoError(f"no identifier column {id_col!r} in the protein table")\n',
     [f"{STUDY}::test_reader_lists_absent_columns_together"]),
    ("collapse keeps the first row", "src/ctx/proteo/tables.py",
     "if name not in best or present[i] > present[best[name]]:", "if name not in best:",
     [f"{STUDY}::test_collapse_after_mapping_keeps_the_best_quantified_row"]),
    ("several genes resolved to one", I, "elif len(found) > 1:", "elif len(found) > 2:", [MAP]),
    ("shared symbol taken", I, "if len(rows) != 1:", "if not rows:", [MAP]),
    ("a member without a map ignored", I, "elif None in keys:", "elif all(k is None for k in keys):",
     [MAP]),
    ("isoform suffix kept", I, 'ISOFORM.sub("", fasta.group(1) if fasta else text)',
     "(fasta.group(1) if fasta else text)", [MAP]),
    ("Ensembl version kept", I, 'VERSION.sub("", value) if kind == "ensembl"',
     'value if kind == "ensembl"', [MAP]),
    ("MaxQuant flags counted twice", A, "hit = keep & (pl.col(column) == FLAGGED)",
     "hit = (pl.col(column) == FLAGGED)", [MQ]),
    ("pre-1.6 Contaminant ignored", A, '("Potential contaminant", "Contaminant")',
     '("Potential contaminant",)', [MQ]),
    ("flagged rows kept", A, "frame = frame.filter(keep).drop(present)",
     "frame = frame.drop(present)", [MQ, E2E]),
    ("DIA-NN run named by its path", A, "    return RUN_ENDINGS.sub(\"\", name)",
     "    return column", [f"{STUDY}::test_diann_names_runs_by_file_and_refuses_twins"]),
    ("plain table read unquoted", A, "names = header(data, separator, quoted=True)",
     "names = header(data, separator, quoted=False)",
     [f"{STUDY}::test_plain_table_honours_quotes"]),
    ("empty cell read as null", Q, 'return frame.with_columns(pl.all().fill_null(""))',
     "return frame", [f"{STUDY}::test_text_frame_keeps_empty_cells"]),
    ("zip read as tar", Q, "    if not zipfile.is_zipfile(path):\n", "    if True:\n",
     [f"{STUDY}::test_archive_member_from_zip"]),
    ("sample table twin accepted", Q, "Counter(v for v in values if v).items() if k > 1",
     "Counter(v for v in values if v).items() if k > 2",
     [f"{STUDY}::test_sample_table_problems_together"]),
    ("RNA sample not gathered accepted", S, "    if not unknown:\n        return []",
     "    if True:\n        return []", [f"{STUDY}::test_study_rna_sample_not_gathered"]),
    ("unpaired protein sample correlated", S, "pairs = [(j, r) for j, r in enumerate(samples.rna) if r]",
     "pairs = [(j, r or samples.rna[0]) for j, r in enumerate(samples.rna)]", [E2E]),
    ("F1 coloured by 24 runs", S, "if name and annotation[name].n_unique() <= len(COLOURS):",
     "if name:", [f"{PREP}::test_template_block_runs_ctx_proteo"]),
    ("UniProt without a map accepted", "src/ctx/dataset.py",
     'if kind == "uniprot" and (id_map not in known', 'if False and (id_map not in known',
     [f"{STUDY}::test_dataset_table_problems_together"]),
    ("table.file beside a resource accepted", "src/ctx/dataset.py", "    if own and \"protein\" in raw:",
     "    if False:", [f"{STUDY}::test_dataset_table_problems_together"]),
    ("PRIDE FTP kept", "src/ctx/proteo/prepare.py", 'ftp[0].replace(FTP, HTTPS, 1)', "ftp[0]",
     [f"{PREP}::test_fetch_checks_size_and_sha1"]),
    ("PRIDE SHA-1 ignored", "src/ctx/proteo/prepare.py",
     'or (entry["sha1"] and sha1.hexdigest() != entry["sha1"])', "",
     [f"{PREP}::test_fetch_checks_size_and_sha1"]),
    ("guide block typed, not the template's", "src/ctx/proteo/guide.py",
     '"block": tuple(pipeline.block("[proteo]", slug=spec.slug).splitlines()),',
     '"block": ("[proteo]",),', [f"{PREP}::test_guide_part_lines_parse"]),
]  # fmt: skip


def clean() -> None:
    """Remove bytecode, so no mutated .pyc outlives its source (lessons.md)."""
    for d in W.rglob("__pycache__"):
        shutil.rmtree(d, ignore_errors=True)


def main() -> int:
    """Plant, test, restore each mutation; the exit code is the number that survived."""
    survived = 0
    for name, rel, old, new, tests in MUTATIONS:
        path = W / rel
        source = path.read_text()
        assert source.count(old) == 1, (name, source.count(old))
        path.write_text(source.replace(old, new))
        clean()
        try:
            done = subprocess.run(
                [f"{BIN}/pytest", "-q", "-x", "-p", "no:cacheprovider", *tests],
                cwd=W, env=ENV, capture_output=True, text=True, check=False,
            )  # fmt: skip
        finally:
            path.write_text(source)
            clean()
        caught = done.returncode != 0
        survived += not caught
        last = (done.stdout.strip().splitlines() or ["(no output)"])[-1]
        print(f"{'caught' if caught else 'SURVIVED'}\t{name}\t{last}", flush=True)
    print(f"{len(MUTATIONS) - survived} of {len(MUTATIONS)} caught")
    return survived


if __name__ == "__main__":
    sys.exit(main())
