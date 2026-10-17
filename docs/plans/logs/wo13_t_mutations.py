"""wo13-t: plant each mutation in place, run the tests named for it, restore; caught or not."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[3]
BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
ENV = {
    **os.environ,
    "PATH": f"{BIN}:{os.environ['PATH']}",
    "PYTHONPATH": str(W / "src"),
    "PYTHONDONTWRITEBYTECODE": "1",
}
DOOR = "tests/test_gdc_door.py"
KIT = "tests/test_kit_datasets.py::test_readme_explains_a_git_show_line_exactly_where_one_is"
GDC = "src/ctx/gather/gdc.py"
BUILD = "scripts/python/build_kit.py"
QUANT_IF = "    if refused := gdc.refuse_sheet(sheet, path, MODULE):"
SENTENCE_IF = '    if entry.get("dataset_file_commit"):'
SENTENCE_PATH = '        path = entry["dataset_file_path"]'
MUTATIONS = [
    ("chosen absent passes", GDC, '    if "chosen" not in sheet.columns:\n        found',
     '    if "chosen" not in sheet.columns:\n        return sheet\n        found', [DOOR]),
    ("marker file_id dropped", GDC, '"drop_reason", "file_id", "md5sum"',
     '"drop_reason", "md5sum"', [DOOR]),
    ("marker case_id dropped", GDC, '"md5sum", "case_id")', '"md5sum")', [DOOR]),
    ("sample type ignored on fetchngs", GDC, "        if sample_type is not None:\n",
     "        if False:\n", [DOOR]),
    ("fetchngs sheet filtered", GDC, "        return sheet\n    if",
     "        return sheet.head(1)\n    if", [DOOR]),
    ("quant refusal dropped", "src/ctx/quant/run.py", QUANT_IF,
     QUANT_IF.replace("if refused", "if False and (refused").replace("):", ")):"), [DOOR]),
    ("process refusal dropped", "src/ctx/process/run.py", "    " + QUANT_IF,
     "    " + QUANT_IF.replace("if refused", "if False and (refused").replace("):", ")):"),
     [DOOR]),
    ("genecard preview bypasses door", "scripts/python/preview_genecard.py",
     "    return gdc.chosen(pl.read_parquet(sheet_file))",
     "    return (pl.read_parquet(sheet_file))", [DOOR]),
    ("proteo preview bypasses door", "scripts/python/preview_proteo.py",
     "    sheet = gdc.chosen(pl.read_parquet(design_sheet(root)))",
     "    sheet = (pl.read_parquet(design_sheet(root)))", [DOOR]),
    ("kit sentence never written", BUILD, SENTENCE_IF, "    if False:", [KIT]),
    ("kit sentence on every entry", BUILD, f"{SENTENCE_IF}\n{SENTENCE_PATH}",
     f'    if True:\n{SENTENCE_PATH} or "conf/datasets/x.toml"', [KIT]),
    ("qc door bypassed", "src/ctx/qc/run.py",
     '    sheet = gdc.chosen(pl.read_parquet(inputs["samplesheet"]))',
     '    sheet = pl.read_parquet(inputs["samplesheet"])', ["tests/test_qc_gdc_door.py"]),
    ("de door bypassed", "src/ctx/de/run.py", "    sheet = gdc.chosen(sheet)\n",
     "    sheet = sheet\n", ["tests/test_ncrna_tcga.py"]),
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
    return survived


if __name__ == "__main__":
    sys.exit(main())
