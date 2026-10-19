"""wo13-f1: plant each mutation, run the tests that should catch it, restore the file.

Usage: wo13_f1_py.sh python docs/plans/logs/wo13_f1/mutations.py <module> [name ...]
Each mutation is one exact replacement; the file is restored even when the run fails.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[4]
FIG = "tests/test_figures.py"
EVERY = f"{FIG}::test_every_spec_redraws_its_data_byte_for_byte"
MUTATIONS = {
    "m5m6": [
        ("dispersions not read back", "src/ctx/figures/draw_m5.py",
         "dispersions=pl.read_parquet(dispersions) if dispersions else pl.DataFrame(),",
         "dispersions=pl.DataFrame(),", EVERY),
        ("dispersions not written", "src/ctx/de/run.py",
         '    writer.table(fits[spec.de.featured].dispersions, f"{spec.de.featured}_dispersions",\n'
         "                 DISPERSION_COLUMNS)  # fmt: skip\n", "", EVERY),
        ("pca sheet past the door", "src/ctx/figures/draw_m5.py",
         'sheet = gdc.chosen(pl.read_parquet(inputs["samplesheet"]))',
         'sheet = pl.read_parquet(inputs["samplesheet"])',
         f"{FIG}::test_a_pca_redraw_reads_the_sheet_through_the_gdc_door"),
        ("pca name not checked", "src/ctx/figures/kinds_m5.py",
         'lambda p: f"pca_{p.pca}", _pca_inputs, _pca_check),',
         'lambda p: f"pca_{p.pca}", _pca_inputs),',
         f"{FIG}::test_a_pca_the_dataset_file_does_not_declare_is_refused"),
        ("box genes out of the file's order", "src/ctx/figures/draw_m5.py",
         "symbols = spec.de.box_genes if spec.de else ()",
         "symbols = spec.de.box_genes[::-1] if spec.de else ()", EVERY),
        ("running GSEA with another seed", "src/ctx/figures/draw_m6.py",
         'enriched = {r.name: r for r in m6._gsea(results, lite.read_gmt(inputs["gene_sets"]))}',
         "from ctx.function import running as rank  # noqa: PLC0415\n"
         "    enriched = {r.name: r for r in lite.prerank(*rank.ranking(results), "
         'lite.read_gmt(inputs["gene_sets"]), seed=lite.SEED + 1)}', EVERY),
        ("background of the down calls", "src/ctx/figures/draw_m6.py",
         "m6_background.DIRECTIONS[0]", "m6_background.DIRECTIONS[1]", EVERY),
        ("resource found under the working folder", "src/ctx/figures/record.py",
         "path = root / found[0] if found else None", "path = Path(found[0]) if found else None",
         EVERY),
    ],
}  # fmt: skip


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
