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
M2 = "tests/test_figures_m2.py"
M3 = "tests/test_figures_m3.py::test_every_figure_redraws_byte_for_byte"
M4 = "tests/test_figures_m4.py::test_fates_and_estimates_redraw_byte_for_byte"
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
    "parsed": [
        ("parsed description not read back", "src/ctx/figures/kinds.py",
         "if params is not None and k.what(params) == what:", "if params is not None:",
         f"{FIG}::test_a_kind_parsed_from_its_description_must_read_back"),
    ],
    "m3": [
        ("bar outcomes in reverse", "src/ctx/process/figure_tables.py",
         "for label, r in reads.items() for k, v in r.items()]",
         "for label, r in reversed(list(reads.items())) for k, v in r.items()]", M3),
        ("settings from the samples' table", "src/ctx/figures/draw_m3.py",
         'settings = _table(inputs, "setting_outcomes")', 'settings = _table(inputs, "outcomes")',
         M3),
        ("gene sums left out", "src/ctx/figures/draw_m3.py",
         'EMPTY, _table(inputs, "sums")),', "EMPTY, EMPTY),", M3),
        ("groups dropped", "src/ctx/figures/draw_m3.py", '_table(inputs, "groups"), row',
         '_table(inputs, "groups").select("sample"), row', M3),
        ("registry of another module", "src/ctx/figures/kinds_m3.py",
         '"registry": registry_of("process")', '"registry": registry_of("qc")', M3),
        ("registry always qc's", "src/ctx/figures/record.py",
         'options = newest(root / "outputs" / FOLDERS[what], spec.slug, "figure_options")',
         'options = newest(root / "outputs" / FOLDERS["qc"], spec.slug, "figure_options")', M3),
    ],
    "m4": [
        ("contrast levels the other way round", "src/ctx/figures/draw_m4.py",
         "(used.height - tested, tested)", "(tested, used.height - tested)", M4),
        ("fold changes written after the figure", "src/ctx/quant/page.py",
         '                                        lambda f: writer.table(f, "fold_changes", '
         "described(f)))\n",
         "                                        None)\n"
         '        writer.table(_, "fold_changes", described(_))\n', M4),
        ("the whole sheet as the contrast's samples", "src/ctx/quant/figure_tables.py",
         'named["contrast_samples"] = blk.contrast_samples(ctx.inputs)',
         "import polars as pl  # noqa: PLC0415\n"
         '        named["contrast_samples"] = ctx.inputs.samples.select("sample", '
         'group=pl.lit("tested"))',
         "tests/test_figures_m4.py::test_the_contrast_samples_are_m5s"),
        ("sheet table without its groups", "src/ctx/quant/figure_tables.py",
         'named = {"figure_sheet": found.sheet}',
         'named = {"figure_sheet": found.sheet.select("sample")}', M4),
    ],
    "m2": [
        ("pca shares swapped", "src/ctx/figures/draw_m2.py", "(share[0], share[1])",
         "(share[1], share[0])", f"{M2}::test_every_count_figure_redraws_byte_for_byte"),
        ("k from the registry file", "src/ctx/figures/draw_m2.py",
         'registry = replace(reg.load(inputs["registry"]), mad_multiplier=row["mad"])',
         'registry = reg.load(inputs["registry"])',
         f"{M2}::test_the_registry_named_by_the_run_is_the_one_copied"),
        ("sheet without read-length groups", "src/ctx/figures/draw_m2.py",
         'sheet = checkpoint.length_groups(gdc.chosen(pl.read_parquet(inputs["samplesheet"])),',
         'sheet = (gdc.chosen(pl.read_parquet(inputs["samplesheet"])),',
         f"{M2}::test_every_count_figure_redraws_byte_for_byte"),
        ("strip spec where reads take its stem", "src/ctx/figures/kinds_m2.py",
         '_with("metrics", "bands", "flags"), unless=READS),',
         '_with("metrics", "bands", "flags")),',
         f"{M2}::test_reads_figures_redraw_and_shared_stems_have_no_spec"),
        ("theoretical GC curve dropped", "src/ctx/qc/figure_tables.py",
         "None if theoretical.is_empty() else theoretical", "None",
         f"{M2}::test_fastqc_plots_redraw_from_their_tables"),
        ("registry resolved afresh", "src/ctx/figures/record.py",
         '            found[name] = Path(pl.read_parquet(options)["registry"][0])\n',
         "            found[name] = registry.resolve(None)\n",
         f"{M2}::test_a_run_with_its_own_registry_names_that_registry"),
        ("grid run not recorded", "src/ctx/qc/run.py",
         "seen.grid_run if seen is not None else None", "None",
         f"{M2}::test_reads_figures_redraw_and_shared_stems_have_no_spec"),
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
