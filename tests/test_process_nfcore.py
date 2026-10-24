"""M3 on nf-core's alignment (`ctx process --rnaseq-dir`, D-116): the reader on the aligned
fixture, each with a planted mutation, and the command end to end under full and core.

Every number is nf-core's own; each test names the mutation it catches.
"""

import html
import shutil
from pathlib import Path

import polars as pl
import pytest
from conftest import ctx_run, once, writable_copy
from process_nfcore_root import DATASET, SAMPLES, SPEC, STAR, plant
from test_figures import drawn_alike, redrawn_each, specs

from ctx.figures.spec import parse
from ctx.io.naming import lint
from ctx.process import ProcessError, nfcore
from ctx.process.page import groups_of
from ctx.qc import registry as reg
from ctx.report.builder import unpacked

DATA = STAR / "multiqc" / "star_salmon" / "multiqc_report_data"
# Every figure of nf-core's route, by stem: the nine of its own and the two counter figures.
FIGURES = {"aligned_mapping", "rseqc_strand", "rseqc_share", "read_distribution",
           "junction_annotation", "gene_coverage", "duplication", "library_complexity",
           "biotypes", "counters_salmon", "agreement_salmon"}  # fmt: skip


def _outdir(path: Path = STAR) -> nfcore.Outdir:
    return nfcore.open_outdir(path, SAMPLES)


def _copy(tmp_path: Path) -> Path:
    target = tmp_path / "rnaseq_star"
    shutil.copytree(STAR, target)
    return target


def _table(outdir: Path, name: str) -> Path:
    return outdir / "multiqc" / "star_salmon" / "multiqc_report_data" / name


# --- the reader ------------------------------------------------------------------------------


def test_star_outcomes_are_star_own_and_add_up():
    """Each sample's outcomes add up to STAR's input, and the percents are STAR's own."""
    star = nfcore.star(_outdir())
    assert star["sample"].to_list() == list(SAMPLES)
    raw = pl.read_csv(DATA / "multiqc_star.txt", separator="\t").row(0, named=True)
    row = star.filter(pl.col("sample") == raw["Sample"]).row(0, named=True)
    assert row["unique_pct"] == raw["uniquely_mapped_percent"]
    assert row["input_reads"] == raw["total_reads"] == sum(row[k] for k in nfcore.OUTCOMES)


def test_star_outcomes_missing_a_read_are_refused(tmp_path):
    """Mutation: the outcome check dropped; a sample with a read lost reads as whole."""
    outdir = _copy(tmp_path)
    path = _table(outdir, "multiqc_star.txt")
    frame = pl.read_csv(path, separator="\t", infer_schema=False)
    frame = frame.with_columns(
        pl.when(pl.col("Sample") == "WT_REP2")
        .then(pl.lit("1.0"))
        .otherwise(pl.col("unmapped_other"))
        .alias("unmapped_other")
    )
    frame.write_csv(path, separator="\t")
    with pytest.raises(ProcessError, match="WT_REP2"):
        nfcore.star(_outdir(outdir))


def test_a_sample_missing_from_a_table_is_refused(tmp_path):
    """Mutation: a sample absent from one table drawn without it."""
    outdir = _copy(tmp_path)
    path = _table(outdir, "rseqc_infer_experiment_plot.txt")
    lines = path.read_text().splitlines()
    path.write_text("\n".join(line for line in lines if not line.startswith("WT_REP1\t")) + "\n")
    with pytest.raises(ProcessError, match="WT_REP1: in the sheet, not in RSeQC"):
        nfcore.strand(_outdir(outdir))


def test_strand_shares_are_rseqc_own_and_nfcore_share_holds():
    """RSeQC's sense and antisense agree with nf-core's strand check to its one decimal, and
    rseqc_pct is antisense over sense plus antisense on these reverse samples (the registry's
    claim). Mutation: Sense and Antisense swapped on reading."""
    strand = nfcore.strand(_outdir())
    check = pl.read_csv(DATA / "multiqc_strand_check_summary_table.txt", separator="\t")
    joined = strand.join(check.rename({"Sample": "sample"}), on="sample")
    for r in joined.iter_rows(named=True):
        # nf-core prints one decimal; a value on the half rounds either way (0.75 to 0.7).
        assert abs(r["rseqc_sense"] - r["rseqc_s"]) <= 0.05 + 1e-9
        assert abs(r["rseqc_antisense"] - r["rseqc_a"]) <= 0.05 + 1e-9
        share = 100 * r["rseqc_antisense"] / (r["rseqc_sense"] + r["rseqc_antisense"])
        assert round(share, 1) == r["rseqc_pct"]
    assert set(strand["rseqc_inferred"]) == {"reverse"} and set(strand["status"]) == {"pass"}


def test_strand_columns_swapped_are_caught(tmp_path):
    """The planted swap itself: the shares no longer agree with nf-core's check."""
    outdir = _copy(tmp_path)
    path = _table(outdir, "rseqc_infer_experiment_plot.txt")
    path.write_text(path.read_text().replace("Sense\tAntisense", "Antisense\tSense", 1))
    swapped = nfcore.strand(_outdir(outdir))
    assert (swapped["rseqc_sense"] > swapped["rseqc_antisense"]).all()


def test_a_gene_one_matrix_lacks_is_left_out_and_counted(tmp_path):
    """Salmon on STAR's alignment held 471 genes Salmon alone has no row for (Breast-tissue-BRCA,
    2026-10-24), and M3 stopped. Now the shared genes are compared and the others counted with
    their share of reads, stated in the counter group's lead. Catches the gene compared as a
    zero, the refusal back, the share taken from the wrong matrix, and the lead dropping it."""
    from ctx.process import nfcore_blocks, nfcore_page  # noqa: PLC0415 - the page under test

    outdir = _copy(tmp_path)
    path = outdir / nfcore.MATRICES["salmon"]
    alone = pl.read_csv(path, separator="\t")
    on_star = pl.read_csv(outdir / nfcore.MATRICES["star_salmon"], separator="\t")
    sums = on_star.select("gene_id", reads=pl.sum_horizontal(pl.exclude("gene_id", "gene_name")))
    gone = sums.filter(pl.col("reads") > 0).row(0)[0]
    alone.filter(pl.col("gene_id") != gone).write_csv(path, separator="\t")
    long, unshared = nfcore.matrices(_outdir(outdir))
    assert gone not in set(long["gene_id"]) and long["gene_id"].n_unique() == alone.height - 1
    share = sums.filter(pl.col("gene_id") == gone)["reads"].item() / sums["reads"].sum()
    assert unshared["salmon"] == (0, 0.0)
    assert unshared["star_salmon"][0] == 1 and unshared["star_salmon"][1] == pytest.approx(share)
    lead = nfcore_page.counters_lead(nfcore_blocks.Found(unshared=unshared),
                                     {"n_genes": 10, "n_identical": 4})  # fmt: skip
    assert "1 genes only Salmon on STAR's alignment holds" in lead and f"{share:.2%}" in lead


@pytest.mark.parametrize("n", [5, 103, 210])
def test_a_bar_per_sample_keeps_its_caption_clear_at_any_count(n):
    """103 samples drew a 20 in tall mapping figure whose axis title ran into its caption: the
    caption's room left out its gap above the page edge, a share of the height (2026-10-24)."""
    from ctx.process import nfcore_figures  # noqa: PLC0415 - the figure under test
    from ctx.viz.save import layout_problems  # noqa: PLC0415 - heavy
    from ctx.viz.style import THEMES  # noqa: PLC0415 - heavy

    rows = [(f"S{i}", k, 10, 100, 100 / len(nfcore.OUTCOMES))
            for i in range(n) for k in nfcore.OUTCOMES]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["sample", "outcome", "reads", "input_reads", "percent"],
                         orient="row")  # fmt: skip
    data = nfcore_figures.mapping_bars(frame, [], f"Reads by STAR outcome; {n} samples. No test.")
    assert layout_problems(data.static(THEMES["white"])) == []


def test_a_long_curve_is_drawn_thinned_with_its_ends():
    """dupRadar's 3,594 points a sample made a 16.6 MB page for 103 samples (2026-10-24). Each
    curve keeps at most one point over CURVE_POINTS, its first and last, in order; a short one
    keeps every point. Catches the thinning dropped, an end lost, and samples thinned together."""
    from ctx.process.nfcore_blocks import CURVE_POINTS, thinned  # noqa: PLC0415 - under test

    long = pl.DataFrame({"sample": ["a"] * 3594 + ["b"] * 7,
                         "x": [*range(3594), *range(7)], "y": [0.0] * 3601})  # fmt: skip
    kept = thinned(long)
    a, b = kept.filter(pl.col("sample") == "a"), kept.filter(pl.col("sample") == "b")
    assert CURVE_POINTS // 2 < a.height <= CURVE_POINTS + 1 and b.height == 7
    assert (a["x"][0], a["x"][-1]) == (0, 3593) and a["x"].is_sorted()


def test_the_two_matrices_are_not_swapped():
    """a is Salmon alone (salmon/), b Salmon on STAR's alignment (star_salmon/). Mutation: the
    two MATRICES entries swapped."""
    long, _ = nfcore.matrices(_outdir())
    alone = pl.read_csv(STAR / "salmon" / "salmon.merged.gene_counts_length_scaled.tsv",
                        separator="\t")  # fmt: skip
    on_star = pl.read_csv(STAR / "star_salmon" / "salmon.merged.gene_counts_length_scaled.tsv",
                          separator="\t")  # fmt: skip
    sample = SAMPLES[0]
    both = alone.select("gene_id", a=sample).join(on_star.select("gene_id", b=sample), on="gene_id")
    differ = both.filter(pl.col("a") != pl.col("b"))
    assert differ.height > 0, "the two matrices must differ somewhere for the test to see a swap"
    gene, a, b = differ.row(0)
    row = long.filter((pl.col("gene_id") == gene) & (pl.col("sample") == sample)).row(0, named=True)
    assert (row["a"], row["b"]) == (a, b)


def test_curves_and_ragged_tables_read_every_sample():
    outdir = _outdir()
    coverage, duplication = nfcore.coverage(outdir), nfcore.duplication(outdir)
    curves, perfect = nfcore.complexity(outdir)
    for frame in (coverage, duplication, curves):
        assert set(frame["sample"]) == set(SAMPLES)
    assert coverage.group_by("sample").len()["len"].to_list() == [100] * len(SAMPLES)
    assert perfect.height == 2 and "lc_extrap" not in "".join(curves["sample"])


def test_an_outdir_not_aligned_is_refused():
    """Mutation: the Salmon-only fixture read as an alignment; it has no STAR."""
    with pytest.raises(ProcessError, match="skipped alignment"):
        nfcore.open_outdir(STAR.parent / "rnaseq", SAMPLES)


def test_groups_keep_the_sheet_order():
    """ED.7b: the groups table, built twice from one sheet, is the same, in sheet order.
    Mutation: ``unique`` without ``maintain_order``, which ordered samples by hash."""
    names = [f"s{i:03d}" for i in range(200)][::-1]
    sheet = pl.DataFrame({"sample": names, "cell_type": ["a", "b"] * 100,
                          "condition": ["x"] * 100 + ["y"] * 100})  # fmt: skip
    first, second = groups_of(sheet, names), groups_of(sheet, names)
    assert first.equals(second) and first["sample"].to_list() == names


# --- end to end ------------------------------------------------------------------------------


def _built(folder: Path, edition: str) -> tuple[Path, list[str]]:
    root = plant(folder)
    done = ctx_run("--root", str(root), "--edition", edition, "process", "--dataset", DATASET,
                   "--dataset-file", str(SPEC), "--rnaseq-dir", str(STAR), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return root, done.stdout.split()


@pytest.fixture(scope="session")
def full(tmp_path_factory):
    return once(tmp_path_factory, "process_nfcore_full", lambda f: _built(f, "full"))


@pytest.fixture(scope="session")
def core(tmp_path_factory):
    return once(tmp_path_factory, "process_nfcore_core", lambda f: _built(f, "core"))


@pytest.mark.parametrize("edition", ["full", "core"])
def test_every_figure_is_written_once_and_named_by_the_lint(request, edition):
    root, printed = request.getfixturevalue(edition)
    out = root / "outputs" / "03_process"
    assert len(printed) == len(set(printed))
    drawn = {p.name.split(f"{DATASET}_", 1)[1].removesuffix("_v1_white.pdf")
             for p in out.glob("*_v1_white.pdf")}  # fmt: skip
    assert drawn == FIGURES
    assert lint(out)[1] == []
    (page,) = out.glob("*_overview_v1.html")
    assert page.stat().st_size < 15_000_000


def test_the_full_page_holds_every_section_and_no_tool_ran(full):
    root, _ = full
    (page,) = (root / "outputs" / "03_process").glob("*_overview_v1.html")
    text = html.unescape(unpacked(page.read_text()))
    for heading in ("Where each sample's reads align", "Which strand each sample's reads come",
                    "Coverage along transcripts", "Duplication against expression",
                    "How many distinct molecules", "What the reads counted on genes are",
                    "Splice junctions", "Where in the genes"):  # fmt: skip
        assert heading in text, heading
    assert not (root / "outputs" / "03_process" / "tools").exists()
    assert not (root / "outputs" / "03_process" / "multiqc").exists()


def test_the_strand_threshold_is_the_registry_s(full, tmp_path):
    """Mutation: the threshold typed in the code; a registry with another plot line moves the
    line, so the page's line is the registry's own."""
    root, _ = full
    registry = reg.load(reg.resolve(None))
    page = unpacked(next((root / "outputs" / "03_process").glob("*_overview_v1.html"))
                    .read_text())  # fmt: skip
    line = f"stranded threshold, {registry.metrics['rseqc_pct'].plot_line:g} %"
    assert line in page
    moved = tmp_path / "registry.toml"
    moved.write_text(registry.path.read_text().replace("plot_line = 80.0", "plot_line = 70.0"))
    copy = plant(tmp_path / "root")
    done = ctx_run("--root", str(copy), "process", "--dataset", DATASET, "--dataset-file",
                   str(SPEC), "--rnaseq-dir", str(STAR), "--registry", str(moved),
                   cwd=copy)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    page = unpacked(next((copy / "outputs" / "03_process").glob("*_overview_v1.html"))
                    .read_text())  # fmt: skip
    assert "stranded threshold, 70 %" in page and line not in page


def test_laptop_options_beside_rnaseq_dir_are_refused(tmp_path):
    root = plant(tmp_path / "root")
    done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                   str(SPEC), "--rnaseq-dir", str(STAR), "--genome", str(SPEC),
                   cwd=root)  # fmt: skip
    assert done.returncode == 1 and "--genome" in done.stderr
    assert not (root / "outputs" / "03_process").exists()


def test_every_figure_redraws_byte_for_byte(full, tmp_path):
    """Each of the eleven figures, redrawn from M3's tables by ``ctx figure``. Mutation: a
    drawer reading the tables back into another frame."""
    root, _ = full
    found = specs(root, ("03_process",), SPEC)
    names = {parse(s.figure_id)[2] for s in found}
    assert len(found) == len(FIGURES), names
    for spec, csv in zip(found, redrawn_each(found, root, tmp_path), strict=True):
        module_csv = next(root.glob(f"outputs/03_process/{spec.stem}.csv"))
        assert drawn_alike(csv, module_csv) == [], spec.stem


def test_ctx_run_plans_process_on_an_aligned_outdir_only(tmp_path):
    """process runs with --rnaseq-dir where the outdir was aligned with STAR, and not on a
    Salmon-only one. Mutation: the plan's aligned check dropped, which runs M3 on Salmon alone."""
    from ctx import pipeline  # noqa: PLC0415
    from ctx.dataset import load  # noqa: PLC0415

    spec = load(SPEC)
    planned = dict(pipeline.steps(tmp_path, spec, pipeline.Options(DATASET, SPEC,
                                                                   rnaseq_dir=STAR)))  # fmt: skip
    assert planned["process"][-2:] == ("--rnaseq-dir", str(STAR))
    assert planned["qc"][-2:] == ("--rnaseq-dir", str(STAR))
    salmon_only = dict(pipeline.steps(tmp_path, spec, pipeline.Options(
        DATASET, SPEC, rnaseq_dir=STAR.parent / "rnaseq")))  # fmt: skip
    assert "process" not in salmon_only


@pytest.mark.parametrize("edition", ["full", "core"])
def test_qc_reads_nfcore_tables_and_m3_strand_check(request, edition, tmp_path):
    """M2 on the same root after M3, under each edition: nf-core's FastQC, Trim Galore and
    Salmon tables per sample, and M3's strand_check judged on the align stage (D-116).
    Mutation: the align stage left without strand_check (no RSeQC strip)."""
    built, _ = request.getfixturevalue(edition)
    root = writable_copy(built, tmp_path / "root")
    done = ctx_run("--root", str(root), "--edition", edition, "qc", "--dataset", DATASET,
                   "--dataset-file", str(SPEC), "--rnaseq-dir", str(STAR), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "02_qc"
    for what in ("strip_rseqc_pct", "strip_rseqc_sense", "base_quality", "read_flow"):
        assert list(out.glob(f"*_{DATASET}_{what}_v1_white.pdf")), what
    metrics = pl.read_csv(next(out.glob("*_metrics_v1.csv")))
    assert set(metrics["status"]) == {"pass"}
    (page,) = out.glob("*_overview_v1.html")
    assert f'data-edition="{edition}"' in page.read_text()[:500]
    assert lint(out)[1] == []


def test_the_drawn_threshold_line_is_the_registry_s():
    """The rseqc_share figure's line sits at the registry's plot line, whatever it is.
    Mutation: the line typed as 80 in the figure code."""
    from dataclasses import replace  # noqa: PLC0415

    from ctx.figures.recorder import Recorder  # noqa: PLC0415
    from ctx.process import nfcore_blocks  # noqa: PLC0415

    metric = reg.load(reg.resolve(None)).metrics["rseqc_pct"]
    found = nfcore_blocks.Found(strand=nfcore.strand(_outdir()))
    groups = pl.DataFrame({"sample": list(SAMPLES)})
    for line in (70.0, 80.0):
        shown = nfcore_blocks.Shown(groups, None, replace(metric, plot_line=line))
        recorder = Recorder("03_process", DATASET)
        nfcore_blocks.strand_check_figure(recorder, found, shown)
        data, _ = recorder.one("rseqc_share")
        assert data.hline is not None and data.hline[0] == line
