"""`ctx qc` end to end on synthetic gather outputs with planted problems."""

import html
import json
import os
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

import polars as pl
import pytest
from conftest import once, writable_copy
from synth import DE_GROUPS, design_counts, write_gather_outputs

from ctx.io.naming import lint
from ctx.report.builder import unpacked

ROOT = Path(__file__).resolve().parents[1]
CTX = Path(sys.executable).parent / "ctx"
SPEC = ROOT / "tests" / "fixtures" / "datasets" / "de_synthetic.toml"
ENV = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
# (gene_id, gene_name, gene_type, chrom) of the genes the mini annotation types specially.
TYPED = {
    "g00001": ("MT-RNR2", "Mt_rRNA", "chrM"), "g00002": ("MT-CO1", "protein_coding", "chrM"),
    "g00003": ("RNA5-8SN1", "rRNA", "chr21"), "g00004": ("HBB", "protein_coding", "chr11"),
}  # fmt: skip


def ctx(*args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(CTX), *args], cwd=cwd, env=ENV, capture_output=True, text=True,
                          check=False)  # fmt: skip


def plant(root: Path) -> dict[str, str]:
    """Gather outputs with a shallow library, a swapped sex label and an MT-rRNA-heavy sample."""
    synthetic = design_counts(DE_GROUPS)
    write_gather_outputs(root, "desynth", synthetic, date(2026, 9, 30))
    gathered = root / "outputs" / "01_gather"
    counts_path = next(gathered.glob("*_counts_v1.parquet"))
    sheet_path = next(gathered.glob("*_samplesheet_v1.parquet"))
    counts, sheet = pl.read_parquet(counts_path), pl.read_parquet(sheet_path)
    shallow, swapped, mito = "S002", "S004", "S020"
    counts = counts.with_columns(
        (pl.col(shallow) // 10).alias(shallow),
        pl.when(pl.col("gene_id") == "g00001").then(pl.col(mito) + 400_000)
        .otherwise(pl.col(mito)).alias(mito),
    )  # fmt: skip
    other = "M" if sheet.filter(pl.col("sample") == swapped)["sex"][0] == "F" else "F"
    swap = pl.when(pl.col("sample") == swapped).then(pl.lit(other)).otherwise("sex")
    sheet = sheet.with_columns(sex=swap)
    counts.write_parquet(counts_path)
    sheet.write_parquet(sheet_path)
    return {"shallow": shallow, "swapped": swapped, "mito": mito}


def with_annotation(root: Path, tmp: Path) -> Path:
    """The synthetic dataset file plus an annotation_gtf resource, and that GTF in place."""
    genes = pl.read_parquet(next((root / "outputs" / "01_gather").glob("*_genes_v1.parquet")))
    lines = []
    for gene_id, symbol in genes.rows():
        name, kind, chrom = TYPED.get(gene_id, (symbol, "protein_coding", "chr1"))
        lines.append(f'{chrom}\tS\tgene\t1\t9\t.\t+\t.\tgene_id "{gene_id}"; '
                     f'gene_type "{kind}"; gene_name "{name}";\n')  # fmt: skip
    gtf = root / "resources" / "gencode" / "mini.gtf"
    gtf.parent.mkdir(parents=True)
    gtf.write_text("".join(lines))
    spec = tmp / "desynth.toml"
    spec.write_text(SPEC.read_text() + '\n[[resources]]\nname = "annotation_gtf"\n'
                    'urls = ["http://127.0.0.1/mini.gtf"]\nbytes = 0\nsha256 = "' + "0" * 64
                    + '"\ndest = "resources/gencode/mini.gtf"\n')  # fmt: skip
    return spec


@pytest.fixture(scope="session")
def ran(tmp_path_factory):
    """ctx qc on the planted gather outputs with a mini annotation: (root, planted, tables,
    page)."""
    return once(tmp_path_factory, "qc_planted", _planted_run)


def _planted_run(folder: Path):
    root = folder / "root"
    root.mkdir()
    planted = plant(root)
    (folder / "spec").mkdir()
    spec = with_annotation(root, folder / "spec")
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(spec),
               cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "02_qc"
    names = ("metrics", "bands", "flags", "status", "confounding", "sex_check")
    tables = {w: pl.read_csv(next(out.glob(f"*_desynth_{w}_v1.csv"))) for w in names}
    (page,) = out.glob("*_qc_desynth_overview_v1.html")
    return root, planted, tables, unpacked(page.read_text())


def test_the_planted_problems_are_flagged(ran):
    _, planted, tables, _ = ran
    flags = {(s, m) for s, m, _ in tables["flags"].rows()}
    assert (planted["shallow"], "library_size") in flags
    assert (planted["swapped"], "sex_check") in flags
    assert (planted["mito"], "mt_rrna_fraction") in flags


def test_one_status_row_per_run_and_nothing_excluded(ran):
    root, planted, tables, _ = ran
    sheet = pl.read_parquet(next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet")))
    status = tables["status"]
    assert sorted(status["run_accession"]) == sorted(sheet["run_accession"])
    assert set(status["qc_status"]) <= {"pass", "flag"}
    flagged = set(status.filter(pl.col("qc_status") == "flag")["sample"])
    assert set(planted.values()) <= flagged


def test_confounding_covers_the_unpaired_contrasts_with_bh(ran):
    _, _, tables, _ = ran
    tests = tables["confounding"].filter(pl.col("p").is_not_null())
    assert tests.height > 0
    assert (tests["p_adj"] >= tests["p"] - 1e-12).all()
    assert "a_case" in " ".join(tables["confounding"]["contrasts"].to_list())


def test_size_against_genes_marks_the_samples_flagged_on_its_axes(ran):
    """The triangles are the samples flagged on library size or genes detected; the lead
    counts them and the samples flagged only on other metrics, such as the planted swap.

    Mutation caught: marking a sample flagged on any metric (the swap and the MT-rRNA sample
    drawn as triangles)."""
    root, planted, tables, page = ran
    out = root / "outputs" / "02_qc"
    drawn = pl.read_csv(next(out.glob("*_desynth_size_genes_v1.csv")))
    flags = tables["flags"]
    axes = flags.filter(pl.col("metric").is_in(["library_size", "genes_detected"]))
    on_axes = set(axes["sample"])
    elsewhere = set(flags["sample"]) - on_axes
    assert planted["shallow"] in on_axes and {planted["swapped"], planted["mito"]} <= elsewhere
    assert set(drawn.filter(pl.col("flagged"))["sample"]) == on_axes
    by = dict(drawn.select("sample", pl.col("flagged_by").fill_null("")).rows())
    assert {s for s, metrics in by.items() if metrics} == set(flags["sample"])
    assert "sex_check" in by[planted["swapped"]].split(", ")
    lead = (f"{len(on_axes)} of {drawn.height} samples are flagged on library size or genes "
            f"detected, drawn as triangles; {len(elsewhere)} more are flagged only on other "
            "metrics.")  # fmt: skip
    assert lead in page


def test_page_and_outputs(ran):
    root, _, tables, page = ran
    assert lint(root / "outputs" / "02_qc")[1] == []
    assert "M2 Quality control" in page
    runs = tables["status"].height
    flagged = tables["status"].filter(pl.col("qc_status") != "pass").height
    assert f"Flagged runs: {flagged} of {runs}" in page
    assert "built-in method" not in page
    assert "mt_rrna_fraction" in tables["metrics"].columns


def test_the_composition_scatter_plots_both_shares(ran):
    """P2Q.3: MT % against nuclear rRNA %, the mito-planted sample far up.

    Mutation caught: the composition figure left out of the identity group."""
    root, planted, tables, page = ran
    assert "Mitochondrial against nuclear rRNA share" in page
    (drawn,) = (root / "outputs" / "02_qc").glob("*_qc_desynth_composition_v1.csv")
    frame = pl.read_csv(drawn)
    shares = dict(frame.select("sample", "mt_fraction").rows())
    metrics = dict(tables["metrics"].select("sample", "mt_fraction").rows())
    assert shares == pytest.approx({s: 100 * v for s, v in metrics.items()})


def test_a_second_run_writes_every_table_byte_for_byte(ran, tmp_path):
    """S.8: the same inputs give the same bytes, row order included.

    linked_pca, size_genes and composition changed row order run to run: a sample table made
    by an unordered ``unique`` met an inner join whose order polars leaves unspecified, and
    a figure's colour levels followed that order too. Every qc join now keeps its left
    table's order. Mutations caught: the unordered table and join put back in
    ``figures.linked_pca``, and in ``figures.density``; either one alone survives in density,
    where the other still holds the order."""
    root, _, _, _ = ran
    again = tmp_path / "again"
    writable_copy(root / "outputs" / "01_gather", again / "outputs" / "01_gather")
    spec = with_annotation(again, tmp_path)
    done = ctx("--root", str(again), "qc", "--dataset", "desynth", "--dataset-file", str(spec),
               cwd=again)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    first = {p.name: p.read_bytes() for p in (root / "outputs" / "02_qc").glob("*.csv")}
    second = {p.name: p.read_bytes() for p in (again / "outputs" / "02_qc").glob("*.csv")}
    assert first.keys() == second.keys() and len(first) > 10  # noqa: PLR2004
    differ = [n for n in sorted(first) if first[n] != second[n]]
    shown = {n: [(a, b) for a, b in zip(first[n].splitlines(), second[n].splitlines(),
                                        strict=False) if a != b][:2] for n in differ}  # fmt: skip
    assert differ == [], shown


def test_the_m3_alignment_table_is_judged_on_the_align_stage(tmp_path):
    """M3.12: with `ctx process`'s tables beside gather's, M2 judges them with dataset bands.

    Mutation caught: align_tables finding nothing (no unique_pct column, S003 unflagged)."""
    root = tmp_path / "root"
    root.mkdir()
    plant(root)
    sheet = pl.read_parquet(next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet")))
    samples = sheet.unique("sample", keep="first", maintain_order=True)["sample"].to_list()
    unique = [1.0 if s == "S003" else 50.0 + (i % 5) * 0.5 for i, s in enumerate(samples)]
    process = root / "outputs" / "03_process"
    process.mkdir(parents=True)
    empty = pl.Series([None] * len(samples), dtype=pl.Float64)
    table = pl.DataFrame({"sample": samples, "unique_pct": unique,
                          "multi_pct": [2.0] * len(samples), "mismatch_rate": empty,
                          "splices_per_100_unique": empty})  # fmt: skip
    table.write_parquet(process / "2026-10-01_process_desynth_alignment_v1.parquet")
    pl.DataFrame({"sample": samples, "forward_share": [0.5] * len(samples)}).write_parquet(
        process / "2026-10-01_process_desynth_strandedness_v1.parquet")  # fmt: skip
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "02_qc"
    metrics = pl.read_csv(next(out.glob("*_metrics_v1.csv")))
    assert {"unique_pct", "forward_share"} <= set(metrics.columns)
    flags = pl.read_csv(next(out.glob("*_flags_v1.csv")))
    assert ("S003", "align_unique_pct") in set(flags.select("sample", "metric").rows())
    assert "Reads mapped once" in unpacked(next(out.glob("*_overview_v1.html")).read_text())


def test_without_an_annotation_the_shares_are_left_out_and_said(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    plant(root)
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "02_qc"
    metrics = pl.read_csv(next(out.glob("*_metrics_v1.csv")))
    assert "mt_fraction" not in metrics.columns and "library_size" in metrics.columns
    assert "top_gene_share" not in metrics.columns
    assert "Composition shares and the top-gene share not measured" in unpacked(
        next(out.glob("*_overview_v1.html")).read_text()
    )


def test_ctx_qc_leaves_the_registry_listed_lncrna_and_mt_rrna_out_of_the_top_gene(tmp_path):
    """End to end: a gene renamed to ENSG00000280614.1 with the largest count, and the planted
    MT-RNR2 surplus, are not S020's top gene; the lncRNA counts as nuclear rRNA.

    Mutation caught: counts_stage passing no exclusion list, as a run that ignores the registry
    would (S020's top gene then the lncRNA)."""
    root = tmp_path / "root"
    root.mkdir()
    planted = plant(root)["mito"]
    gathered = root / "outputs" / "01_gather"
    lnc = "ENSG00000280614.1"
    for what in ("counts", "genes"):
        path = next(gathered.glob(f"*_{what}_v1.parquet"))
        frame = pl.read_parquet(path).with_columns(
            gene_id=pl.when(pl.col("gene_id") == "g00005").then(pl.lit(lnc))
            .otherwise(pl.col("gene_id")))  # fmt: skip
        if what == "counts":
            frame = frame.with_columns(
                pl.when(pl.col("gene_id") == lnc).then(pl.col(planted) + 900_000)
                .otherwise(pl.col(planted)).alias(planted))  # fmt: skip
            counts = frame
        frame.write_parquet(path)
    spec = with_annotation(root, tmp_path)
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(spec),
               cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    metrics = pl.read_csv(next((root / "outputs" / "02_qc").glob("*_metrics_v1.csv")))
    row = metrics.filter(pl.col("sample") == planted).row(0, named=True)
    column = counts.filter(~pl.col("gene_id").is_in([lnc, "g00001"]))[planted]
    size = counts[planted].sum()
    assert row["top_gene_share"] == pytest.approx(column.max() / size)
    lnc_reads = counts.filter(pl.col("gene_id") == lnc)[planted][0]
    rrna = counts.filter(pl.col("gene_id") == "g00003")[planted][0]
    assert row["rrna_fraction"] == pytest.approx((lnc_reads + rrna) / size)


def test_a_strata_column_the_sheet_lacks_stops_the_run(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    plant(root)
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               "--strata", "tissue,nowhere", cwd=root)  # fmt: skip
    assert done.returncode == 1
    assert "nowhere" in done.stderr


RNASEQ = ROOT / "tests" / "fixtures" / "rnaseq_star"
RNASEQ_SAMPLES = ("RAP1_IAA_30M_REP1", "RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2", "WT_REP1",
                  "WT_REP2")  # fmt: skip


@pytest.fixture(scope="session")
def with_reports(tmp_path_factory):
    """Gather outputs named after the rnaseq fixture's samples, then ctx qc on its reports;
    once for this file and test_figures_m2."""
    return once(tmp_path_factory, "qc_reports", _with_reports)


def _with_reports(root: Path):
    from synth import Synthetic  # noqa: PLC0415

    synthetic = design_counts((("control", "short", 3), ("case", "short", 2)))
    names = dict(zip(synthetic.samples["sample"], RNASEQ_SAMPLES, strict=True))
    renamed = Synthetic(synthetic.counts.rename(names),
                        synthetic.samples.with_columns(pl.col("sample").replace_strict(names)),
                        synthetic.truth)  # fmt: skip
    write_gather_outputs(root, "desynth", renamed, date(2026, 9, 30))
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               "--rnaseq-dir", str(RNASEQ), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return root, root / "outputs" / "02_qc"


def test_reports_add_fates_and_fastqc_figures_from_nfcore_tables(with_reports):
    """The read stage is nf-core's own MultiQC tables: the fates, FastQC's plots per sample,
    and no MultiQC run here (D-116). Mutation: MultiQC run again over the outdir."""
    root, out = with_reports
    fates = pl.read_csv(next(out.glob("*_desynth_fates_v1.csv")))
    assert sorted(fates["sample"]) == sorted(RNASEQ_SAMPLES)
    assert (fates["num_processed"] + fates["removed_before_salmon"] == fates["reads_written"]).all()
    metrics = pl.read_csv(next(out.glob("*_desynth_metrics_v1.csv")))
    assert {"salmon-percent_mapped", "total_reads", "Total Sequences"} <= set(metrics.columns)
    for what in ("read_fates", "read_flow", "base_quality", "fastqc_adapter", "fastqc_n_content",
                 "gc_distribution", "overrepresented"):  # fmt: skip
        assert len(list(out.glob(f"*_desynth_{what}_v1*"))) == 7, what
    assert not (out / "multiqc").exists() and not list(out.glob("*multiqc_config*"))
    assert lint(out)[1] == []
    page = html.unescape(unpacked(next(out.glob("*_overview_v1.html")).read_text()))
    assert "Where the reads went" in page and "built-in method" not in page
    assert "FastQC's own plots, sample by sample" in page
    readme = next(out.glob("*_desynth_fates_v1_README.md")).read_text()
    assert "multiqc (nf-core)" in readme and "1.35" in readme


# Every sentence a registry text says about a page or to its reader, with the figure it is
# about and the page controls that figure must carry for the sentence to hold. A registry
# text speaks of the page in no other sentence: D5 found fastp's settings called ones "the
# page lets you move" (2026-10-07), and the page had no such control.
PAGE_CLAIMS: dict[str, tuple[str, list[str]]] = {}
ABOUT_THE_PAGE = re.compile(r"\b(?:pages?|you|your)\b", re.IGNORECASE)


def _figures(page: str) -> dict[str, list[str]]:
    """Each figure on a module page, by its file description: the kinds of its controls."""
    found = {}
    for block in page.split('<figure class="fig"')[1:]:
        (name,) = set(re.findall(r"_desynth_(\w+?)_v1_white\.pdf", block))
        controls = re.findall(r'id="ctx-controls-\d+">(.*?)</script>', block, re.S)
        found[name] = [json.loads(c)["kind"] for c in controls]
    return found


def test_no_registry_text_claims_a_control_the_page_lacks(with_reports):
    """Every registry sentence about the page names a figure the page draws, with the controls
    it has; since the fastp grid went (D-116), no registry text speaks of the page.

    Mutation caught: the old rationale, "which the page lets you move", restored."""
    from ctx.qc import registry  # noqa: PLC0415

    _, out = with_reports
    page = html.unescape(unpacked(next(out.glob("*_overview_v1.html")).read_text()))
    figures = _figures(page)
    assert figures["linked_pca"] == ["brush"], "the reader of controls must see a real one"
    loaded = registry.load(ROOT / "conf" / "qc_metrics.toml")
    claims = [(m.id, sentence.rstrip(".;")) for m in loaded.metrics.values()
              for field in (m.title, m.meaning, m.computed, m.rationale, m.source, m.plot_source,
                            *m.causes)
              for sentence in re.split(r"(?<=[.;])\s+", field)
              if ABOUT_THE_PAGE.search(sentence)]  # fmt: skip
    for metric, claim in claims:
        assert claim in PAGE_CLAIMS, f"{metric} speaks of the page unchecked: {claim!r}"
        figure, controls = PAGE_CLAIMS[claim]
        assert figures[figure] == controls, (metric, figure)
        assert claim in page, "the claim must be printed on the page it speaks of"
