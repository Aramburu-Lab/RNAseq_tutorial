"""`ctx qc` end to end on synthetic gather outputs with planted problems."""

import os
import subprocess
import sys
from datetime import date
from pathlib import Path

import polars as pl
import pytest
from synth import DE_GROUPS, design_counts, write_gather_outputs

from ctx.io.naming import lint

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


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    root = tmp_path_factory.mktemp("qc_root")
    planted = plant(root)
    spec = with_annotation(root, tmp_path_factory.mktemp("spec"))
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(spec),
               cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "02_qc"
    names = ("metrics", "bands", "flags", "status", "confounding", "sex_check")
    tables = {w: pl.read_csv(next(out.glob(f"*_desynth_{w}_v1.csv"))) for w in names}
    (page,) = out.glob("*_qc_desynth_overview_v1.html")
    return root, planted, tables, page.read_text()


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


def test_page_and_outputs(ran):
    root, _, tables, page = ran
    assert lint(root / "outputs" / "02_qc")[1] == []
    assert "M2 Quality control" in page
    runs = tables["status"].height
    flagged = tables["status"].filter(pl.col("qc_status") != "pass").height
    assert f"Flagged runs: {flagged} of {runs}" in page
    assert "built-in method" not in page
    assert "mt_rrna_fraction" in tables["metrics"].columns


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
    assert "Composition shares not measured" in next(out.glob("*_overview_v1.html")).read_text()


def test_a_strata_column_the_sheet_lacks_stops_the_run(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    plant(root)
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               "--strata", "tissue,nowhere", cwd=root)  # fmt: skip
    assert done.returncode == 1
    assert "nowhere" in done.stderr


RNASEQ = ROOT / "tests" / "fixtures" / "rnaseq"
RNASEQ_SAMPLES = ("RAP1_IAA_30M_REP1", "RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2", "WT_REP1",
                  "WT_REP2")  # fmt: skip


@pytest.fixture(scope="module")
def with_reports(tmp_path_factory):
    """Gather outputs named after the rnaseq fixture's samples, then ctx qc on its reports."""
    from synth import Synthetic  # noqa: PLC0415

    root = tmp_path_factory.mktemp("qc_rnaseq")
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


def test_reports_add_fates_figures_and_the_multiqc_report(with_reports):
    root, out = with_reports
    fates = pl.read_csv(next(out.glob("*_desynth_fates_v1.csv")))
    assert sorted(fates["sample"]) == sorted(RNASEQ_SAMPLES)
    assert (fates["num_processed"] == fates["reads_written"]).all()
    metrics = pl.read_csv(next(out.glob("*_desynth_metrics_v1.csv")))
    assert {"salmon-percent_mapped", "trim_galore-tg_total_reads"} <= set(metrics.columns)
    for what in ("read_fates", "read_flow"):
        assert len(list(out.glob(f"*_desynth_{what}_v1*"))) == 7, what
    assert (out / "multiqc" / "multiqc_report.html").is_file()
    assert next(out.glob("*_desynth_multiqc_config_v1.yml")).read_text().startswith("{")
    assert lint(out)[1] == []
    page = next(out.glob("*_overview_v1.html")).read_text()
    assert "Where the reads went" in page and "built-in method" not in page
    assert not (RNASEQ / "multiqc").exists()


SE = ROOT / "tests" / "fixtures" / "chr22_se"
SE_SAMPLES = ("se_a", "se_b", "se_c", "se_d", "se_e")


@pytest.fixture(scope="module")
def with_reads(tmp_path_factory):
    """The chr22_se reads as a fetchngs outdir, gather outputs naming them, ctx qc --reads."""
    import gzip  # noqa: PLC0415
    import shutil  # noqa: PLC0415

    from synth import Synthetic  # noqa: PLC0415

    root = tmp_path_factory.mktemp("qc_reads")
    fastq = root / "data" / "raw" / "desynth" / "fastq"
    fastq.mkdir(parents=True)
    for name in SE_SAMPLES:
        with (
            (SE / f"{name}.fq").open("rb") as src,
            gzip.open(fastq / f"{name}.fastq.gz", "wb") as dst,
        ):
            shutil.copyfileobj(src, dst)
    synthetic = design_counts((("control", "short", 3), ("case", "short", 2)))
    names = dict(zip(synthetic.samples["sample"], SE_SAMPLES, strict=True))
    samples = synthetic.samples.with_columns(pl.col("sample").replace_strict(names))
    renamed = Synthetic(synthetic.counts.rename(names), samples, synthetic.truth)
    write_gather_outputs(root, "desynth", renamed, date(2026, 9, 30))
    sheet_path = next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))
    sheet = pl.read_parquet(sheet_path).with_columns(
        fastq_1=pl.format("fastq/{}.fastq.gz", "sample"), status_1=pl.lit("verified"),
        fastq_2=pl.lit(""), status_2=pl.lit(""))  # fmt: skip
    sheet.write_parquet(sheet_path)
    runs = dict(sheet.select("sample", "run_accession").rows())
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               "--reads", "--grid-run", runs["se_d"], cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return root / "outputs" / "02_qc", runs


def test_reads_are_judged_run_by_run(with_reads):
    out, runs = with_reads
    metrics = pl.read_csv(next(out.glob("*_desynth_read_metrics_v1.csv")))
    assert sorted(metrics["run_accession"]) == sorted(runs.values())
    assert (metrics["fastqc-total_sequences"] == 1000).all()
    flags = {(r, m) for r, m, _ in pl.read_csv(next(out.glob("*_read_flags_v1.csv"))).rows()}
    assert (runs["se_c"], "fastqc_adapter") in flags
    assert (runs["se_d"], "fastqc_per_base_quality") in flags
    status = pl.read_csv(next(out.glob("*_desynth_status_v1.csv")))
    reasons = dict(status.select("run_accession", "qc_reason").rows())
    assert "Adapter content" in reasons[runs["se_c"]]
    assert "Per-base quality" in reasons[runs["se_d"]]


def test_reads_add_curves_the_grid_and_their_own_multiqc_report(with_reads):
    out, runs = with_reads
    for what in ("base_quality", "fastp_grid"):
        assert len(list(out.glob(f"*_desynth_{what}_v1*"))) == 7, what
    grid = pl.read_csv(next(out.glob("*_desynth_grid_v1.csv")))
    assert grid.height == 36 and grid["passed_pct"].n_unique() > 1
    assert (out / "multiqc" / "reads" / "multiqc_report.html").is_file()
    assert lint(out)[1] == []
    page = next(out.glob("*_overview_v1.html")).read_text()
    assert "The reads themselves" in page and "built-in method" not in page


FETCHNGS = ROOT / "tests" / "fixtures" / "fetchngs" / "download"


def test_fetchngs_reads_alone_within_the_track_a_budget(tmp_path):
    """Gather on fetchngs's own test download, then ctx qc --reads with no counts at all."""
    import re  # noqa: PLC0415

    from ctx.runtime import BUDGETS  # noqa: PLC0415

    root = tmp_path / "root"
    spec = tmp_path / "nfcore_test.toml"
    spec.write_text('[dataset]\nslug = "nfcore_test"\n\n[guard]\nids = []\npatterns = []\n')
    for args in (("gather", "--dataset", "nfcore_test", "--fetchngs-dir", str(FETCHNGS)),
                 ("qc", "--dataset", "nfcore_test", "--dataset-file", str(spec), "--reads",
                  "--fetchngs-dir", str(FETCHNGS))):  # fmt: skip
        root.mkdir(exist_ok=True)
        done = ctx("--root", str(root), *args, cwd=root)
        assert done.returncode == 0, (args[0], done.stderr[-3000:])
    out = root / "outputs" / "02_qc"
    page = next(out.glob("*_overview_v1.html")).read_text()
    assert "Counts not judged" in page
    metrics = pl.read_csv(next(out.glob("*_read_metrics_v1.csv")))
    assert metrics.height == 4 and (metrics["fastqc-total_sequences"] > 0).all()
    assert lint(out)[1] == []
    wall = float(re.search(r'<td class="wall">([\d.]+) s</td>', page).group(1))
    rss = [int(re.search(rf'<td class="{c}">[\d.]+ MiB \((\d+) bytes\)</td>', page).group(1))
           for c in ("rss-self", "rss-children")]  # fmt: skip
    assert 0 < wall < BUDGETS["A"].minutes * 60
    assert 0 < sum(rss) < BUDGETS["A"].ram_gb * 10**9
