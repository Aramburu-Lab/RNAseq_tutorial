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
