"""ctx.qc.reads on the chr22_se fixture and fetchngs's paired run: planted problems caught."""

import gzip
import json
import shutil
from pathlib import Path

import polars as pl
import pytest

from ctx.qc import QcError, checkpoint, multiqc_io, reads, registry

ROOT = Path(__file__).resolve().parents[1]
SE = ROOT / "tests" / "fixtures" / "chr22_se"
PAIRED = ROOT / "tests" / "fixtures" / "fetchngs" / "download" / "fastq"
PAIR = "SRX9626017_SRR13191702"
REG = registry.load(ROOT / "conf" / "qc_metrics.toml")
SINGLE = ("se_a", "se_b", "se_c", "se_d", "se_e")


def outdir(root: Path) -> Path:
    """A fetchngs-like outdir: the five single-end samples gzipped, and one paired run."""
    fastq = root / "fastq"
    fastq.mkdir(parents=True)
    for name in SINGLE:
        with (
            (SE / f"{name}.fq").open("rb") as src,
            gzip.open(fastq / f"{name}.fastq.gz", "wb") as dst,
        ):
            shutil.copyfileobj(src, dst)
    for mate in (1, 2):
        shutil.copy(PAIRED / f"{PAIR}_{mate}.fastq.gz", fastq)
    return root


def sheet() -> pl.DataFrame:
    rows = [(f"R_{n}", n, f"fastq/{n}.fastq.gz", "verified", "", "") for n in SINGLE]
    rows += [("SRR13191702", "pair", f"fastq/{PAIR}_1.fastq.gz", "no_md5",
              f"fastq/{PAIR}_2.fastq.gz", "verified"),
             ("R_bad", "bad", "fastq/bad.fastq.gz", "mismatch", "", ""),
             ("R_meta", "meta", "", "not_downloaded", "", "")]  # fmt: skip
    return pl.DataFrame(rows, orient="row", schema=["run_accession", "sample", "fastq_1",
                                                    "status_1", "fastq_2", "status_2"])  # fmt: skip


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    root = outdir(tmp_path_factory.mktemp("fetchngs"))
    files, skipped = reads.fastq_files(sheet(), root)
    out = tmp_path_factory.mktemp("outputs") / "02_qc" / "multiqc"
    folders = reads.run_tools(files, out / "logs", threads=2, fastp=reads.Fastp())
    config = out / "config.json"
    config.write_text(json.dumps(registry.multiqc_config(REG, "test")))
    mqc = multiqc_io.run(folders, out, config=config, modules=reads.MODULES)
    return files, skipped, reads.run_values(mqc, files, REG), out


def test_unreadable_runs_are_listed_not_read(ran):
    files, skipped, _, _ = ran
    assert skipped.rows() == [("R_bad", "mismatch"), ("R_meta", "not_downloaded")]
    assert files.filter(pl.col("run_accession") == "SRR13191702")["mate"].to_list() == [1, 2]


def test_a_missing_file_stops_the_run(tmp_path):
    root = outdir(tmp_path)
    (root / "fastq" / "se_b.fastq.gz").unlink()
    (root / "fastq" / f"{PAIR}_2.fastq.gz").unlink()
    with pytest.raises(QcError, match=r"2 FASTQ files .*se_b.fastq.gz.*_2.fastq.gz"):
        reads.fastq_files(sheet(), root)


def test_values_are_the_tools_own(ran):
    _, _, values, out = ran
    stats = pl.read_csv(out / "logs" / "seqkit" / "seqkit_stats.tsv", separator="\t")
    q30 = {reads.stem(f): q for f, q in stats.select("file", "Q30(%)").rows()}
    got = {r["run_accession"]: r for r in values.iter_rows(named=True)}
    for name in SINGLE:
        assert got[f"R_{name}"]["fastqc-total_sequences"] == 1000
        assert got[f"R_{name}"]["seqkit_seqkit-Q30_pct"] == q30[name]
        assert got[f"R_{name}"]["fastp-filtering_result_passed_filter_reads"] == 1000
    # The paired run: FastQC and seqkit from read 1, fastp from the pair's own row.
    assert got["SRR13191702"]["seqkit_seqkit-Q30_pct"] == q30[f"{PAIR}_1"]
    assert got["SRR13191702"]["fastqc-total_sequences"] == 91
    assert got["SRR13191702"]["fastp-pct_surviving_reads"] is not None
    assert set(values["run_accession"]) == {f"R_{n}" for n in SINGLE} | {"SRR13191702"}


def test_planted_problems_are_flagged(ran):
    _, _, values, _ = ran
    judged = values.drop("sample").rename({"run_accession": "sample"})
    units = judged.select("sample", stratum=pl.lit("all"))
    _, flags = checkpoint.judge(judged, units, REG, ["stratum"])
    planted = {"fastqc_adapter", "fastqc_per_base_quality", "q30_rate"}
    found = {(s, m) for s, m, _ in flags.rows() if m in planted}
    # The paired run's read 1 is fetchngs's own test data, at 70 % Q30.
    assert found == {("R_se_c", "fastqc_adapter"), ("R_se_d", "fastqc_per_base_quality"),
                     ("R_se_d", "q30_rate"), ("SRR13191702", "q30_rate")}  # fmt: skip


def test_fastp_grid_runs_each_cell_with_its_own_settings(tmp_path):
    path = SE / "se_d.fq"
    grid = reads.fastp_grid([str(path)], tmp_path / "logs", tmp_path / "multiqc" / "grid",
                            grid=((15, 38), (15, 102)), fastp=reads.Fastp(), threads=2)  # fmt: skip
    kept = {(q, n): r for q, n, r, _ in grid.rows()}
    # Every base is below Q38 ("F" is Q37), and no read is 102 nt long.
    assert kept == {(15, 15): 1000, (15, 102): 0, (38, 15): 0, (38, 102): 0}
    assert grid.filter((pl.col("q") == 15) & (pl.col("l") == 15))["passed_pct"][0] == 100.0
