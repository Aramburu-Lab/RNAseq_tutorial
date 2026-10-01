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


def test_a_relative_outdir_still_reaches_the_tools(tmp_path, monkeypatch):
    """Each tool runs in its own log folder; a relative FASTQ path once pointed nowhere."""
    outdir(tmp_path / "fetchngs")
    monkeypatch.chdir(tmp_path)
    one = sheet().filter(pl.col("run_accession") == "R_se_a")
    files, _ = reads.fastq_files(one, Path("fetchngs"))
    assert Path(files["path"][0]).is_absolute()
    reads.run_tools(files, tmp_path / "logs", threads=1, fastp=reads.Fastp())
    assert (tmp_path / "logs" / "fastp" / "R_se_a.fastp.json").is_file()


def fastq_lengths(path: Path) -> list[int]:
    """Each read's length in a gzipped FASTQ, read here without the tools."""
    with gzip.open(path, "rt") as handle:
        return [len(line.rstrip("\n")) for i, line in enumerate(handle) if i % 4 == 1]


def test_run_fastp_writes_the_reads_it_passes(tmp_path):
    """M3 trims with M2's filter and keeps the reads: as many as fastp's report says passed,
    as long as it says. se_d's planted tail trims every read from 101 to 81 nt.

    Mutations caught: -o or -O left out (no file); read 2 written to read 1's name."""
    root = outdir(tmp_path / "fetchngs")
    trimmed = tmp_path / "trimmed"
    single = reads.run_fastp([str(root / "fastq" / "se_d.fastq.gz")], tmp_path, "se_d",
                             reads.Fastp(cut_tail=True), 2, trimmed=trimmed)  # fmt: skip
    pair = reads.run_fastp([str(root / "fastq" / f"{PAIR}_{m}.fastq.gz") for m in (1, 2)],
                           tmp_path, "pair", reads.Fastp(), 2, trimmed=trimmed)  # fmt: skip
    assert single == [trimmed / "se_d.fastq.gz"]
    assert pair == [trimmed / "pair_1.fastq.gz", trimmed / "pair_2.fastq.gz"]
    for title, files in (("se_d", single), ("pair", pair)):
        report = json.loads((tmp_path / f"{title}.fastp.json").read_text())
        passed = report["filtering_result"]["passed_filter_reads"]
        lengths = [fastq_lengths(f) for f in files]
        assert {len(n) for n in lengths} == {passed // len(files)}, title
        after = report["summary"]["after_filtering"]["read1_mean_length"]
        assert round(sum(lengths[0]) / len(lengths[0])) == after, title
    assert set(fastq_lengths(single[0])) == {81}
    # The mates are equally long, so only their names tell read 1's file from read 2's.
    for mate, path in zip((1, 2), pair, strict=True):
        with gzip.open(path, "rt") as handle:
            names = [line.rstrip("\n") for i, line in enumerate(handle) if i % 4 == 0]
        assert names and all(name.endswith(f"/{mate}") for name in names), path.name


def test_run_fastp_reaches_relative_folders(tmp_path, monkeypatch):
    """fastp runs in its report folder, so a path relative to the caller pointed nowhere: the
    trimmed reads failed the run, and a missed report path passed silently.

    Mutations caught: the input, the report folder or the trimmed folder passed as given."""
    monkeypatch.chdir(tmp_path)
    shutil.copy(SE / "se_a.fq", tmp_path)
    Path("logs").mkdir()
    written = reads.run_fastp(["se_a.fq"], Path("logs"), "se_a", reads.Fastp(), 1,
                              trimmed=Path("trimmed"))  # fmt: skip
    assert written == [tmp_path / "trimmed" / "se_a.fastq.gz"]
    assert len(fastq_lengths(written[0])) == 1000
    assert (tmp_path / "logs" / "se_a.fastp.json").is_file()


def test_m2_writes_no_reads(ran):
    """M2's fastp runs filter and report; they never write reads.

    Mutation caught: M2's run_tools passing trimmed= (its fastp folder gains FASTQ files)."""
    out = ran[3]
    assert list((out / "logs" / "fastp").glob("*.json"))
    assert not list((out / "logs" / "fastp").glob("*.f*q*"))
    for report in (out / "logs" / "fastp").glob("*.json"):
        assert " -o " not in json.loads(report.read_text())["command"]
