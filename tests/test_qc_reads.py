"""ctx.qc.reads on the chr22_se fixture and fetchngs's paired run: the files M3's laptop route
reads and fastp, planted problems caught (M2's read stage reads nf-core's tables, D-116)."""

import gzip
import json
import shutil
from pathlib import Path

import polars as pl
import pytest

from ctx.qc import QcError, reads, registry

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


def test_unreadable_runs_are_listed_not_read(tmp_path):
    files, skipped = reads.fastq_files(sheet(), outdir(tmp_path))
    assert skipped.rows() == [("R_bad", "mismatch"), ("R_meta", "not_downloaded")]
    assert files.filter(pl.col("run_accession") == "SRR13191702")["mate"].to_list() == [1, 2]


def test_a_missing_file_stops_the_run(tmp_path):
    root = outdir(tmp_path)
    (root / "fastq" / "se_b.fastq.gz").unlink()
    (root / "fastq" / f"{PAIR}_2.fastq.gz").unlink()
    with pytest.raises(QcError, match=r"2 FASTQ files .*se_b.fastq.gz.*_2.fastq.gz"):
        reads.fastq_files(sheet(), root)


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
