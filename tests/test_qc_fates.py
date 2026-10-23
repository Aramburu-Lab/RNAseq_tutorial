"""ctx.qc.fates on the aligned fixture: values and fates from nf-core's own MultiQC tables, held
to the tools' own reports (D-116). Each test names the mutation it catches."""

import json
import shutil
from pathlib import Path

import polars as pl
import pytest

from ctx.qc import QcError, fates, multiqc_io, registry

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "rnaseq_star"
DATA = Path("multiqc") / "star_salmon" / "multiqc_report_data"
SAMPLES = ["RAP1_IAA_30M_REP1", "RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2", "WT_REP1", "WT_REP2"]
PAIRED = {"RAP1_IAA_30M_REP1", "WT_REP1", "WT_REP2"}
REG = registry.load(ROOT / "conf" / "qc_metrics.toml")


def copy(tmp: Path) -> Path:
    out = tmp / "rnaseq_star"
    shutil.copytree(FIXTURE, out)
    return out


@pytest.fixture(scope="module")
def ran():
    run = fates.checked(multiqc_io.nfcore(FIXTURE), SAMPLES)
    return fates.values(run, SAMPLES, REG), fates.fates(run, SAMPLES)


def report(sample: str) -> dict:
    """The sample's read-1 Trim Galore JSON report, or its only one."""
    folder = FIXTURE / "trimgalore"
    paths = sorted(folder.glob(f"{sample}_trimmed_1.*.json")) or sorted(
        folder.glob(f"{sample}_trimmed.*.json"))  # fmt: skip
    return json.loads(paths[0].read_text())["read_processing"]


def meta(sample: str) -> dict:
    return json.loads((FIXTURE / "salmon" / sample / "aux_info" / "meta_info.json").read_text())


def test_fates_are_the_reports_own_and_add_up(ran):
    """Trim Galore's fates and Salmon's equal their own reports; every read has one fate.
    Mutation: read 2 ranked before read 1."""
    _, table = ran
    got = {r["sample"]: r for r in table.iter_rows(named=True)}
    assert sorted(got) == sorted(SAMPLES)
    for sample in SAMPLES:
        tg, salmon, row = report(sample), meta(sample), got[sample]
        for column in fates.TRIM:
            assert row[column] == tg[column], (sample, column)
        assert row["num_processed"] == salmon["num_processed"]
        assert row["num_mapped"] == salmon["num_mapped"]
        # The test profile runs BBSplit between Trim Galore and Salmon.
        assert row["removed_before_salmon"] == tg["reads_written"] - salmon["num_processed"] >= 0
        assert row["not_assigned"] == row["num_processed"] - row["num_mapped"] - row[
            "num_decoy_fragments"] >= 0  # fmt: skip


def test_values_take_read_1_and_are_unscaled(ran):
    """General Statistics scales counts to millions; the read counts come from the module
    tables, whole. Mutation: tg_total_reads read from General Statistics (0.05, not 50,000)."""
    values, _ = ran
    got = {r["sample"]: r for r in values.iter_rows(named=True)}
    for sample in SAMPLES:
        assert got[sample]["total_reads"] == report(sample)["total_reads"]
        assert got[sample]["salmon-percent_mapped"] == pytest.approx(meta(sample)["percent_mapped"])
    raw = pl.read_csv(FIXTURE / DATA / "multiqc_fastqc_fastqc_raw.txt", separator="\t")
    first = dict(raw.select("Sample", "Total Sequences").rows())
    for sample in SAMPLES:
        name = f"{sample}_1" if sample in PAIRED else sample
        assert got[sample]["Total Sequences"] == first[name]
    assert {m.column for m in REG.metrics.values() if m.module == "fastqc"} <= set(values.columns)
    verdicts = set(values["per_base_sequence_quality"])
    assert verdicts and verdicts <= {"pass", "warn", "fail"}


def test_a_missing_salmon_row_stops_the_run(tmp_path):
    """Mutation: a sample absent from one table drawn without it."""
    outdir = copy(tmp_path)
    path = outdir / DATA / "multiqc_salmon.txt"
    path.write_text("".join(line for line in path.read_text().splitlines(keepends=True)
                            if not line.startswith("WT_REP2\t")))  # fmt: skip
    with pytest.raises(QcError, match=r"no salmon report for \['WT_REP2'\]"):
        fates.checked(multiqc_io.nfcore(outdir), SAMPLES)


def test_a_missing_raw_fastqc_row_stops_the_run(tmp_path):
    """Mutation: FastQC left out of the tables checked for completeness."""
    outdir = copy(tmp_path)
    path = outdir / DATA / "multiqc_fastqc_fastqc_raw.txt"
    path.write_text("".join(line for line in path.read_text().splitlines(keepends=True)
                            if not line.startswith("RAP1_UNINDUCED_REP2\t")))  # fmt: skip
    with pytest.raises(QcError, match=r"no fastqc_fastqc_raw report for \[.RAP1_UNINDUCED_REP2"):
        fates.checked(multiqc_io.nfcore(outdir), SAMPLES)


def test_reads_without_a_fate_stop_the_run(tmp_path):
    """Mutation: the fate checks dropped; a Trim Galore total that misses a fate, and Salmon
    processing more than was written, are each named."""
    outdir = copy(tmp_path)
    path = outdir / DATA / "multiqc_trim_galore.txt"
    frame = pl.read_csv(path, separator="\t", infer_schema=False)
    frame = frame.with_columns(
        pl.when(pl.col("Sample") == "WT_REP1_1")
        .then(pl.lit("10"))
        .otherwise(pl.col("reads_written"))
        .alias("reads_written")
    )
    frame.write_csv(path, separator="\t")
    run = fates.checked(multiqc_io.nfcore(outdir), SAMPLES)
    said = r"miss its total for \['WT_REP1'\].*written for \['WT_REP1'\]"
    with pytest.raises(QcError, match=said):
        fates.fates(run, SAMPLES)
