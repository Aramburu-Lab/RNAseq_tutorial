"""ctx.qc.fates on the rnaseq fixture: values and fates against the tools' own reports."""

import json
import shutil
from pathlib import Path

import pytest

from ctx.qc import QcError, fates, registry

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "rnaseq"
SAMPLES = ["RAP1_IAA_30M_REP1", "RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2", "WT_REP1", "WT_REP2"]
REG = registry.load(ROOT / "conf" / "qc_metrics.toml")


def config(tmp: Path) -> Path:
    path = tmp / "config.json"
    path.write_text(json.dumps(registry.multiqc_config(REG, "test")))
    return path


def copy(tmp: Path) -> Path:
    out = tmp / "rnaseq"
    for part in ("trimgalore", "salmon"):
        shutil.copytree(FIXTURE / part, out / part)
    return out


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("fates")
    run = fates.run_multiqc(FIXTURE, tmp / "multiqc", config(tmp), SAMPLES)
    return fates.values(run, SAMPLES, REG), fates.fates(run, SAMPLES)


def report(sample: str) -> dict:
    """The sample's read-1 Trim Galore report, or its only one."""
    folder = FIXTURE / "trimgalore"
    paths = sorted(folder.glob(f"{sample}_trimmed_1.*.json")) or sorted(
        folder.glob(f"{sample}_trimmed.*.json"))  # fmt: skip
    return json.loads(paths[0].read_text())["read_processing"]


def meta(sample: str) -> dict:
    return json.loads((FIXTURE / "salmon" / sample / "aux_info" / "meta_info.json").read_text())


def test_fates_are_the_reports_own_and_add_up(ran):
    _, table = ran
    got = {r["sample"]: r for r in table.iter_rows(named=True)}
    assert sorted(got) == sorted(SAMPLES)
    for sample in SAMPLES:
        tg, salmon, row = report(sample), meta(sample), got[sample]
        for column in fates.TRIM:
            assert row[column] == tg[column], (sample, column)
        assert row["num_processed"] == salmon["num_processed"] == tg["reads_written"]
        assert row["num_mapped"] == salmon["num_mapped"]
        assert row["num_decoy_fragments"] == salmon["num_decoy_fragments"]
        assert row["not_assigned"] == row["num_processed"] - row["num_mapped"] - row[
            "num_decoy_fragments"] >= 0  # fmt: skip


def test_values_take_read_1_and_stay_raw(ran):
    values, _ = ran
    got = {r["sample"]: r for r in values.iter_rows(named=True)}
    for sample in SAMPLES:
        assert got[sample]["trim_galore-tg_total_reads"] == report(sample)["total_reads"]
        assert got[sample]["salmon-percent_mapped"] == pytest.approx(meta(sample)["percent_mapped"])


def test_a_missing_salmon_report_stops_the_run(tmp_path):
    outdir = copy(tmp_path)
    shutil.rmtree(outdir / "salmon" / "WT_REP2")
    with pytest.raises(QcError, match=r"no salmon report for \['WT_REP2'\]"):
        fates.run_multiqc(outdir, tmp_path / "multiqc", config(tmp_path), SAMPLES)


def test_reads_without_a_fate_stop_the_run(tmp_path):
    outdir = copy(tmp_path)
    path = next((outdir / "trimgalore").glob("WT_REP1_trimmed_1.*.json"))
    data = json.loads(path.read_text())
    data["read_processing"]["reads_written"] -= 5
    path.write_text(json.dumps(data))
    run = fates.run_multiqc(outdir, tmp_path / "multiqc", config(tmp_path), SAMPLES)
    with pytest.raises(
        QcError, match=r"miss its total for \['WT_REP1'\].*written for \['WT_REP1'\]"
    ):
        fates.fates(run, SAMPLES)
