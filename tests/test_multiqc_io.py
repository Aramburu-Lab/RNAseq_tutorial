"""ctx.qc.multiqc_io on the rnaseq fixture: MultiQC's tables against the reports themselves."""

import json
import shutil
from pathlib import Path

import polars as pl
import pytest

from ctx.io.naming import lint
from ctx.qc import QcError, multiqc_io

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "rnaseq"
SAMPLES = ("RAP1_IAA_30M_REP1", "RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2", "WT_REP1", "WT_REP2")
FATES = ("reads_written", "reads_too_short", "reads_too_long", "reads_too_many_n",
         "reads_discarded_untrimmed")  # fmt: skip


def aliases() -> dict[str, str]:
    """Every name the fixture's samples appear under in MultiQC's tables."""
    names = {}
    for sample in SAMPLES:
        for suffix in ("", "_trimmed", "_trimmed_1", "_trimmed_2"):
            names[sample + suffix] = sample
    return names


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    inputs = tmp_path_factory.mktemp("inputs")
    for part in ("trimgalore", "salmon"):
        shutil.copytree(FIXTURE / part, inputs / part)
    out = tmp_path_factory.mktemp("outputs") / "02_qc" / "multiqc"
    return multiqc_io.run([inputs], out, modules=("trim_galore", "salmon")), out


def unit(name: str) -> tuple[str, str]:
    """(sample, mate) of a report or table name: MultiQC drops `_trimmed` from single-end ones."""
    return aliases()[name], "2" if name.endswith("_trimmed_2") else "1"


def reports() -> dict[tuple[str, str], dict]:
    """The Trim Galore JSON reports, read here without MultiQC."""
    found = {}
    for path in (FIXTURE / "trimgalore").glob("*_trimming_report.json"):
        stem = path.name.removesuffix(".fastq.gz_trimming_report.json")
        found[unit(stem)] = json.loads(path.read_text())["read_processing"]
    return found


def test_trim_galore_table_equals_the_reports(ran):
    run, _ = ran
    table = multiqc_io.numeric(multiqc_io.table(run, "trim_galore"), ("total_reads", *FATES))
    by_name = {unit(r["Sample"]): r for r in table.iter_rows(named=True)}
    source = reports()
    assert set(by_name) == set(source)
    for name, values in source.items():
        for column in ("total_reads", *FATES):
            assert by_name[name][column] == values[column], (name, column)


def test_every_read_has_a_fate(ran):
    """total = written + too short + too long + too many N + discarded, exactly, per report."""
    run, _ = ran
    table = multiqc_io.numeric(multiqc_io.table(run, "trim_galore"), ("total_reads", *FATES))
    gap = table.select(pl.col("total_reads") - pl.sum_horizontal(FATES))
    assert gap.to_series().to_list() == [0.0] * table.height


def test_reads_written_reach_salmon(ran):
    """What Trim Galore wrote is what Salmon processed: read 1 for paired samples (D-022)."""
    run, _ = ran
    trim = multiqc_io.table(run, "trim_galore")
    names = multiqc_io.map_names(trim["Sample"], aliases())
    written = {}
    for row in trim.iter_rows(named=True):
        if not row["Sample"].endswith("_trimmed_2"):
            written[names[row["Sample"]]] = float(row["reads_written"])
    for sample in SAMPLES:
        meta = json.loads((FIXTURE / "salmon" / sample / "aux_info" / "meta_info.json").read_text())
        assert written[sample] == meta["num_processed"], sample
    salmon = multiqc_io.numeric(multiqc_io.table(run, "salmon"), ["num_processed"])
    assert dict(salmon.select("Sample", "num_processed").rows()) == {
        s: float(json.loads((FIXTURE / "salmon" / s / "aux_info" / "meta_info.json").read_text())
                 ["num_processed"]) for s in SAMPLES
    }  # fmt: skip


def test_unmapped_names_are_listed(ran):
    with pytest.raises(QcError, match=r"\['STRANGER', 'ZED_trimmed_1'\]"):
        multiqc_io.map_names(["WT_REP1", "STRANGER", "ZED_trimmed_1"], aliases())


def test_a_missing_report_stops_the_run(ran, tmp_path):
    run, _ = ran
    multiqc_io.require_sources(run, "trim_galore", SAMPLES, aliases())
    copy = tmp_path / "multiqc_data"
    shutil.copytree(run.data, copy)
    sources = (copy / "multiqc_sources.txt").read_text().splitlines()
    kept = [line for line in sources if "WT_REP2_trimmed" not in line]
    (copy / "multiqc_sources.txt").write_text("\n".join(kept) + "\n")
    with pytest.raises(QcError, match=r"no trim_galore report for \['WT_REP2'\]"):
        multiqc_io.require_sources(multiqc_io.MultiqcRun(copy, run.version), "trim_galore",
                                   SAMPLES, aliases())  # fmt: skip


def test_multiqc_folder_is_exempt_from_the_naming_lint(ran):
    _, out = ran
    scanned, violations = lint(out.parent)
    assert violations == []
    assert any(out.rglob("multiqc_trim_galore.txt"))
    (out.parent / "badly named.csv").write_text("x")
    assert lint(out.parent)[1] == ["badly named.csv"]
    assert scanned == 0
