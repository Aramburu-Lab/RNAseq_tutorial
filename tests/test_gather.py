"""ctx.gather.fetchngs.read on the real fixtures, and on tmp copies with one fault injected."""

import csv
import re
import shutil
from collections import Counter
from pathlib import Path

import polars as pl
import pytest

from ctx.gather.curation import library_replicates, rank_replicates
from ctx.gather.fetchngs import FETCHNGS_VERSION, VERSIONS_FILE, FetchngsError, read
from ctx.gather.harmonise import sheet

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "fetchngs"


def copy(tmp_path, mode="download"):
    target = tmp_path / mode
    shutil.copytree(FIXTURES / mode, target)
    return target


def rows(outdir, name="samplesheet"):
    with (outdir / "samplesheet" / f"{name}.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def rewrite(outdir, name, table):
    with (outdir / "samplesheet" / f"{name}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(table)


def test_reads_download_fixture():
    outdir = FIXTURES / "download"
    result = read(outdir)
    with_md5 = {p.name.removesuffix(".md5") for p in (outdir / "fastq" / "md5").glob("*.md5")}
    statuses = Counter(result.runs["status_1"])
    assert statuses == {"verified": 2, "no_md5": 2}
    assert {Path(p).name for p in result.runs["fastq_1"] if p} >= {
        n for n in with_md5 if "_2." not in n
    }
    assert result.unlisted_fastq == ()
    assert result.provenance["fetchngs"] == FETCHNGS_VERSION
    assert re.fullmatch(r"[0-9a-f]{7}", result.provenance["fetchngs_commit"])


def test_original_paths_unreachable():
    written = [r[k] for r in rows(FIXTURES / "download") for k in ("fastq_1", "fastq_2") if r[k]]
    assert written
    # The sheet's own paths must point nowhere, or reading by file name would prove nothing.
    assert not any(Path(p).exists() or (ROOT / p).exists() for p in written)
    assert set(read(FIXTURES / "download").runs["status_1"]) == {"verified", "no_md5"}


def test_md5_mismatch_fails(tmp_path):
    outdir = copy(tmp_path)
    target = outdir / "fastq" / "SRX9626017_SRR13191702_1.fastq.gz"
    data = bytearray(target.read_bytes())
    data[-1] ^= 1
    target.write_bytes(bytes(data))
    with pytest.raises(FetchngsError, match=r"SRX9626017_SRR13191702_1\.fastq\.gz is mismatch"):
        read(outdir)


def test_missing_fastq_fails(tmp_path):
    outdir = copy(tmp_path)
    (outdir / "fastq" / "ERX1234253_ERR1160846.fastq.gz").unlink()
    with pytest.raises(FetchngsError, match=r"ERX1234253_ERR1160846\.fastq\.gz is missing"):
        read(outdir)


def test_missing_metadata_fails(tmp_path):
    outdir = copy(tmp_path)
    (outdir / "metadata" / "ERR1160846.runinfo_ftp.tsv").unlink()
    with pytest.raises(FetchngsError, match="input ERR1160846: no metadata"):
        read(outdir)


def test_missing_input_list_fails(tmp_path):
    outdir = copy(tmp_path)
    (outdir / "ids.csv").unlink()
    with pytest.raises(FetchngsError, match="input list ids.csv is not in the outdir"):
        read(outdir)


def test_unlisted_fastq_reported(tmp_path):
    outdir = copy(tmp_path)
    (outdir / "fastq" / "SRX0_SRR0_3.fastq.gz").write_bytes(b"x")
    assert read(outdir).unlisted_fastq == ("SRX0_SRR0_3.fastq.gz",)


def test_two_params_refused(tmp_path):
    outdir = copy(tmp_path)
    (params,) = (outdir / "pipeline_info").glob("params_*.json")
    shutil.copyfile(params, params.with_name("params_2099-01-01_00-00-00.json"))
    with pytest.raises(FetchngsError, match="2 params_"):
        read(outdir)


def test_missing_strandedness_fails(tmp_path):
    outdir = copy(tmp_path)
    table = [{k: v for k, v in r.items() if k != "strandedness"} for r in rows(outdir)]
    rewrite(outdir, "samplesheet", table)
    with pytest.raises(FetchngsError, match="--nf_core_pipeline rnaseq"):
        read(outdir)


@pytest.mark.parametrize("change", ["other_release", "no_versions_file"])
def test_other_version_refused(tmp_path, change):
    outdir = copy(tmp_path)
    versions = outdir / "pipeline_info" / VERSIONS_FILE
    if change == "no_versions_file":
        versions.unlink()
    else:
        versions.write_text(versions.read_text().replace(f"v{FETCHNGS_VERSION}", "v1.13.0"))
    with pytest.raises(FetchngsError, match="missing|1.13.0"):
        read(outdir)


def test_metadata_only():
    outdir = FIXTURES / "metadata_only"
    runs = read(outdir).runs
    assert set(runs["status_1"]) == {"not_downloaded"}
    paired = runs.filter(runs["library_layout"] == "PAIRED")
    assert set(paired["status_2"]) == {"not_downloaded"}
    assert set(runs.filter(runs["library_layout"] == "SINGLE")["status_2"]) == {""}
    assert set(runs["fastq_1"]) == {""}


def test_multi_run():
    outdir = FIXTURES / "metadata_only"
    counted = Counter(r["sample"] for r in rows(outdir))
    multi = {s: k for s, k in counted.items() if k > 1}
    assert multi == {"SRX9504942": 4}
    runs = read(outdir).runs.filter(read(outdir).runs["sample"] == "SRX9504942").sort("run_index")
    assert runs["n_runs"].to_list() == [4] * 4
    assert runs["run_index"].to_list() == [1, 2, 3, 4]
    numbers = [int(re.search(r"\d+$", r).group()) for r in runs["run_accession"]]
    assert numbers == sorted(numbers)


def test_idmap_disagreement_fails(tmp_path):
    outdir = copy(tmp_path)
    table = rows(outdir, "id_mappings")
    table[0]["sample_title"] = "edited by hand"
    rewrite(outdir, "id_mappings", table)
    with pytest.raises(FetchngsError, match=r"id_mappings.csv disagrees on \['sample_title'\]"):
        read(outdir)


def test_every_problem_listed_together(tmp_path):
    outdir = copy(tmp_path)
    (outdir / "metadata" / "ERR1160846.runinfo_ftp.tsv").unlink()
    (outdir / "fastq" / "SRX11047067_SRR14709033.fastq.gz").unlink()
    with pytest.raises(FetchngsError) as raised:
        read(outdir)
    assert "input ERR1160846" in str(raised.value)
    assert "SRX11047067_SRR14709033.fastq.gz is missing" in str(raised.value)


def test_preconditions_stop_first(tmp_path):
    # Order: a wrong release is reported alone, before files it would describe differently.
    outdir = copy(tmp_path)
    versions = outdir / "pipeline_info" / VERSIONS_FILE
    versions.write_text(versions.read_text().replace(f"v{FETCHNGS_VERSION}", "v1.13.0"))
    (outdir / "fastq" / "ERX1234253_ERR1160846.fastq.gz").unlink()
    with pytest.raises(FetchngsError) as raised:
        read(outdir)
    assert "1.13.0" in str(raised.value)
    assert "missing" not in str(raised.value)


def test_counts_keep_int64_when_they_do_not_fit():
    checked = read(FIXTURES / "download")
    deep = checked.runs.with_columns(
        base_count=pl.when(pl.col("run_accession") == "ERR1160846")
        .then(pl.lit("23415996610"))
        .otherwise(pl.col("base_count"))
    )
    frame = sheet(type(checked)(deep, checked.unlisted_fastq, checked.provenance))
    assert frame.schema["base_count"] == pl.Int64
    assert frame.schema["read_count"] == pl.Int32
    assert 23415996610 in frame["base_count"].to_list()
    assert sheet(checked).schema["base_count"] == pl.Int32


def test_runs_of_one_experiment_are_one_replicate():
    """An ENA experiment is one library, so a fetchngs experiment's runs share one replicate
    number, and a second experiment of the same group is the next (D-063). The fixture's one
    experiment of four runs joins the group of an experiment whose run number is lower. Catches
    each run numbered as a replicate (2, 3, 4, 5), and the experiments ranked by a later run."""
    runs = read(FIXTURES / "metadata_only").runs.select("run_accession", "experiment_accession")
    (several,) = [e for e, n in runs["experiment_accession"].value_counts().rows() if n > 1]
    number = pl.col("run_accession").str.extract(r"(\d+)$").cast(pl.Int64)
    lowest = runs.filter(pl.col("experiment_accession") != several).sort(number)["run_accession"][0]
    grouped = (pl.col("experiment_accession") == several) | (pl.col("run_accession") == lowest)
    frame = runs.with_columns(
        donor=pl.when(grouped).then(pl.lit("d")).otherwise(pl.col("run_accession")),
        cell_type=pl.lit("t"),
        condition=pl.lit("c"),
    )
    numbered = frame.with_columns(replicate=library_replicates(frame, {}))
    expected = pl.when(pl.col("experiment_accession") == several).then(2).otherwise(1).alias("e")
    assert numbered["replicate"].to_list() == numbered.select(expected)["e"].to_list()
    first = runs.filter(pl.col("experiment_accession") == several).select(number.min()).item()
    assert runs.filter(pl.col("run_accession") == lowest).select(number).item() < first


def test_library_ranked_by_its_first_run_then_by_name():
    """A library takes its place in the group by its first run, and a tie on the run number
    is decided by the name, whatever the row order. Library x's runs 1 and 5 bracket y's run 3;
    the rows list z_7 before its twin w_7. Catches a library ranked by its last run, and a tie
    left to the row order."""
    frame = pl.DataFrame({
        "run_accession": ["x_5", "y_3", "x_1", "z_7", "w_7"],
        "experiment_accession": ["x", "y", "x", "z", "w"],
        "donor": ["d", "d", "d", "e", "e"],
        "cell_type": ["t"] * 5,
        "condition": ["c"] * 5,
    })  # fmt: skip
    assert library_replicates(frame, {}).to_list() == [1, 2, 1, 2, 1]
    # M5 ranks runs with rank_replicates itself, so the tie is decided there too.
    twins = frame.filter(pl.col("donor") == "e")
    assert rank_replicates(twins, ("donor",), {}).to_list() == [2, 1]
