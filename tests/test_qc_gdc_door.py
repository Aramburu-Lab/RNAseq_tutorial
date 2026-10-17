"""`ctx qc` reads a GDC sheet through `gdc.chosen`, the door M5 reads it through (TODO S.1)."""

from datetime import date

import polars as pl
from test_qc_run import SPEC, ctx, plant

from ctx.qc.run import read_gathered

UNCHOSEN = "S900"


def unchosen_file(root):
    """The planted sheet as a GDC sheet: every file chosen, plus a second file of S001's case.

    The second file is not chosen; its counts are S001's, so only the door tells them apart.
    """
    gathered = root / "outputs" / "01_gather"
    counts_path = next(gathered.glob("*_counts_v1.parquet"))
    sheet_path = next(gathered.glob("*_samplesheet_v1.parquet"))
    counts, sheet = pl.read_parquet(counts_path), pl.read_parquet(sheet_path)
    twin = sheet.filter(pl.col("sample") == "S001").with_columns(
        sample=pl.lit(UNCHOSEN), run_accession=pl.lit("R900")
    )
    sheet = pl.concat([sheet.with_columns(chosen=pl.lit(value=True)),
                       twin.with_columns(chosen=pl.lit(value=False))])  # fmt: skip
    counts.with_columns(pl.col("S001").alias(UNCHOSEN)).write_parquet(counts_path)
    sheet.write_parquet(sheet_path)
    return {"samplesheet": sheet_path, "counts": counts_path}, sheet


def test_the_reader_keeps_the_chosen_files_and_their_counts_alone(tmp_path):
    plant(tmp_path)
    inputs, sheet = unchosen_file(tmp_path)
    gathered = read_gathered(inputs)
    kept = sheet.filter(pl.col("chosen"))["sample"].to_list()
    assert gathered["samplesheet"]["sample"].to_list() == kept
    assert gathered["counts"].columns == ["gene_id", *kept]


def test_qc_counts_the_samples_m5_compares(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    plant(root)
    _, sheet = unchosen_file(root)
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "02_qc"
    day = date.today().isoformat()
    metrics = pl.read_csv(out / f"{day}_qc_desynth_metrics_v1.csv")
    status = pl.read_csv(out / f"{day}_qc_desynth_status_v1.csv")
    tests = pl.read_csv(out / f"{day}_qc_desynth_confounding_v1.csv")
    assert UNCHOSEN not in metrics["sample"].to_list()
    assert UNCHOSEN not in status["sample"].to_list()
    assert metrics.height == sheet.filter(pl.col("chosen")).height
    # S001 is a short type_a control: a_case compares 6 controls, never 7.
    assert set(tests.filter(pl.col("contrasts").str.contains("a_case"))["n_reference"]) == {6}
