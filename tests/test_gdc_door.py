"""Every reader of gather's sheet goes through `gdc.chosen`; quant and process refuse GDC (S.6).

Each test names the mutation it catches.
"""

import importlib.util
import sys
from datetime import date
from pathlib import Path
from types import SimpleNamespace

import polars as pl
import pytest
from gdc_fixture import aliquot, build, spec_file
from synth import DE_GROUPS, design_counts, write_gather_outputs

from ctx import pipeline
from ctx.dataset import load as load_dataset
from ctx.gather import gdc, gdc_check
from ctx.process import ProcessError
from ctx.process import run as process_run
from ctx.quant import QuantError
from ctx.quant import run as quant_run

ROOT = Path(__file__).resolve().parents[1]
DATASET = "desynth"
UNCHOSEN = "S900"


def load(name: str):
    """A preview script as a module, registered so its dataclasses resolve their annotations."""
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / "python" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def gathered(root, *, gdc_sheet: bool) -> pl.DataFrame:
    """A gather root; as a GDC sheet, every file chosen plus an unchosen twin of S001."""
    write_gather_outputs(root, DATASET, design_counts(DE_GROUPS), date(2026, 9, 30))
    path = next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))
    sheet = pl.read_parquet(path)
    if gdc_sheet:
        twin = sheet.filter(pl.col("sample") == "S001").with_columns(sample=pl.lit(UNCHOSEN))
        sheet = pl.concat([sheet.with_columns(chosen=pl.lit(value=True)),
                           twin.with_columns(chosen=pl.lit(value=False))])  # fmt: skip
        sheet = sheet.with_columns(file_id=pl.col("sample"))
        sheet.write_parquet(path)
    return sheet


def test_a_fetchngs_sheet_passes_the_door_unchanged():
    """Mutation caught: the door filtering on a column a fetchngs sheet does not have."""
    sheet = pl.DataFrame({"sample": ["a", "b"], "run_accession": ["r1", "r2"]})
    assert gdc.chosen(sheet).equals(sheet)
    assert not gdc.is_sheet(sheet)


@pytest.mark.parametrize("column", ["drop_reason", "file_id", "md5sum", "case_id"])
def test_a_gdc_sheet_without_chosen_is_refused(column):
    """A lost `chosen` column must not pass every file as a fetchngs sheet would. Mutations
    caught: the door returning any sheet without `chosen` unchanged; a marker dropped."""
    sheet = pl.DataFrame({"sample": ["a", "b"], column: ["x", "y"]})
    with pytest.raises(gdc.GdcError, match="without its chosen column"):
        gdc.chosen(sheet)


def test_a_sample_type_is_refused_on_a_fetchngs_sheet():
    """Mutation caught: the type ignored, every run returned as if of that type."""
    with pytest.raises(gdc.GdcError, match="no GDC sample type"):
        gdc.chosen(pl.DataFrame({"sample": ["a"]}), "Primary Tumor")


@pytest.mark.parametrize("gdc_sheet", [True, False])
def test_quant_refuses_a_gdc_sheet(tmp_path, gdc_sheet):
    """Mutation caught: quant reading a GDC sheet as runs (the refusal dropped)."""
    gathered(tmp_path, gdc_sheet=gdc_sheet)
    path, sheet, problems = quant_run._sheet(tmp_path, DATASET)
    assert path is not None
    if gdc_sheet:
        assert sheet is None and len(problems) == 1
        assert f"{path.name} is a GDC sheet: ctx quant reads" in problems[0]
    else:
        assert problems == [] and sheet is not None and UNCHOSEN not in sheet["sample"]


def test_process_refuses_a_gdc_sheet(tmp_path):
    """Mutation caught: process reading a GDC sheet's files as runs (the refusal dropped)."""
    gathered(tmp_path, gdc_sheet=True)
    found, problems = process_run._sheet(tmp_path, DATASET, SimpleNamespace(fetchngs_dir=None))
    assert found is None and len(problems) == 1
    assert "is a GDC sheet: ctx process reads a fetchngs dataset's runs" in problems[0]


@pytest.mark.parametrize("script", ["preview_genecard", "preview_proteo"])
def test_the_previews_count_the_chosen_files_alone(tmp_path, script):
    """S.6: both previews counted every file of a GDC sheet. Mutation caught: either read
    bypassing `gdc.chosen` (the unchosen twin counted)."""
    sheet = gathered(tmp_path, gdc_sheet=True)
    module = load(script)
    if script == "preview_genecard":
        path = next((tmp_path / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))
        samples = module.samples(path)
    else:
        samples = module.design(tmp_path).samples
    kept = sheet.filter(pl.col("chosen"))["sample"].unique()
    assert UNCHOSEN not in samples["sample"]
    assert samples.height == kept.len()


# Cases 1 and 2 paired; case 3's normal is FFPE, which `chosen` drops, so case 3 is unpaired.
PAIRED = (aliquot(1), aliquot(1, "11"), aliquot(2), aliquot(2, "11"), aliquot(3),
          aliquot(3, "11", preservation="FFPE"))  # fmt: skip
CONTRAST = """
[[contrasts]]
name = "tumour_normal"
variable = "condition"
reference = "normal"
tested = "tumour"
design = "~donor + group"
paired = true
n = [2, 2]
"""


@pytest.mark.parametrize("contrasts", [True, False])
def test_gdc_check_counts_the_chosen_files_alone(tmp_path, contrasts):
    """R1.3: `ctx run --check`'s n is the pairs of chosen files, [2, 2] here, and [3, 3] if the
    FFPE normal is counted. Mutation caught: either `gdc.chosen()` in gdc_check dropped, the
    declared contrast's n or the n printed into the template's lines for a file with none."""
    build(tmp_path / "gdc", PAIRED)
    path = spec_file(tmp_path / "fake_brca.toml")
    if contrasts:
        path.write_text(path.read_text() + CONTRAST)
    spec = load_dataset(path)
    options = pipeline.Options(spec.slug, path, gdc_dirs=(tmp_path / "gdc",))
    problems = gdc_check.check(spec, options, tmp_path / "root")
    named = [p for p in problems if "[[contrasts]]" in p.text]
    if contrasts:
        assert named == [], [p.text for p in named]
    else:
        (missing,) = named
        assert missing.text.startswith("no [[contrasts]]") and "n = [2, 2]" in missing.lines


@pytest.mark.parametrize(("script", "error"), [("preview_quant", QuantError),
                                               ("preview_process", ProcessError)])  # fmt: skip
def test_the_quant_and_process_previews_refuse_a_gdc_sheet(tmp_path, script, error):
    """Both previews read gather's sheet past `gdc.chosen`, every file a run, where their modules
    refuse a GDC sheet. Mutation caught: either preview's sheet read bypassing the module's door."""
    gathered(tmp_path, gdc_sheet=True)
    with pytest.raises(error, match="is a GDC sheet"):
        load(script).page(tmp_path, DATASET, date(2026, 10, 21))
