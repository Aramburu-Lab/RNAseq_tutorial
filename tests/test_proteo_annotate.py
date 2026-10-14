"""The declared annotation join of ctx gather --gdc-dir (ctx.proteo.annotate) on a made-up
freeze; each test names the mutation it catches. Nothing here comes from a real dataset."""

import io
import tarfile
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from gdc_fixture import SPEC
from proteo_root import gct, resource_toml

from ctx.dataset import DatasetError, load
from ctx.proteo import ProteoError
from ctx.proteo.annotate import sheet_columns

# Three tumours of three cases, one unchosen copy, and a normal of the first case.
SHEET = pl.DataFrame({
    "sample": ["A1", "A2", "A3", "A1b", "N1"],
    "case": ["C1", "C2", "C3", "C1", "C1"],
    "sample_id": ["S1", "S2", "S3", "S1b", "SN1"],
    "sample_type": ["Primary Tumor"] * 4 + ["Solid Tissue Normal"],
    "chosen": [True, True, True, False, True],
})  # fmt: skip
PROTEO = """
[proteo]
design = "across_samples"
protein = "freeze"
member = "p.gct"
sample_prefix = "X"
plex = "Plex"
"""
COLUMNS = """
[columns.er_status]
source = "ER"
meaning = "made-up ER call"
values = {{ positive = "ER+", negative = "ER-"{extra} }}
"""


def root_with(tmp_path: Path, er: list[str], *, extra: str = "", source: str = "ER") -> tuple:
    """A root holding a freeze of columns XC1, S2 (by sample ID) and XC9, and its dataset."""
    body = gct(["G1"], ["XC1", "S2", "XC9"], np.zeros((1, 3)),
               {"Plex": ["1", "1", "2"], "ER": er})  # fmt: skip
    (tmp_path / "resources").mkdir(parents=True)
    with tarfile.open(tmp_path / "resources" / "freeze.tar.gz", "w:gz") as tar:
        info = tarfile.TarInfo("freeze/p.gct")
        info.size = len(body)
        tar.addfile(info, io.BytesIO(body))
    spec = tmp_path / "fake_brca.toml"
    spec.write_text(SPEC.format(protocol="unstranded", cdr="")
                    + resource_toml("freeze", "resources/freeze.tar.gz") + PROTEO
                    + COLUMNS.format(extra=extra).replace('"ER"', f'"{source}"'))  # fmt: skip
    return tmp_path, load(spec)


def test_values_land_on_the_paired_chosen_tumours_only(tmp_path):
    """Mutation caught: a join by case alone (S2's column, named by sample ID, left out); the
    values left unmapped (``positive`` in the sheet, where contrasts name ER+); every file of
    a case labelled (the unchosen copy and the normal of C1 taking the tumour's call)."""
    root, spec = root_with(tmp_path, ["negative", "positive", "positive"])
    out = sheet_columns(SHEET, spec, root)
    assert out["er_status"].to_list() == ["ER-", "ER+", None, None, None]
    assert out.select(SHEET.columns).equals(SHEET)


def test_an_unmapped_value_stops_gather(tmp_path):
    """Mutation caught: an unmapped value passed through or made null in silence, which drops
    its tumours from every contrast unseen."""
    root, spec = root_with(tmp_path, ["negative", "unknown", "positive"])
    with pytest.raises(ProteoError, match="'unknown' is not in values"):
        sheet_columns(SHEET, spec, root)
    root2, spec2 = root_with(tmp_path / "b", ["negative", "unknown", "positive"],
                             extra=', unknown = "unknown"')  # fmt: skip
    assert sheet_columns(SHEET, spec2, root2)["er_status"][1] == "unknown"


def test_a_missing_source_or_design_is_refused(tmp_path):
    """Mutation caught: a source the freeze lacks read as an all-null column; columns accepted
    on a dataset whose protein table has no sample annotation."""
    root, spec = root_with(tmp_path, ["negative", "positive", "positive"], source="ER.Status")
    with pytest.raises(ProteoError, match="no annotation 'ER.Status'"):
        sheet_columns(SHEET, spec, root)
    bare = tmp_path / "bare.toml"
    bare.write_text(SPEC.format(protocol="unstranded", cdr="") + COLUMNS.format(extra=""))
    with pytest.raises(DatasetError, match="across-samples"):
        sheet_columns(SHEET, load(bare), root)
