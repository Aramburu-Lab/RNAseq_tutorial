# ruff: noqa: D103 - a proof record, run once and kept as run
"""wo-z5: every session root's pages copied to $Z5_OUT/<fixture>/, one test per fixture."""

import os
import shutil
from pathlib import Path

import pytest

FIXTURES = [
    "setup_run",
    "gather_run",
    "pipeline",
    "qc_run",
    "clinical_run",
    "proteo_run",
    "genecard_run",
    "ncrna_run",
    "quant_run",
    "process_root",
    "guide_run",
]


@pytest.mark.parametrize("name", FIXTURES)
def test_copy_pages(name, request):
    got = request.getfixturevalue(name)
    root = Path(got[0] if isinstance(got, tuple) else got)
    out = Path(os.environ["Z5_OUT"]) / name
    out.mkdir(parents=True, exist_ok=True)
    pages = sorted(root.rglob("*.html"))
    for page in pages:
        target = out / page.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(page, target)
    assert pages
