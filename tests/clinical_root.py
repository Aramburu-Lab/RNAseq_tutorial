"""A made-up GDC project gathered into a ctx root, a TCGA-CDR copy beside it, for ``ctx clinical``.

Used by tests/test_clinical_page.py and tests/test_cli.py. Nothing here comes from a real
dataset: barcodes start with FAKE (tests/gdc_fixture.py).
"""

import random
import subprocess
import sys
from datetime import date
from pathlib import Path

import polars as pl
from gdc_fixture import SPEC, aliquot, build

from ctx.dataset import ENDPOINTS

CTX = Path(sys.executable).parent / "ctx"
PATIENTS = 80
CLINICAL = '\n[clinical]\ngenes = ["G1", "G3", "G5"]\nendpoint = "PFI"\n'
STAGES = ("Stage I", "Stage IIA", "Stage IIB", "Stage IIIA", "Stage IV", "Stage X")


def ctx(*args: str) -> subprocess.CompletedProcess[str]:
    """Run ctx; refuse a failed run with its stderr."""
    done = subprocess.run([str(CTX), *args], capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr[-3000:]
    return done


def cdr_rows(cases: list[str], seed: int) -> pl.DataFrame:
    """A made-up TCGA-CDR copy, as text, its endpoints drawn per patient."""
    rng = random.Random(seed)
    rows = []
    for case in cases:
        row = {"bcr_patient_barcode": case, "type": "BRCA", "gender": "FEMALE",
               "age_at_initial_pathologic_diagnosis": str(rng.randint(35, 85)),
               "ajcc_pathologic_tumor_stage": rng.choice(STAGES),
               "histological_type": rng.choice(("Infiltrating Ductal Carcinoma",
                                                 "Infiltrating Lobular Carcinoma")),
               "menopause_status": rng.choice(("Pre (<6 months since LMP)",
                                               "Post (prior bilateral ovariectomy)")),
               "Redaction": ""}  # fmt: skip
        for e in ENDPOINTS:
            row[e] = str(int(rng.random() < 0.35))
            row[f"{e}.time"] = str(rng.randint(30, 3000))
        rows.append(row)
    return pl.DataFrame(rows)


def make(tmp: Path) -> tuple[Path, Path]:
    """Gather a made-up project of ``PATIENTS`` patients into ``tmp/root`` and write its CDR.

    Returns:
        The root and the dataset file.
    """
    rng = random.Random(4)
    items = [aliquot(i) for i in range(1, PATIENTS + 1)] + [aliquot(i, code="11")
                                                           for i in range(1, 6)]  # fmt: skip
    results = ("Positive", "Positive", "Negative")
    cases = {aliquot(i).case: {"tests": [{"gene_symbol": "ESR1", "molecular_analysis_method": "IHC",
                                          "test_result": rng.choice(results)}]}
             for i in range(1, PATIENTS + 1)}  # fmt: skip
    build(tmp / "gdc", items, cases=cases, seed=9)
    spec = tmp / "fake_brca.toml"
    spec.write_text(SPEC.format(protocol="unstranded", cdr="") + CLINICAL)
    root = tmp / "root"
    ctx("--root", str(root), "gather", "--dataset", "fake_brca", "--gdc-dir", str(tmp / "gdc"),
        "--dataset-file", str(spec))  # fmt: skip
    gather = root / "outputs" / "01_gather"
    name = f"{date.today().isoformat()}_gather_fake_brca_cdr_v1.parquet"
    cdr_rows([aliquot(i).case for i in range(1, PATIENTS + 1)], 5).write_parquet(gather / name)
    return root, spec
