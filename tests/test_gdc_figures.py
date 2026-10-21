"""Which figures the registry offers a GDC dataset, and `ctx figure` on a gathered GDC root
(R1.3b). Each test names the mutation it catches."""

import subprocess
import sys
from pathlib import Path

import pytest
from clinical_root import CLINICAL
from gdc_fixture import SPEC, aliquot, build

from ctx import pipeline
from ctx.dataset import load
from ctx.figures import kinds
from ctx.figures.params import ParamsError, load_params

CTX = Path(sys.executable).parent / "ctx"
DE = ('\n[[contrasts]]\nname = "tumour_normal"\nvariable = "condition"\nreference = "normal"\n'
      'tested = "tumour"\ndesign = "~donor + group"\npaired = true\nn = [2, 2]\n'
      '\n[de]\nfeatured = "tumour_normal"\n')  # fmt: skip
# Measured 2026-10-21 by drawing each with `ctx figure` from a gathered made-up GDC root
# after de and clinical ran: every one exits 0.
OFFERED = {
    "de": ["dispersion", "ma", "phist", "volcano"],
    "clinical": ["competing", "confounding", "er", "expression", "km_score", "power",
                 "scan_score"],
}  # fmt: skip


def gdc_spec(tmp: Path) -> Path:
    """A GDC dataset file with contrasts, [de] and [clinical], and no [function]."""
    path = tmp / "fake_brca.toml"
    path.write_text(SPEC.format(protocol="unstranded", cdr='attribution = "tcga"') + CLINICAL + DE)
    return path


def test_the_registry_offers_a_gdc_dataset_only_what_its_modules_draw(tmp_path):
    """For each module the GDC route runs, the kinds made are the listed ones. Catches M1's
    composition and depth bars made for every dataset, which `ctx figure` drew from a GDC sheet
    into a ColumnNotFoundError, and M6's kinds made from [de] alone, where `ctx function`, which
    needs [function], never runs and the background lesson stopped on a KeyError."""
    spec = load(gdc_spec(tmp_path))
    options = pipeline.Options(spec.slug, spec.path, gdc_dirs=(tmp_path,))
    route = [module for module, _ in pipeline.steps(spec, options)]
    assert route == ["gather", "de", "clinical", "report"]
    made: dict[str, list[str]] = {}
    for k in kinds.KINDS:
        if k.module in route and k.parsed is None and k.made(spec):
            made.setdefault(k.module, []).append(k.figure)
    assert {m: sorted(f) for m, f in made.items()} == OFFERED
    assert not [k for k in kinds.KINDS if k.module == "function" and k.made(spec)]


@pytest.fixture(scope="module")
def gathered(tmp_path_factory) -> tuple[Path, Path]:
    """A made-up GDC project gathered into a root: the root's gather folder and the dataset file."""
    tmp = tmp_path_factory.mktemp("gdc_figures")
    build(tmp / "gdc", [aliquot(1), aliquot(2), aliquot(1, code="11"), aliquot(2, code="11")])
    spec = gdc_spec(tmp)
    done = subprocess.run([str(CTX), "--root", str(tmp / "root"), "gather", "--dataset",
                           "fake_brca", "--gdc-dir", str(tmp / "gdc"), "--dataset-file", str(spec)],
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return tmp / "root" / "outputs" / "01_gather", spec


@pytest.mark.parametrize("figure", ["composition", "library_sizes"])
def test_ctx_figure_refuses_gathers_figures_for_a_gdc_sheet(gathered, tmp_path, figure):
    """`ctx figure breast_m1_<figure>` on a GDC sheet exits 1 naming the figure, with no
    traceback. Catches the refusal left out of the params check, which reaches the drawer and
    stops on a ColumnNotFoundError for cell_type or read_count."""
    out, spec = gathered
    sheet = next(out.glob("*_fake_brca_samplesheet_v1.parquet"))
    (tmp_path / sheet.name).write_bytes(sheet.read_bytes())
    (tmp_path / spec.name).write_bytes(spec.read_bytes())
    registry = "" if figure == "composition" else 'qc_registry = "qc_metrics.toml"\n'
    if registry:
        (tmp_path / "qc_metrics.toml").write_bytes(
            (Path(__file__).resolve().parents[1] / "conf" / "qc_metrics.toml").read_bytes()
        )
    params = tmp_path / "params.toml"
    params.write_text(
        f'figure = "breast_m1_{figure}"\ndataset = "fake_brca"\n\n[inputs]\n'
        f'samplesheet = "{sheet.name}"\ndataset_file = "{spec.name}"\n{registry}'
    )
    done = subprocess.run([str(CTX), "figure", f"breast_m1_{figure}", "--params", str(params),
                           "--out", str(tmp_path / "drawn")], capture_output=True, text=True,
                          check=False)  # fmt: skip
    assert done.returncode == 1, done.stdout
    assert f"{figure}: fake_brca's dataset file makes no such figure" in done.stderr, done.stderr
    assert "Traceback" not in done.stderr and "Error:" not in done.stderr, done.stderr[-2000:]
    assert not (tmp_path / "drawn").exists()


def test_m6s_background_without_function_is_refused_by_name(tmp_path):
    """A params file for M6's background lesson of a dataset file with [de] and no [function]
    is refused before drawing. Catches M6's kinds made from [de] alone, whose inputs then omit
    the gene sets the drawer reads."""
    spec = gdc_spec(tmp_path)
    (tmp_path / "results.parquet").write_bytes(b"")
    params = tmp_path / "params.toml"
    params.write_text('figure = "breast_m6_background"\ndataset = "fake_brca"\n'
                      'contrast = "tumour_normal"\n\n[inputs]\nresults = "results.parquet"\n'
                      f'dataset_file = "{spec.name}"\n')  # fmt: skip
    with pytest.raises(ParamsError, match="background: fake_brca's dataset file makes no such"):
        load_params(params, "breast_m6_background")
