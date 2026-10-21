"""F.1: M10's figures, each redrawn from its spec's inputs alone in a clean folder.

Every figure with a spec writes M10's data table, PDF and SVG byte for byte, and the score's
curves carry the module's live panel. Each test names its mutation.
"""

from pathlib import Path

import polars as pl
import pytest
from clinical_root import ctx, make
from conftest import ctx_run, once
from test_clinical_page import controls
from test_figures import copied, drawn_alike, redrawn_each, specs

from ctx.clinical import analyses
from ctx.figures.params import ParamsError, load_params
from ctx.figures.spec import parse
from ctx.report.builder import unpacked

CLINICAL = ("12_clinical",)
# Every figure the made-up root draws: ER curves under the four endpoints, the ER box, the
# score's curves and scan, the competing risks, power and the confounding map.
EXPECTED = {"er:OS", "er:DSS", "er:DFI", "er:PFI", "expression:er", "km_score", "scan_score",
            "competing", "power", "confounding"}  # fmt: skip
# The options a run can change, each away from its default.
OPTIONS = ("--endpoint", "OS", "--split", "tertiles", "--horizon", "36", "--covariates", "age",
           "--cox-strata", "histology")  # fmt: skip


def _kinds(found) -> set[str]:
    """The kinds, an ER curve and the box named by their parameter."""
    names = set()
    for spec in found:
        extra = getattr(spec.params, "endpoint", None) or getattr(spec.params, "box", None)
        kind = parse(spec.figure_id)[2]
        names.add(f"{kind}:{extra}" if extra else kind)
    return names


def _drawn(root: Path) -> int:
    """How many figures M10 wrote."""
    return len(list((root / "outputs" / "12_clinical").glob("*_v1_white.pdf")))


def _redraws_alike(root: Path, found, tmp_path: Path) -> None:
    """Each spec redrawn in its own clean folder equals M10's files; the live panel too."""
    for spec, csv in zip(found, redrawn_each(found, root, tmp_path), strict=True):
        module_csv = next(root.glob(f"outputs/12_clinical/{spec.stem}.csv"))
        assert drawn_alike(csv, module_csv) == [], spec.stem
        if spec.figure_id.endswith("_km_score"):
            pages = (unpacked(p.with_suffix(".html").read_text()) for p in (module_csv, csv))
            module, again = (controls(page, "survival") for page in pages)
            assert module == again, spec.stem


def test_every_figure_redraws_byte_for_byte(clinical_run, tmp_path):
    """Every figure M10 draws has a spec, and each redraws byte for byte, the score's live
    panel equal to the module's. Catches the panel values read in another order, the build
    list taken from the wrong rows, the screen left out of the live panel, and a kind whose
    stem no longer matches."""
    root = clinical_run[0]
    found = specs(root, CLINICAL, root.parent / "fake_brca.toml")
    assert _kinds(found) == EXPECTED
    assert len(found) == _drawn(root)
    _redraws_alike(root, found, tmp_path)


@pytest.fixture(scope="session")
def optioned(tmp_path_factory):
    """ctx clinical with every figure option away from its default and its own gene list."""
    return once(tmp_path_factory, "m10_options", _optioned)


def _optioned(tmp: Path) -> tuple[Path, Path]:
    root, spec = make(tmp)
    listed = tmp / "mine.txt"
    listed.write_text("# source: x:y\n# count: 2\nG2\nG7\n")
    ctx("--root", str(root), "clinical", "--dataset", "fake_brca", "--dataset-file", str(spec),
        "--genes-from", str(listed), *OPTIONS)  # fmt: skip
    return root, spec


def test_the_runs_options_are_the_ones_redrawn(optioned, tmp_path):
    """A run given another endpoint, split, horizon, covariates, strata and gene list redraws
    byte for byte, its options table holding each. Catches an option taken from its default
    in the redraw: the median split, five years, the dataset file's endpoint and list."""
    root, spec = optioned
    found = specs(root, CLINICAL, spec)
    row = pl.read_parquet(next(root.glob("outputs/12_clinical/*_figure_options_v1.parquet")))
    assert row.row(0, named=True) == {"endpoint": "OS", "split": "tertiles", "horizon": 36.0,
                                      "covariates": "age", "strata": "histology",
                                      "source": "mine.txt"}  # fmt: skip
    chosen = [
        s
        for s in found
        if parse(s.figure_id)[2] in ("km_score", "scan_score", "power", "competing")
    ]
    assert len(chosen) == 4, _kinds(found)  # fmt: skip
    _redraws_alike(root, chosen, tmp_path)


def test_values_of_other_samples_are_refused(clinical_run, tmp_path):
    """Panel values whose sample columns are not the patient table's, in its order, stop the
    redraw, named. Catches the values matched to patients by position alone, which would
    score each patient with another's expression."""
    root = clinical_run[0]
    spec = next(s for s in specs(root, CLINICAL, root.parent / "fake_brca.toml")
                if s.figure_id.endswith("_power"))  # fmt: skip
    folder = copied(spec, root, tmp_path / "swapped")
    values = folder / next(i.file for i in spec.inputs if i.name == "panel_values")
    frame = pl.read_parquet(values)
    first, second = frame.columns[3:5]
    frame.select(*frame.columns[:3], second, first, *frame.columns[5:]).write_parquet(values)
    raw = (folder / "params.toml").read_text().splitlines()
    (folder / "params.toml").write_text(
        "\n".join(x for x in raw if not x.startswith("panel_values = ") or ".parquet" in x) + "\n"
    )
    done = ctx_run("figure", spec.figure_id, "--params", "params.toml", cwd=folder)
    assert done.returncode != 0
    assert "not the patient table's" in done.stderr, done.stderr[-2000:]


def test_cramers_v_is_the_same_every_call(clinical_run):
    """Cramér's V of every confounding pair is one value over 50 calls on one table. Catches the
    counts grouped in hashed order, which changed the chi-squared sum's last bits from call to
    call, so the confounding map's table differed between the run and its redraw."""
    root = clinical_run[0]
    patients = pl.read_parquet(next(root.glob("outputs/12_clinical/*_patients_v1.parquet")))
    table = patients.with_columns(analyses.age_group(pl.col("age")).alias("age_group"),
                                  pl.col("stage").cast(pl.String))  # fmt: skip
    names = ("er", "stage", "histology", "menopause", "age_group", "plate", "site")
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            seen = {analyses.cramers_v(table, a, b)[0] for _ in range(50)}
            assert len(seen) == 1, (a, b, seen)


def test_an_endpoint_or_box_the_file_does_not_draw_is_refused(clinical_run, tmp_path):
    """A params file naming an endpoint the CDR lacks, or a box the dataset file does not
    declare, is refused with the name. Catches either check left out, which would draw a
    figure M10 never drew or fail inside the drawer."""
    root = clinical_run[0]
    found = specs(root, CLINICAL, root.parent / "fake_brca.toml")
    for kind, old, new, said in (("er", 'endpoint = "OS"', 'endpoint = "XYZ"',
                                  "endpoint: 'XYZ' is not an endpoint"),
                                 ("expression", 'box = "er"', 'box = "stage_group"',
                                  "box: 'stage_group' is not the clinical.box")):  # fmt: skip
        spec = next(s for s in found if parse(s.figure_id)[2] == kind
                    and getattr(s.params, "endpoint", "OS") == "OS")  # fmt: skip
        folder = copied(spec, root, tmp_path / kind)
        raw = (folder / "params.toml").read_text()
        assert raw.count(old) == 1, old
        (folder / "params.toml").write_text(raw.replace(old, new))
        with pytest.raises(ParamsError, match=said):
            load_params(folder / "params.toml", spec.figure_id)
