"""R1: the GDC dataset template a student copies for any TCGA project, and `ctx run --check` on it.

The template's placeholders are read from the template itself, so no project ID is typed here;
the folders are tests/gdc_fixture.py's made-up project. Each test names the mutation it catches.
"""

import importlib.util
import tomllib
import zipfile
from pathlib import Path

import gdc_mirna_fixture as mirna_fixture
import pytest
from gdc_fixture import GENES, PROJECT, RELEASE, aliquot, build, write_cdr

from ctx import pipeline
from ctx.dataset import load
from ctx.gather import gdc_check, gdc_mirna

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "conf" / "dataset_template_gdc.toml"
KIRC = ROOT / "conf" / "datasets" / "tcga_kirc.toml"
XLSX = importlib.util.find_spec("openpyxl") is not None
STAR = (aliquot(1), aliquot(1, "11"), aliquot(2), aliquot(2, "11"), aliquot(3))
CDR_NAME = "TCGA-CDR-SupplementalTableS1.xlsx"
CDR_COLUMNS = ("bcr_patient_barcode", "type", "OS", "OS.time", "PFI", "PFI.time", "DFI",
               "DFI.time", "DSS", "DSS.time")  # fmt: skip
# The made-up CDR's flags, per patient of STAR, then their Table 3 count, typed apart.
FLAGS = {"OS": (1, 0, 1), "PFI": (0, 0, 1), "DFI": (0, "#N/A", 0), "DSS": (1, 0, 0)}
TABLE3 = {"OS": (2, 1), "PFI": (1, 2), "DFI": (0, 2), "DSS": (1, 2)}
# TargetScan's family file as M8 reads it: human rows by taxon, the miRBase name and accession.
FAMILIES = ("miR family\tSeed+m8\tSpecies ID\tMiRBase ID\tMature sequence\tFamily Conservation?"
            "\tMiRBase Accession\nmiR-375\tUUGUUCG\t9606\thsa-miR-375\tUUUGUUCG\t2\tMIMAT0000728"
            "\nmiR-375\tUUGUUCG\t10090\tmmu-miR-375-3p\tUUUGUUCG\t2\tMIMAT0000739\n")  # fmt: skip


def _raw(path: Path) -> dict:
    return tomllib.loads(path.read_text())


def test_template_loads_with_every_module_m10_needs():
    """The template as copied is a dataset file: [gdc] with the CDR, [clinical] with an
    endpoint, one paired tumour-against-normal contrast, [de], [function], [ncrna] with its
    miRNA contrast, and [sanity] with its registration time. Catches a table dropped from it,
    the contrast unpaired, and a key ctx.dataset refuses."""
    spec = load(TEMPLATE)
    assert spec.gdc is not None and spec.gdc.cdr_file == CDR_NAME
    assert spec.clinical is not None and spec.clinical.endpoint == "OS"
    (contrast,) = spec.contrasts
    assert (contrast.reference, contrast.tested, contrast.paired) == ("normal", "tumour", True)
    assert contrast.design == "~donor + group"
    assert spec.de and spec.function and spec.ncrna and spec.ncrna.mirna
    assert spec.ncrna.mirna.paired and "registered_on" in spec.sanity
    assert set(_raw(TEMPLATE)["clinical"]["published"]["table3"]) == {"source", *TABLE3}


def test_template_resources_are_the_kidney_files():
    """Every resource the template declares is the kidney dataset's entry of that name, byte
    for byte, and the [gdc] fields every TCGA project shares are the kidney file's. Catches a
    URL, size or hash edited in one copy alone, and the workflow, gene model or protocol."""
    kirc = {r["name"]: r for r in _raw(KIRC)["resources"]}
    mine = {r["name"]: r for r in _raw(TEMPLATE)["resources"]}
    targets = {"targetscan_predictions", "targetscan_families", "targetscan_genes"}
    assert set(mine) == {"annotation_gtf", "msigdb_hallmark", "mirtarbase_strong", *targets}
    assert all(mine[n] == kirc[n] for n in mine)
    shared = ("workflow", "gene_model", "protocol", "tumour", "normal", "mirna_workflow",
              "cdr_file", "attribution")  # fmt: skip
    mine_gdc, kirc_gdc = _raw(TEMPLATE)["gdc"], _raw(KIRC)["gdc"]
    assert {k: mine_gdc[k] for k in shared} == {k: kirc_gdc[k] for k in shared}
    assert _raw(TEMPLATE)["function"]["licence"] == _raw(KIRC)["function"]["licence"]


def test_template_table_is_its_own_text():
    """The lines --check prints for a missing table are the template's table, n set to the
    files'. Catches the lines drifting from the template, and n left at its placeholder."""
    lines = gdc_check.table("[[contrasts]]", reference=7, tested=9)
    (entry,) = tomllib.loads(lines)["contrasts"]
    assert entry == {**_raw(TEMPLATE)["contrasts"][0], "n": [7, 9]}
    assert tomllib.loads(gdc_check.table("[de]")) == {"de": _raw(TEMPLATE)["de"]}


def _copied(tmp_path: Path, **gdc_values: str) -> Path:
    """The template copied, its placeholders replaced by the made-up project's values."""
    raw = _raw(TEMPLATE)
    text = TEMPLATE.read_text().replace(raw["gdc"]["project"], PROJECT)
    # The made-up counts files hold genes G1 to G12, not the template's MKI67.
    genes = ", ".join(f'"{g}"' for g in raw["clinical"]["genes"])
    text = text.replace(f"genes = [{genes}]", f'genes = ["{GENES[0][1]}"]')
    for key, value in gdc_values.items():
        text = text.replace(f'{key} = "{raw["gdc"][key]}"', f'{key} = "{value}"')
    path = tmp_path / "mine.toml"
    path.write_text(text)
    return path


@pytest.fixture
def folders(tmp_path):
    """The made-up project's three folders, the STAR one with a CDR of three patients."""
    star = tmp_path / "raw" / "fake"
    build(star, STAR)
    small = tuple(mirna_fixture.mirna_aliquot(a) for a in STAR[:4])
    for kind in (gdc_mirna.MIRNA, gdc_mirna.ISOFORM):
        mirna_fixture.build(Path(f"{star}_{kind}"), kind, small)
    cases = sorted({a.case for a in STAR})
    rows = [{"bcr_patient_barcode": c, "type": "BRCA", **{e: v[i] for e, v in FLAGS.items()},
             **{f"{e}.time": 100 for e in FLAGS}} for i, c in enumerate(cases)]  # fmt: skip
    if XLSX:
        write_cdr(star / CDR_NAME, rows, CDR_COLUMNS)
    return star


def _check(tmp_path, star, path):
    """pipeline.check on the three folders, every resource in place under the root."""
    spec = load(path)
    for resource in spec.resources:
        (tmp_path / "root" / resource.dest).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / "root" / resource.dest).touch()
        if resource.name == "targetscan_families":
            with zipfile.ZipFile(tmp_path / "root" / resource.dest, "w") as archive:
                archive.writestr("miR_Family_Info.txt", FAMILIES)
    options = pipeline.Options(spec.slug, path, gdc_dirs=(star,),
                               mirna_dir=Path(f"{star}_mirna"),
                               isoform_dir=Path(f"{star}_isoform"))  # fmt: skip
    problems, modules = pipeline.check(tmp_path / "root", options)
    return problems, modules


@pytest.mark.skipif(not XLSX, reason="openpyxl ==3.1.5 is pinned in pixi.toml; run pixi install")
def test_check_names_each_edit_the_copy_needs(tmp_path, folders):
    """--check on the template as copied names the release the folder holds, the CDR type its
    rows are of, Table 3 still at zeros with the CDR's own counts, and each contrast's n, both
    the STAR pairs and the miRNA pairs; with those set, nothing stops the run. Catches a
    release, type or n the check does not read from the files, and Table 3's zeros passed."""
    path = _copied(tmp_path)
    first, modules = _check(tmp_path, folders, path)
    stops = [p for p in first if p.blocks]
    assert [p.lines for p in stops if p.lines.startswith("release")] == [f'release = "{RELEASE}"']
    assert modules == ["gather", "de", "function", "ncrna", "clinical", "report"]
    path = _copied(tmp_path, release=RELEASE)
    second = [p for p in _check(tmp_path, folders, path)[0] if p.blocks]
    texts = "\n".join(p.text for p in second)
    assert [p.lines for p in second if p.lines] == ['cdr_type = "BRCA"'], texts
    assert "[[contrasts]] tumour_normal: n = [0, 0] declared, [2, 2] in the files" in texts
    assert "[ncrna.mirna] mirna_tumour_normal: n = [0, 0] declared, [2, 2] in the files" in texts
    text = _copied(tmp_path, release=RELEASE, cdr_type="BRCA").read_text()
    text = text.replace("n = [0, 0]", "n = [2, 2]")
    path.write_text(text)
    third = [p for p in _check(tmp_path, folders, path)[0] if p.blocks]
    (table3,) = third
    shown = "; ".join(f"{e} {n[0]} events, {n[1]} censored" for e, n in TABLE3.items())
    assert table3.text.endswith(shown) and "Table 3" in table3.text
    for endpoint, (events, censored) in TABLE3.items():
        text = text.replace(f"{endpoint} = [0, 0]", f"{endpoint} = [{events}, {censored}]")
    path.write_text(text)
    last = _check(tmp_path, folders, path)[0]
    assert [p for p in last if p.blocks] == [], [p.text for p in last if p.blocks]
    assert {p.text.split(":")[0] for p in last} == {"[ncrna] mirnas is empty",
                                                    "no [[sanity.rules]]"}  # fmt: skip


@pytest.mark.skipif(not XLSX, reason="openpyxl ==3.1.5 is pinned in pixi.toml; run pixi install")
def test_check_says_when_table3_differs_and_when_rules_lack_a_time(tmp_path, folders):
    """A Table 3 row that differs from the CDR's counts is named, endpoint by endpoint, without
    stopping the run; rules with a malformed registration time are named too. Catches a
    mismatch passed in silence, and a rule set registered with no time."""
    text = _copied(tmp_path, release=RELEASE, cdr_type="BRCA").read_text()
    text = text.replace("n = [0, 0]", "n = [2, 2]").replace("OS = [0, 0]", "OS = [2, 1]")
    text = text.replace("PFI = [0, 0]", "PFI = [1, 3]").replace("DFI = [0, 0]", "DFI = [0, 2]")
    text = text.replace("DSS = [0, 0]", "DSS = [1, 2]")
    text = text.replace('registered_on = ""', 'registered_on = "today"\n\n[[sanity.rules]]\n'
                        'id = "x"\ntext = "x"\nsource = "x"\nconditions = []')  # fmt: skip
    path = tmp_path / "mine.toml"
    path.write_text(text)
    problems = _check(tmp_path, folders, path)[0]
    assert [p for p in problems if p.blocks] == []
    (differ,) = [p for p in problems if p.text.startswith("[clinical.published.table3]")]
    assert differ.text.startswith("[clinical.published.table3]: PFI differ")
    assert any(p.text.startswith("[sanity] registered_on is 'today'") for p in problems)


@pytest.mark.skipif(not XLSX, reason="openpyxl ==3.1.5 is pinned in pixi.toml; run pixi install")
def test_check_names_a_misspelt_gene_or_mirna(tmp_path, folders):
    """A gene symbol the counts files lack, in [clinical] or a rule, and a miRNA name TargetScan's
    human family rows lack, in [ncrna] or its rules, each stop the run with the closest names.
    Catches a misspelt name left to drop out of a page, or a rule NOT MEASURED, in silence."""
    text = _copied(tmp_path, release=RELEASE, cdr_type="BRCA").read_text()
    text = text.replace("n = [0, 0]", "n = [2, 2]")
    for endpoint, (events, censored) in TABLE3.items():
        text = text.replace(f"{endpoint} = [0, 0]", f"{endpoint} = [{events}, {censored}]")
    text = text.replace(f'genes = ["{GENES[0][1]}"]', 'genes = ["G1", "G13"]')
    text = text.replace("mirnas = []", 'mirnas = ["hsa-miR-375-3p"]')
    text += ('\n[[sanity.rules]]\nid = "x"\ntext = "x"\nsource = "x"\nconditions = [{ metric = '
             '"gene.tumour_normal.G2X.lfc", op = ">", value = 0 }]\n')  # fmt: skip
    path = tmp_path / "mine.toml"
    path.write_text(text)
    stops = sorted(p.text for p in _check(tmp_path, folders, path)[0] if p.blocks)
    heads = ("[[sanity.rules]]: 'G2X' is not a gene symbol of the counts files; the closest are G2",
             "[clinical] genes: 'G13' is not a gene symbol of the counts files; the closest are G3",
             "[ncrna] mirnas and rules: 'hsa-miR-375-3p' is not a human miRBase name in "
             "TargetScan's family file; the closest are hsa-miR-375")  # fmt: skip
    assert len(stops) == len(heads) and all(
        s.startswith(h) for s, h in zip(stops, heads, strict=True)
    ), stops
