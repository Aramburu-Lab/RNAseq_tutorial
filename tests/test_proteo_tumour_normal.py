"""M7's tumour-against-normal design (H.1) on a made-up root, its readers, and the registered rules.

The root: a made-up GDC project of tumours and normals with a PDC biospecimen record, a TSV
protein table named by aliquot, a TMT mapping workbook, and made-up M5 results. Planted: the
"oxphos" class falls 2 log2 units in tumours at protein only, G9 rises 3 at protein, one case is
excluded by the dataset file, and one protein cell is missing. Each test names its mutation.
Nothing here comes from a real dataset.
"""

import io
import json
import re
import subprocess
import tarfile
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from clinical_root import ctx
from gdc_fixture import RELEASE, SPEC, aliquot, build
from proteo_root import resource_toml
from scipy import stats

from ctx.dataset import DatasetError, load
from ctx.proteo import ProteoError, folds, measure, tables
from ctx.report import text

ROOT = Path(__file__).resolve().parents[1]
CASES = 8
OXPHOS, GLYCOLYSIS = ["G1", "G2", "G3", "G4"], ["G5", "G6", "G7", "G8"]
GMT = "SET_OXPHOS\tx\t" + "\t".join(OXPHOS) + "\nSET_GLYCOLYSIS\tx\t" + "\t".join(GLYCOLYSIS)
DESCRIPTORS = ("Index", "NumberPSM")
RULES = """
[proteo.sanity]
registered_at = "2099-01-01 00:00"

[[proteo.sanity.rules]]
id = "T3"
source = "planted"
conditions = [
    { metric = "class.oxphos.protein_lfc.median", op = "<", value = 0 },
    { metric = "class.oxphos_below_other.p", op = "<", value = 0.05 },
]

[[proteo.sanity.rules]]
id = "T5"
source = "planted"
conditions = [{ metric = "protein.tumour_normal.G9.lfc", op = ">=", value = 1.0 }]
"""


def tsv(genes: list[str], samples: list[str], values: np.ndarray) -> bytes:
    """A CPTAC-shaped TSV: quoted text cells, two descriptors, NA for missing."""
    lines = ["\t".join(f'"{h}"' for h in (*DESCRIPTORS, *samples))]
    for i, gene in enumerate(genes):
        cells = ["NA" if np.isnan(v) else f"{v:.6f}" for v in values[i]]
        lines.append("\t".join([f'"{gene}"', "7", *cells]))
    return ("\n".join(lines) + "\n").encode()


def workbook(path: Path, plexes: dict[str, list[str]]) -> None:
    """A TMT mapping workbook: one proteome and one phospho row per plex, a pooled channel."""
    import openpyxl  # noqa: PLC0415 - pinned with M10; only this root writes a workbook

    book = openpyxl.Workbook()
    sheet = book.active
    width = max(len(v) for v in plexes.values()) + 1
    sheet.append(
        [
            "PCC",
            "Folder Name",
            *(h for k in range(width) for h in (f"T{k} Participant ID", f"T{k} Specimen Label")),
        ]
    )
    for name, labels in plexes.items():
        for kind in ("Proteome", "Phosphoproteome"):
            cells = [x for label in [*labels, "pooled sample"] for x in ("case", label)]
            sheet.append(["X", f"{name}_{kind}_X", *cells])
    book.save(path)


def made_root(tmp: Path) -> tuple[Path, Path]:
    """Gather a made-up tumour-normal project, write the protein side and M5's results."""
    items = [aliquot(c, code=code) for c in range(1, CASES + 1) for code in ("01", "11")]
    folder = build(tmp / "gdc", items, seed=4).path
    rows = [
        {
            "aliquot_submitter_id": f"CPT{c}{code}",
            "case_submitter_id": f"FAKE-A1-{c:04d}",
            "sample_type": {"01": "Primary Tumor", "11": "Solid Tissue Normal"}[code],
        }
        for c in range(1, CASES + 1)
        for code in ("01", "11")
    ]
    rows.append(
        {"aliquot_submitter_id": "QC1", "case_submitter_id": "QC", "sample_type": "Not Reported"}
    )
    (folder / "pdc.json").write_text(json.dumps({"data": {"biospecimenPerStudy": rows}}))
    with (folder / "_retrieval_log.tsv").open("a") as log:
        log.write(f"2099-01-01T00:00:00Z\t{RELEASE}\tpdc.json\tmade up\n")
    root = tmp / "root"
    (root / "resources" / "msigdb").mkdir(parents=True)
    (root / "resources" / "msigdb" / "classes.gmt").write_text(GMT + "\n")
    excluded = aliquot(CASES).case
    spec = tmp / "fake_brca.toml"
    spec.write_text(
        SPEC.format(protocol="unstranded", cdr="").replace(
            "[guard]",
            f'biospecimens = "pdc.json"\nexclude_cases = ["{excluded}"]\n'
            'exclude_reason = "made up."\n\n[guard]',
        )
        + resource_toml("table", "resources/table.tar.gz")
        + resource_toml("plexes", "resources/plexes.xlsx")
        + resource_toml("classes", "resources/msigdb/classes.gmt")
        + '\n[[contrasts]]\nname = "tumour_normal"\nvariable = "condition"\nreference = "normal"\n'
        'tested = "tumour"\ndesign = "~donor + group"\npaired = true\n'
        f"n = [{CASES - 1}, {CASES - 1}]\n"
        '\n[proteo]\ndesign = "tumour_normal"\nprotein = "table"\nmember = "t.tsv"\n'
        f'descriptors = {json.dumps(list(DESCRIPTORS))}\nunit = "log2 abundance"\nplex = "plex"\n'
        'plex_map = "plexes"\nplex_rows = "_Proteome_"\ncontrast = "tumour_normal"\n'
        'genes = ["G9"]\npositive_fdr = 0.01\n'
        '[proteo.classes.oxphos]\nresource = "classes"\nsets = ["SET_OXPHOS"]\n'
        '[proteo.classes.glycolysis]\nresource = "classes"\nsets = ["SET_GLYCOLYSIS"]\n' + RULES
    )
    ctx(
        "--root",
        str(root),
        "gather",
        "--dataset",
        "fake_brca",
        "--gdc-dir",
        str(folder),
        "--dataset-file",
        str(spec),
    )
    genes = [f"G{i}" for i in range(1, 13)]
    names = [r["aliquot_submitter_id"] for r in rows]
    rng = np.random.default_rng(11)
    base = rng.normal(20, 1, (len(genes), 1)) + rng.normal(0, 0.3, (len(genes), len(names)))
    tumour = np.array([n.endswith("01") for n in names])
    for gene, shift in (*((g, -2.0) for g in OXPHOS), ("G9", 3.0)):
        base[genes.index(gene), tumour] += shift
    base[genes.index("G9"), names.index("CPT101")] = np.nan
    body = tsv(genes, names, base)
    with tarfile.open(root / "resources" / "table.tar.gz", "w:gz") as tar:
        info = tarfile.TarInfo("x/t.tsv")
        info.size = len(body)
        tar.addfile(info, io.BytesIO(body))
    workbook(
        root / "resources" / "plexes.xlsx",
        {"01P": names[: len(names) // 2], "02P": names[len(names) // 2 :]},
    )
    de = root / "outputs" / "05_de"
    de.mkdir(parents=True)
    pl.DataFrame(
        {
            "gene_id": [f"ENSG{i}" for i in range(12)],
            "gene_name": genes,
            "log2FoldChange": [0.1] * 12,
            "padj": [0.5] * 12,
        }
    ).write_parquet(
        de / f"{date.today().isoformat()}_de_fake_brca_tumour_normal_results_v1.parquet"
    )
    ctx("--root", str(root), "proteo", "--dataset", "fake_brca", "--dataset-file", str(spec))
    return root, spec


@pytest.fixture(scope="module")
def made(tmp_path_factory) -> tuple[Path, Path]:
    """One made-up root for every end-to-end test here."""
    return made_root(tmp_path_factory.mktemp("tn"))


def out(root: Path, what: str) -> pl.DataFrame:
    """One table ctx proteo wrote."""
    return pl.read_parquet(next((root / "outputs" / "07_proteo").glob(f"*_{what}_v1.parquet")))


def test_planted_shifts_found_on_the_paired_cases(made):
    """Catches: Welch in place of the paired t (n would not be the pairs); the excluded case's
    columns used (8 pairs, not 7); the oxphos test turned round (above, not below)."""
    root, _ = made
    metrics = dict(out(root, "metrics").select("metric", "value").rows())
    assert metrics["class.oxphos.protein_lfc.median"] == pytest.approx(-2.0, abs=0.3)
    assert metrics["class.oxphos_below_other.p"] < 0.05
    assert metrics["protein.tumour_normal.G9.lfc"] == pytest.approx(3.0, abs=0.3)
    details = dict(out(root, "metrics").select("metric", "detail").rows())
    assert details["protein.tumour_normal.G9.lfc"] == f"paired t over {CASES - 1} pairs"
    verdicts = dict(out(root, "sanity").select("rule", "verdict").rows())
    assert verdicts == {"T3": "PASS", "T5": "PASS"}


def test_both_correlations_on_their_own_samples(made):
    """Catches: the normals' correlation computed on the tumours (equal columns), and the
    excluded case's tumour counted (n_quantified 8)."""
    root, _ = made
    genes = out(root, "genes")
    assert {"normal_rho", "normal_rho_padj"} <= set(genes.columns)
    assert genes["n_quantified"].max() == CASES - 1
    assert not genes["rho"].equals(genes["normal_rho"])
    metrics = dict(out(root, "metrics").select("metric", "value").rows())
    assert {"pergene.tumour.positive_share", "pergene.normal.positive_share"} <= set(metrics)


def test_normals_draw_their_own_class_figure(made):
    """S.4: the normals' per-gene rho by class is drawn from the normal_ columns.

    Catches: the normals' figure drawn from the tumours' rho; the normals' class medians
    taken over every gene rather than those tested in normals; the figure left out."""
    root, _ = made
    genes = out(root, "genes").filter(pl.col("normal_tested") & (pl.col("gene_class") != ""))
    (data,) = (root / "outputs" / "07_proteo").glob("*_rho_by_class_normal_v1.csv")
    drawn = pl.read_csv(data).sort("gene")
    expected = genes.select("gene", rho="normal_rho").drop_nulls().sort("gene")
    assert drawn["gene"].to_list() == expected["gene"].to_list()
    assert np.allclose(drawn["rho"].to_numpy(), expected["rho"].to_numpy())
    metrics = dict(out(root, "metrics").select("metric", "value").rows())
    oxphos = genes.filter(pl.col("gene_class") == "oxphos")["normal_rho"].drop_nulls()
    assert metrics["normal.class.oxphos.median"] == pytest.approx(float(np.median(oxphos)))


def test_normal_class_medians_read_the_genes_tested_in_normals():
    """Catches: the normals' class medians taken over the genes tested in tumours (the fixture
    above tests the same genes in both, so it cannot tell them apart)."""
    from types import SimpleNamespace  # noqa: PLC0415 - this test alone fakes a spec

    from ctx.proteo import run  # noqa: PLC0415

    genes = pl.DataFrame({
        "gene": [f"G{i}" for i in range(8)], "gene_class": ["oxphos"] * 4 + [measure.OTHER] * 4,
        "tested": [True, True, False, False] * 2, "normal_tested": [False, False, True, True] * 2,
        "normal_rho": [0.9, 0.8, 0.1, 0.2, 0.5, 0.4, 0.3, 0.6], "normal_r": [0.0] * 8,
    })  # fmt: skip
    measured = run.Measured(None, [], genes, np.zeros(8))  # type: ignore[arg-type]
    run.normal_metrics(measured, SimpleNamespace(classes={"oxphos": None}), "spearman")  # type: ignore[arg-type]
    metrics = {m: v for m, v, _ in measured.metrics}
    assert metrics["normal.class.oxphos.median"] == pytest.approx(0.15)
    assert metrics["normal.pergene.spearman.median"] == pytest.approx(0.25)


def test_page_text_is_written_for_tumour_against_normal(made):
    """S.4: the Interpretation, Questions and Pitfalls are the design's, not M7's generic text.

    Catches: the page built from the generic sections (granules, complexes)."""
    root, _ = made
    (page,) = (root / "outputs" / "07_proteo").glob("*_overview_v1.html")
    html = page.read_text()
    own = text.raw("proteo", "variants.tumour_normal")
    assert own["pitfalls"][1].split(":")[0] in html
    assert own["questions"][1]["prompt"].split(" the same")[0] in html
    assert "granule" not in html.lower()
    assert "Per-gene RNA-protein correlation in normals" in html


def test_paired_calls_leave_a_missing_value_out_never_zero():
    """Catches: NaN read as 0 (G's difference then 20 units off); the share counted over one
    group; differences taken unpaired (column order shuffled)."""
    rng = np.random.default_rng(2)
    tumour = rng.normal(10, 1, (3, 6))
    normal = tumour - np.array([[1.0], [0.0], [2.0]]) + rng.normal(0, 0.1, (3, 6))
    tumour[0, 0] = np.nan
    found = folds.paired_calls(["A", "B", "C"], tumour, normal, 0.7, 1.0)
    keep = ~np.isnan(tumour[0])
    expected = stats.ttest_rel(tumour[0, keep], normal[0, keep])
    assert found["protein_lfc"][0] == pytest.approx(np.mean(tumour[0, keep] - normal[0, keep]))
    assert found["protein_p"][0] == pytest.approx(expected.pvalue)
    assert found["quantified_share"][0] == pytest.approx(11 / 12)
    normal[1, 0] = normal[1, 1] = normal[1, 2] = normal[1, 3] = np.nan
    assert not folds.paired_calls(["A", "B", "C"], tumour, normal, 0.7, 1.0)["tested"][1]


def test_positive_share_counts_defined_genes_only():
    """Catches: untested genes (null rho) in the denominator; a negative rho counted."""
    frame = pl.DataFrame({"rho": [0.5, -0.6, 0.4, None], "rho_padj": [0.001, 0.001, 0.2, None]})
    assert measure.positive_share(frame, "", 0.01) == (pytest.approx(1 / 3), 1, 3)


def test_tsv_reader_refuses_another_header_and_reads_na(tmp_path):
    """Catches: a table read whatever its header (samples shifted onto descriptors); NA as 0."""
    body = tsv(["A", "B"], ["S1", "S2"], np.array([[1.0, np.nan], [2.0, 3.0]]))
    archive = tmp_path / "t.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("t.tsv")
        info.size = len(body)
        tar.addfile(info, io.BytesIO(body))
    layer = tables.read_tsv(archive, "t.tsv", DESCRIPTORS)
    assert layer.samples == ("S1", "S2") and np.isnan(layer.values[0, 1])
    with pytest.raises(ProteoError, match="header"):
        tables.read_tsv(archive, "t.tsv", ("Index",) + ("Other",))


def test_plex_map_keeps_the_proteome_rows_only(tmp_path):
    """Catches: phospho rows read as plexes (every label in two); the pooled label in one."""
    workbook(tmp_path / "w.xlsx", {"01P": ["a", "b"], "02P": ["c"]})
    plexes = tables.read_plexes(tmp_path / "w.xlsx", "_Proteome_")
    assert plexes["a"] == {"01P_Proteome_X"} and plexes["c"] == {"02P_Proteome_X"}
    assert len(plexes["pooled sample"]) == 2


def test_tumour_normal_needs_its_keys(tmp_path):
    """Catches: the design accepted without a plex map or the biospecimens gather reads."""
    path = tmp_path / "x.toml"
    path.write_text(
        SPEC.format(protocol="unstranded", cdr="")
        + resource_toml("t", "t")
        + '\n[proteo]\ndesign = "tumour_normal"\nprotein = "t"\n'
    )
    with pytest.raises(DatasetError, match=r"plex_map.*|biospecimens"):
        load(path)


def _rules_block(text: str) -> str:
    """The [proteo.sanity] tables of a TOML text, from the header to the end of the last rule."""
    found = re.search(r"^\[proteo\.sanity\]$", text, re.M)
    assert found is not None
    block = text[found.start() :]
    end = block.find("\n```")
    return (block if end < 0 else block[:end]).strip()


def test_ccrcc_rules_in_force_equal_the_registered_ones():
    """The rules C0 to C6 were registered in the plan's first commit before any M7 run on this
    cohort (lessons.md: compare content, not only commit order). Catches any edit to a rule,
    threshold or metric name after the registration."""
    plan = "docs/plans/2026-10-14_cptac_ccrcc_m7_rules.md"
    first = subprocess.run(
        ["git", "log", "--diff-filter=A", "--format=%H", "--", plan],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    ).stdout.split()
    if not first:
        pytest.skip("no git history here to read the registration from")
    registered = subprocess.run(
        ["git", "show", f"{first[-1]}:{plan}"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout
    in_force = (ROOT / "conf" / "datasets" / "cptac_ccrcc.toml").read_text()
    assert re.sub(r"\s+", " ", _rules_block(registered)) == re.sub(
        r"\s+", " ", _rules_block(in_force)
    )
