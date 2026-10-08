"""``ctx proteo`` end to end on made-up roots, both designs, and the rules' registration order.

Across samples: a made-up GDC project gathered by ``ctx gather --gdc-dir`` (tests/gdc_fixture.py)
and a made-up freeze archive whose protein follows RNA for one planted class and not another.
Across studies: made-up M5 results and a made-up copy-number workbook with a planted granule
class. Each test names the mutation it catches.
"""

import io
import re
import subprocess
import tarfile
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from clinical_root import ctx
from gdc_fixture import SPEC, aliquot, build

from ctx.de.filters import log2cpm
from ctx.report.builder import MAX_BYTES

ROOT = Path(__file__).resolve().parents[1]
TUMOURS = 30
COMPLEX, METABOLIC = ["G1", "G2", "G3", "G4"], ["G5", "G6", "G7", "G8"]
RULES = """
[proteo.sanity]
registered = "2099-01-01 00:00"

[[proteo.sanity.rules]]
id = "K2"
source = "planted"
conditions = [{{ metric = "class.complex_below_metabolic.p", op = "<", value = 0.05 }}]

[[proteo.sanity.rules]]
id = "K3"
source = "planted"
conditions = [{{ metric = "page.missingness_first", op = "==", value = 1 }}]
{extra}"""
GMT = "SET_COMPLEX\tx\t" + "\t".join(COMPLEX) + "\nSET_METABOLIC\tx\t" + "\t".join(METABOLIC)


def gct(genes: list[str], samples: list[str], values: np.ndarray, annotation: dict) -> bytes:
    """A GCT 1.3 table: one row descriptor (geneSymbol), the column annotations, NA for NaN."""
    lines = ["#1.3", f"{len(genes)}\t{len(samples)}\t1\t{len(annotation)}",
             "\t".join(["id", "geneSymbol", *samples])]  # fmt: skip
    lines += ["\t".join([name, "na", *map(str, vals)]) for name, vals in annotation.items()]
    for i, gene in enumerate(genes):
        cells = ["NA" if np.isnan(v) else f"{v:.6f}" for v in values[i]]
        lines.append("\t".join([f"NP_{i}", gene, *cells]))
    return "\n".join(lines).encode()


def resources(root: Path) -> None:
    """The class GMT where the dataset file says ``ctx setup`` puts it."""
    (root / "resources" / "msigdb").mkdir(parents=True, exist_ok=True)
    (root / "resources" / "msigdb" / "classes.gmt").write_text(GMT + "\n")


def resource_toml(name: str, dest: str) -> str:
    """A [[resources]] entry; ctx proteo reads only the name and dest."""
    return (f'\n[[resources]]\nname = "{name}"\nurls = ["https://example.org/{name}"]\n'
            f'bytes = 1\nsha256 = "{"0" * 64}"\ndest = "{dest}"\n')  # fmt: skip


@pytest.fixture(scope="module")
def across_root(tmp_path_factory):
    """A gathered made-up project and a freeze archive; ctx proteo run on it."""
    tmp = tmp_path_factory.mktemp("m7_across")
    items = [aliquot(i) for i in range(1, TUMOURS + 1)]
    build(tmp / "gdc", items, seed=3)
    root = tmp / "root"
    proteo = (resource_toml("freeze", "resources/freeze.tar.gz")
              + resource_toml("classes", "resources/msigdb/classes.gmt")
              + '\n[proteo]\ndesign = "across_samples"\nprotein = "freeze"\nmember = "p.gct"\n'
              'sample_prefix = "X"\nplex = "TMT.Plex"\ntracks = ["PAM50"]\n'
              'amplified = "ERBB2.Gene.Amplified"\ngenes = ["G9"]\n'
              '[proteo.classes.complex]\nresource = "classes"\nsets = ["SET_COMPLEX"]\n'
              '[proteo.classes.metabolic]\nresource = "classes"\nsets = ["SET_METABOLIC"]\n'
              + RULES.format(extra=""))  # fmt: skip
    spec = tmp / "fake_brca.toml"
    spec.write_text(SPEC.format(protocol="unstranded", cdr="") + proteo)
    ctx("--root", str(root), "gather", "--dataset", "fake_brca", "--gdc-dir", str(tmp / "gdc"),
        "--dataset-file", str(spec))  # fmt: skip
    gather = root / "outputs" / "01_gather"
    sheet = pl.read_parquet(next(gather.glob("*_samplesheet_v1.parquet")))
    counts = pl.read_parquet(next(gather.glob("*_counts_v1.parquet")))
    names = pl.read_parquet(next(gather.glob("*_genes_v1.parquet")))["gene_name"].to_list()
    samples = sheet["sample"].to_list()
    rna = log2cpm(counts.select(samples).to_numpy(), sheet["library_size"].to_numpy())
    rng = np.random.default_rng(8)
    # Metabolic genes and G9 follow RNA; complex subunits do not; G10 is missing in half.
    follows = np.array([n in METABOLIC or n == "G9" for n in names])
    protein = np.where(follows[:, None], rna - rna.mean(axis=1, keepdims=True), 0.0)
    protein += rng.normal(0, np.where(follows, 0.02, 1.0)[:, None], rna.shape)
    protein[names.index("G10"), ::2] = np.nan
    annotation = {"TMT.Plex": [f"P{j % 3}" for j in range(len(samples))],
                  "PAM50": [("LumA", "Basal")[j % 2] for j in range(len(samples))],
                  "ERBB2.Gene.Amplified": [str(j % 2) for j in range(len(samples))]}  # fmt: skip
    # The freeze lists its samples in another order than gather's sheet: joined by case.
    cases = ["X" + c for c in sheet["case"].to_list()][::-1]
    protein = np.vstack([protein, np.full(len(samples), np.nan)])[:, ::-1]
    body = gct([*names, "G1"], cases, protein, annotation)
    (root / "resources").mkdir(exist_ok=True)
    with tarfile.open(root / "resources" / "freeze.tar.gz", "w:gz") as tar:
        info = tarfile.TarInfo("freeze/p.gct")
        info.size = len(body)
        tar.addfile(info, io.BytesIO(body))
    resources(root)
    ctx("--root", str(root), "proteo", "--dataset", "fake_brca", "--dataset-file", str(spec))
    return root


@pytest.fixture(scope="module")
def cross_root(tmp_path_factory):
    """Made-up M5 results and copy numbers with a planted granule class; ctx proteo run."""
    openpyxl = pytest.importorskip("openpyxl")
    tmp = tmp_path_factory.mktemp("m7_cross")
    root = tmp / "root"
    rng = np.random.default_rng(5)
    n = 120
    genes = [f"G{i}" for i in range(n)]
    rna_lfc = rng.normal(0, 2, n)
    granule = set(genes[:20])
    de = root / "outputs" / "05_de"
    de.mkdir(parents=True)
    pl.DataFrame({"gene_id": [f"ENSG{i}" for i in range(n)], "gene_name": genes,
                  "log2FoldChange": rna_lfc, "padj": rng.uniform(0, 1, n),
                  "baseMean": rng.uniform(10, 1000, n)}).write_parquet(
        de / f"{date.today().isoformat()}_de_fake_cells_a_vs_b_results_v1.parquet")  # fmt: skip
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.append(["Majority protein IDs", "Gene names",
                  *[f"CopyNumber_A_{d:02d}" for d in range(1, 5)],
                  *[f"CopyNumber_B_{d:02d}" for d in range(1, 5)]])  # fmt: skip
    for i, gene in enumerate(genes):
        shift = rna_lfc[i] + (3.0 if gene in granule else 0.0)
        a = 2 ** (10 + shift / 2 + rng.normal(0, 0.1, 4))
        b = 2 ** (10 - shift / 2 + rng.normal(0, 0.1, 4))
        row = [*a.round(1), *b.round(1)]
        if i == 30:
            row[0] = row[1] = row[4] = 0  # 5 of 8: below the 70 % share, untested
        sheet.append([f"P{i}", gene, *row])
    (root / "resources").mkdir()
    book.save(root / "resources" / "copies.xlsx")
    resources(root)
    granule_gmt = "SET_GRANULE\tx\t" + "\t".join(sorted(granule))
    (root / "resources" / "msigdb" / "classes.gmt").write_text(granule_gmt + "\n")
    proteo = (resource_toml("copies", "resources/copies.xlsx")
              + resource_toml("classes", "resources/msigdb/classes.gmt")
              + '\n[[contrasts]]\nname = "a_vs_b"\nvariable = "cell_type"\nreference = "B"\n'
              'tested = "A"\nn = [3, 3]\n'
              '\n[proteo]\ndesign = "cross_dataset"\nprotein = "copies"\ncontrast = "a_vs_b"\n'
              "tested = '^CopyNumber_A_'\nreference = '^CopyNumber_B_'\n"
              '[proteo.classes.granule]\nresource = "classes"\nsets = ["SET_GRANULE"]\n'
              '\n[proteo.sanity]\nregistered = "2099-01-01 00:00"\n'
              '[[proteo.sanity.rules]]\nid = "P1"\nsource = "planted"\n'
              'conditions = [{ metric = "fold.spearman.low", op = ">", value = 0 }]\n'
              '[[proteo.sanity.rules]]\nid = "P2"\nsource = "planted"\n'
              'conditions = [{ metric = "class.granule.p", op = "<", value = 0.05 }]\n'
              )  # fmt: skip
    spec = tmp / "fake_cells.toml"
    spec.write_text('[dataset]\nslug = "fake_cells"\n[guard]\nids = []\npatterns = []\n' + proteo)
    ctx("--root", str(root), "proteo", "--dataset", "fake_cells", "--dataset-file", str(spec))
    return root


def output(root: Path, what: str) -> Path:
    """The one proteo output of this kind."""
    (path,) = (root / "outputs" / "07_proteo").glob(f"*_proteo_*_{what}")
    return path


@pytest.mark.parametrize("which", ["across_root", "cross_root"])
def test_page_offline_and_every_table_described(which, request):
    """Catches a page over the cap or fetching a script, and a table without its README."""
    root = request.getfixturevalue(which)
    html = output(root, "overview_v1.html").read_text()
    assert len(html.encode()) <= MAX_BYTES
    assert not re.search(r"<script[^>]+\bsrc\s*=", html, re.IGNORECASE)
    for table in (root / "outputs" / "07_proteo").glob("*.csv"):
        assert table.with_name(table.stem + "_README.md").is_file(), table.name


@pytest.mark.parametrize("which", ["across_root", "cross_root"])
def test_missingness_first_and_rules_pass(which, request):
    """Catches the class or fold figure drawn before the missingness figures (K3 fails), and
    the planted class shift or fold agreement lost (K2, P1 or P2 fails)."""
    root = request.getfixturevalue(which)
    verdicts = pl.read_parquet(output(root, "sanity_v1.parquet"))
    assert set(verdicts["rule_verdict"]) == {"PASS"}, verdicts.select("rule", "left").rows()
    html = output(root, "overview_v1.html").read_text()
    stems = re.findall(r"_proteo_fake_\w+?_(\w+?)_v1\.csv", html)
    assert stems[:2] == ["quantified", "coverage"], stems


def test_across_missing_never_zero_and_one_row_per_gene(across_root):
    """Catches a missing protein read as 0 (G10 would count every sample), the extra G1 row
    kept as a second gene, and samples joined by position instead of case."""
    genes = pl.read_parquet(output(across_root, "genes_v1.parquet"))
    assert genes["gene"].n_unique() == genes.height
    g10 = genes.filter(pl.col("gene") == "G10").row(0, named=True)
    assert g10["n_quantified"] == TUMOURS // 2 and not g10["tested"]
    g9 = genes.filter(pl.col("gene") == "G9").row(0, named=True)
    assert g9["rho"] > 0.95


def test_cross_untested_below_the_share(cross_root):
    """Catches the 70 % filter skipped: G30 has 5 of 8 values."""
    genes = pl.read_parquet(output(cross_root, "genes_v1.parquet"))
    g30 = genes.filter(pl.col("gene") == "G30").row(0, named=True)
    assert g30["n_samples"] == 5 and not g30["tested"]
    html = output(cross_root, "overview_v1.html").read_text()
    assert "Cross-study comparison: RNA from fake_cells" in html


def git(*args: str) -> str:
    """git in this checkout."""
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True,
                          check=True).stdout.strip()  # fmt: skip


@pytest.mark.parametrize("name", ["tier1", "cptac_brca"])
def test_rules_registered_before_the_protein_reader(name):
    """The registration commit of [proteo.sanity] precedes the commit that added the reader.

    Catches a rule edited after protein values could be read without registering again."""
    if git("rev-parse", "--is-shallow-repository") == "true":
        pytest.skip("shallow clone: the history is absent")
    path = f"conf/datasets/{name}.toml"
    registered = git("log", "--format=%H", "-G", r"^registered_at = \"2026-10", "--", path).split()
    reader = git("log", "--reverse", "--diff-filter=A", "--format=%H", "--",
                 "src/ctx/proteo/tables.py").split()  # fmt: skip
    assert registered and reader
    assert registered[0] != reader[0]
    ancestor = subprocess.run(["git", "-C", str(ROOT), "merge-base", "--is-ancestor",
                               registered[0], reader[0]], check=False)  # fmt: skip
    assert ancestor.returncode == 0
