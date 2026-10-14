"""Made-up roots for ``ctx proteo``, both designs, for tests/test_proteo_page.py and test_cli.py.

Across samples: a made-up GDC project gathered by ``ctx gather --gdc-dir`` (tests/gdc_fixture.py)
and a made-up freeze archive whose protein follows RNA for one planted class and not another.
Across studies: made-up M5 results and a made-up copy-number workbook with a planted granule
class. Nothing here comes from a real dataset.
"""

import io
import tarfile
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
from clinical_root import ctx
from gdc_fixture import SPEC, aliquot, build

from ctx.de.filters import log2cpm

TUMOURS = 30
COMPLEX, METABOLIC = ["G1", "G2", "G3", "G4"], ["G5", "G6", "G7", "G8"]
RULES = """
[proteo.sanity]
registered_at = "2099-01-01 00:00"

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


def subtypes(sheet: pl.DataFrame) -> list[str]:
    """Each sheet tumour's made-up PAM50 call, alternating in sheet order."""
    return [("LumA", "Basal")[i % 2] for i in range(sheet.height)]


def results(root: Path, names: list[str]) -> None:
    """Made-up M5 results for basal_vs_luma: G9 and G11 up at RNA, G10 up but untested at
    protein, G5 up with |log2 FC| under 1, the rest unchanged."""
    planted = {"G9": (2.0, 1e-6), "G11": (2.5, 1e-6), "G10": (2.0, 1e-6), "G5": (0.9, 1e-6)}
    lfc = [planted.get(n, (0.1, 0.0))[0] for n in names]
    padj = [planted.get(n, (0.0, 0.6))[1] for n in names]
    de = root / "outputs" / "05_de"
    de.mkdir(parents=True, exist_ok=True)
    name = f"{date.today().isoformat()}_de_fake_brca_basal_vs_luma_results_v1.parquet"
    pl.DataFrame({"gene_id": [f"ENSG{i}" for i in range(len(names))], "gene_name": names,
                  "log2FoldChange": lfc, "padj": padj}).write_parquet(de / name)  # fmt: skip


def across(tmp: Path) -> Path:
    """A gathered made-up project and a freeze archive under ``tmp``; ctx proteo run on it.

    The sheet gains ``pam50`` after gather, as the declared annotation join is to write it, and
    M5's results for basal_vs_luma are made up (``results``); the freeze names one column by its
    GDC sample ID, as CPTAC's names two.
    """
    items = [aliquot(i) for i in range(1, TUMOURS + 1)]
    build(tmp / "gdc", items, seed=3)
    root = tmp / "root"
    proteo = (resource_toml("freeze", "resources/freeze.tar.gz")
              + resource_toml("classes", "resources/msigdb/classes.gmt")
              + '\n[[contrasts]]\nname = "basal_vs_luma"\nvariable = "pam50"\n'
              f'reference = "LumA"\ntested = "Basal"\nn = [{TUMOURS // 2}, {TUMOURS // 2}]\n'
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
    sheet_file = next(gather.glob("*_samplesheet_v1.parquet"))
    sheet = pl.read_parquet(sheet_file)
    pam50 = subtypes(sheet)
    sheet.with_columns(pam50=pl.Series(pam50)).write_parquet(sheet_file)
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
    # G11's protein is 3 log2 units higher in Basal tumours: called up at both layers.
    protein[names.index("G11")] += 3.0 * (np.array(pam50) == "Basal")
    annotation = {"TMT.Plex": [f"P{j % 3}" for j in range(len(samples))][::-1],
                  "PAM50": pam50[::-1],
                  "ERBB2.Gene.Amplified": [str(j % 2) for j in range(len(samples))]}  # fmt: skip
    # The freeze lists its samples in another order than gather's sheet: joined by case, and
    # its first column by GDC sample ID alone.
    cases = ["X" + c for c in sheet["case"].to_list()][::-1]
    cases[0] = sheet["sample_id"][-1]
    protein = np.vstack([protein, np.full(len(samples), np.nan)])[:, ::-1]
    results(root, names)
    body = gct([*names, "G1"], cases, protein, annotation)
    (root / "resources").mkdir(exist_ok=True)
    with tarfile.open(root / "resources" / "freeze.tar.gz", "w:gz") as tar:
        info = tarfile.TarInfo("freeze/p.gct")
        info.size = len(body)
        tar.addfile(info, io.BytesIO(body))
    resources(root)
    ctx("--root", str(root), "proteo", "--dataset", "fake_brca", "--dataset-file", str(spec))
    return root


def cross(tmp: Path) -> Path:
    """Made-up M5 results and copy numbers with a planted granule class; ctx proteo run."""
    import openpyxl  # noqa: PLC0415 - pinned with M10; only this root writes a workbook

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
              '\n[proteo.sanity]\nregistered_at = "2099-01-01 00:00"\n'
              '[[proteo.sanity.rules]]\nid = "P1"\nsource = "planted"\n'
              'conditions = [{ metric = "fold.spearman.low", op = ">", value = 0 }]\n'
              '[[proteo.sanity.rules]]\nid = "P2"\nsource = "planted"\n'
              'conditions = [{ metric = "class.granule.p", op = "<", value = 0.05 }]\n'
              )  # fmt: skip
    spec = tmp / "fake_cells.toml"
    spec.write_text('[dataset]\nslug = "fake_cells"\n[guard]\nids = []\npatterns = []\n' + proteo)
    ctx("--root", str(root), "proteo", "--dataset", "fake_cells", "--dataset-file", str(spec))
    return root
