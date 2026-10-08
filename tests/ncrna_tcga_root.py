"""A made-up GDC ctx root for the TCGA half of ``ctx ncrna``, read by tests/test_ncrna_tcga.py.

Gather's STAR tables shaped as ``ctx gather --gdc-dir`` writes them, its miRNA tables as
``--mirna-dir`` writes them, M5's results, the annotation, the four target files and a
dataset file. STAR counts are tests/ncrna_root.py's, with its planted lncRNA pair, for twelve
patients with a tumour and a normal each, plus a second tumour aliquot of the first patient
that gather's rule did not take and that sorts first by run accession: a reader that bypasses
``gdc.chosen`` takes it in place of the patient's chosen tumour. Planted: hsa-miR-1-3p rises
in tumours and its validated target G00010 falls as it rises, tumour by tumour; hsa-miR-2-5p
falls. Nothing comes from a real dataset.
"""

import hashlib
import zipfile
from datetime import date
from pathlib import Path

import ncrna_root
import numpy as np
import polars as pl
from ncrna_root import GTF, N_GENES, gene_type, write_gtf

from ctx.de.stats import classify
from ctx.io.naming import stem
from ctx.io.readme_companion import write_table

DATASET = "fakegdc"
PAIRS = 12
GROUPS = (("normal", "short", PAIRS, "tissue", "p"), ("tumour", "short", PAIRS, "tissue", "p"))
MIR_UP, MIR_DOWN = "hsa-miR-1-3p", "hsa-miR-2-5p"
MIMATS = [f"MIMAT9{i:06d}" for i in range(1, 31)]
TARGET, PREDICTED, WEAK = "g00010", "g00011", "g00012"
EXTRA = "S999"
TARGETS = "data/raw/mirna_targets"
# GTF miRNA genes are ncrna_root's 352 to 371, on chr2 at 10,000 i + 1.
HAIRPINS = range(352, 372)


def _write(root: Path, module: str, what: str, frame: pl.DataFrame, day: date) -> None:
    """One table as the module writes it, with its README."""
    out = root / "outputs" / {"gather": "01_gather", "de": "05_de"}[module]
    write_table(frame, out, stem(module, f"{DATASET}_{what}", 1, day=day),
                columns={c: ("none", c) for c in frame.columns}, data_releases={})  # fmt: skip


def star(day: date, root: Path) -> tuple[pl.DataFrame, pl.DataFrame, list[str]]:
    """The STAR sheet, counts and genes, the extra aliquot in; returns sheet, counts, names."""
    frame, synthetic = ncrna_root.counts(GROUPS)
    sheet = synthetic.samples
    tumour = (sheet["condition"] == "tumour").to_numpy()
    rank = np.cumsum(tumour) - 1
    # ncrna_root's miRNA genes are near silent; here they vary, as hairpin loci do in STAR
    # counts, so each sample ranks them and the agreement with the small-RNA library is measured.
    rng = np.random.default_rng(9)
    hairpin_genes = [f"g{i + 1:05d}" for i in HAIRPINS]
    frame = frame.with_columns(
        pl.when(pl.col("gene_id") == TARGET).then(pl.lit(int(4000 // (rank[j] + 1) + 50)
                                                       if tumour[j] else 600))
        .when(pl.col("gene_id").is_in(hairpin_genes))
        .then(pl.Series(rng.integers(0, 40, frame.height)))
        .otherwise(pl.col(c)).cast(pl.Int32).alias(c)
        for j, c in enumerate(sheet["sample"])
    )  # fmt: skip
    gdc = sheet.drop("cell_type").with_columns(
        case=pl.col("donor"),
        sample_id=pl.col("donor") + pl.when(pl.col("condition") == "tumour")
        .then(pl.lit("-01A")).otherwise(pl.lit("-11A")),
        sample_type=pl.when(pl.col("condition") == "tumour").then(pl.lit("Primary Tumor"))
        .otherwise(pl.lit("Solid Tissue Normal")),
        preservation=pl.lit("OCT"), redacted=pl.lit(value=False), chosen=pl.lit(value=True),
        drop_reason=pl.lit(""),
    )  # fmt: skip
    first = gdc.filter(pl.col("condition") == "tumour").head(1)
    extra = first.with_columns(sample=pl.lit(EXTRA), run_accession=pl.lit("XRR000000"),
                               chosen=pl.lit(value=False),
                               drop_reason=pl.lit("not_first"))  # fmt: skip
    gdc = pl.concat([gdc, extra])
    frame = frame.with_columns(pl.col(first["sample"][0]).alias(EXTRA) // 3)
    ids = frame["gene_id"].to_list()
    names = [g.upper() for g in ids]
    genes = pl.DataFrame({"gene_id": ids, "gene_id_version": [f"{g}.1" for g in ids],
                          "gene_name": names,
                          "gene_type": [gene_type(i) for i in range(N_GENES)]})  # fmt: skip
    for what, table in (("samplesheet", gdc), ("counts", frame), ("genes", genes)):
        _write(root, "gather", what, table, day)
    return gdc, frame, names


def mirna(day: date, root: Path, sheet: pl.DataFrame) -> None:
    """The miRNA sheet, mature and hairpin counts and the hairpins' loci."""
    rng = np.random.default_rng(5)
    taken = sheet.filter(pl.col("chosen"))
    tumour = (taken["condition"] == "tumour").to_numpy()
    rank = np.cumsum(tumour) - 1
    base = np.exp(rng.normal(np.log(500.0), 1.0, len(MIMATS)))
    mature = rng.poisson(base[:, None] * np.ones(taken.height)).astype(np.int64)
    mature[0] = np.where(tumour, 400 * (rank + 1), 200)
    mature[1] = np.where(tumour, 100, 2000)
    samples = [f"M{s}" for s in taken["sample"]]
    hairpins = rng.poisson(300, (len(HAIRPINS), taken.height))
    mirna_sheet = taken.select(
        "case", "sample_id", "sample_type", "condition", "donor", "preservation", "redacted",
        "chosen", "drop_reason",
        sample=pl.Series(samples), rna_sample="sample",
        run_accession=pl.Series([f"GDC{i + 1:06d}" for i in range(taken.height)]),
        mirna_reads=pl.Series(mature.sum(axis=0) + hairpins.sum(axis=0)),
    )  # fmt: skip
    names = [f"hsa-mir-x{i}" for i in HAIRPINS]
    loci = pl.DataFrame({"miRNA_ID": names, "chrom": "chr2",
                         "start": [10_000 * i + 101 for i in HAIRPINS],
                         "end": [10_000 * i + 180 for i in HAIRPINS], "strand": "+",
                         "mimats": ""})  # fmt: skip
    tables = (
        ("mirna_samplesheet", mirna_sheet),
        ("mirna_mature_counts", pl.DataFrame({"mimat": MIMATS} | {
            s: pl.Series(mature[:, j], dtype=pl.Int32) for j, s in enumerate(samples)})),
        ("mirna_hairpin_counts", pl.DataFrame({"hairpin": names} | {
            s: pl.Series(hairpins[:, j], dtype=pl.Int32) for j, s in enumerate(samples)})),
        ("mirna_hairpins", loci),
    )  # fmt: skip
    for what, table in tables:
        _write(root, "gather", what, table, day)


def _zip(path: Path, name: str, text: str) -> None:
    """A zip archive of one text file, as TargetScan serves its tables."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(name, text)


def _tsv(header: list[str], rows: list[list[object]]) -> str:
    """Tab-separated text with a header."""
    return "\n".join("\t".join(str(v) for v in r) for r in [header, *rows]) + "\n"


def target_files(root: Path) -> dict[str, Path]:
    """TargetScan's three tables and miRTarBase's strong evidence, made up and small."""
    folder = root / TARGETS
    header = ["Gene ID", "Gene Symbol", "Transcript ID", "Gene Tax ID", "miRNA", "Site Type",
              "UTR_start", "UTR end", "context++ score", "context++ score percentile",
              "weighted context++ score", "weighted context++ score percentile",
              "Predicted relative KD"]  # fmt: skip

    def site(gene: str, mirna: str, score: float, tax: int = 9606) -> list[object]:
        """One predicted site."""
        return [f"{gene}.1", gene.upper(), "ENST0.1", tax, mirna, 3, 1, 8, score, 90, score,
                90, "NULL"]  # fmt: skip

    sites = [site(TARGET, MIR_UP, -0.5), site(PREDICTED, MIR_UP, -0.15),
             site(PREDICTED, MIR_UP, -0.15), site(WEAK, MIR_UP, -0.05),
             site("g00014", "mml-miR-1", -0.9, tax=9544)]  # fmt: skip
    sites += [site(f"g{i:05d}", MIR_UP, -0.3) for i in range(20, 60)]
    sites += [site(f"g{i:05d}", MIR_DOWN, -0.3) for i in range(60, 100)]
    paths = {"targetscan_predictions": folder / "predictions.zip",
             "targetscan_families": folder / "families.zip",
             "targetscan_genes": folder / "genes.zip",
             "mirtarbase_strong": folder / "mirtarbase.csv"}  # fmt: skip
    _zip(paths["targetscan_predictions"], "predictions.txt", _tsv(header, sites))
    families = [
        ["miR-1", "GGAAUGU", 9606, MIR_UP, "U", 2, MIMATS[0]],
        ["miR-2", "AUCACAG", 9606, MIR_DOWN, "U", 2, MIMATS[1]],
        ["miR-1", "GGAAUGU", 10090, "mmu-miR-1a-3p", "U", 2, "MIMAT0000123"],
    ]
    _zip(paths["targetscan_families"], "families.txt", _tsv(
        ["miR family", "Seed+m8", "Species ID", "MiRBase ID", "Mature sequence",
         "Family Conservation?", "MiRBase Accession"], families))  # fmt: skip
    genes = [[f"ENST{i}.1", f"g{i:05d}.1", f"G{i:05d}", "made up", 9606, 1, 1]
             for i in range(1, 201)]  # fmt: skip
    _zip(paths["targetscan_genes"], "genes.txt", _tsv(
        ["Transcript ID", "Gene ID", "Gene symbol", "Gene description", "Species ID",
         "3P-seq tags", "Representative transcript?"], genes))  # fmt: skip
    head = (
        "miRTarBase ID,miRNA,Species (miRNA),Target Gene,Target Gene (Entrez ID),"
        "Species (Target Gene),Experiments,Support Type,References (PMID)"
    )
    lines = [f"MIRT1,{MIR_UP},hsa,G00010,1.0,hsa,Luciferase reporter assay//Western blot,"
             "Functional MTI,1.0", f"MIRT1,{MIR_UP},hsa,G00010,1.0,hsa,qRT-PCR,Functional MTI,2.0",
             f"MIRT2,{MIR_UP},hsa,G00002,2.0,hsa,Western blot,Functional MTI,3.0",
             "MIRT3,mmu-miR-1a-3p,mmu,Cd320,3.0,mmu,Western blot,Functional MTI,4.0"]  # fmt: skip
    paths["mirtarbase_strong"].write_text("﻿" + "\n".join([head, *lines]) + "\n",
                                          encoding="utf-8")  # fmt: skip
    return paths


def spec_text(root: Path, paths: dict[str, Path]) -> str:
    """The dataset file: [gdc], the resources with their hashes, the contrast and [ncrna]."""
    resources = ""
    for name, path in {"annotation_gtf": root / GTF, **paths}.items():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        resources += (f'\n[[resources]]\nname = "{name}"\nurls = ["http://127.0.0.1/{path.name}"]'
                      f'\nbytes = {path.stat().st_size}\nsha256 = "{digest}"\n'
                      f'dest = "{path.relative_to(root)}"\n')  # fmt: skip
    contrast = (
        'variable = "condition"\nreference = "normal"\ntested = "tumour"\n'
        f'design = "~donor + group"\npaired = true\nn = [{PAIRS}, {PAIRS}]\n'
    )
    return f"""[dataset]
slug = "{DATASET}"
label = "Made-up GDC"

[gdc]
project = "FAKE-GDC"
release = "Data Release 99.0"
workflow = "STAR - Counts"
gene_model = "GENCODE v36"
protocol = "unstranded"
tumour = "Primary Tumor"
normal = "Solid Tissue Normal"
mirna_workflow = "BCGSC miRNA Profiling"

[guard]
ids = []
patterns = []
{resources}
[[contrasts]]
name = "tumour_normal"
{contrast}
[de]
featured = "tumour_normal"

[ncrna]
featured = "tumour_normal"
pairs = 5
mirnas = ["{MIR_UP}", "{MIR_DOWN}"]

[ncrna.mirna]
name = "mirna_tumour_normal"
{contrast}
[ncrna.sanity]
registered_on = "2026-10-08 00:00"

[[ncrna.sanity.rules]]
id = "up"
text = "The planted miRNA rises."
source = "test"
conditions = [{{ metric = "mirna.{MIR_UP}.lfc", op = ">=", value = 1.0 }}]

[[ncrna.sanity.rules]]
id = "pair"
text = "Its target falls as it rises, within tumours."
source = "test"
conditions = [{{ metric = "pair.{MIR_UP}.G00010.rho_tumour", op = "<", value = 0.0 }}]

[[ncrna.sanity.rules]]
id = "absent"
text = "STAR's miRNA genes disagree with the small-RNA libraries."
source = "test"
conditions = [{{ metric = "agreement.median_rho", op = "<", value = 0.3 }}]

[[ncrna.sanity.rules]]
id = "unmeasured"
text = "A rule on a metric no run measures."
source = "test"
conditions = [{{ metric = "library.histone_pairs", op = ">=", value = 1 }}]
"""


def make(tmp: Path, day: date) -> tuple[Path, Path]:
    """Write the root and the dataset file under ``tmp``.

    Returns:
        The root and the dataset file.
    """
    root = tmp / "root"
    sheet, frame, names = star(day, root)
    mirna(day, root, sheet)
    write_gtf(root / GTF, frame["gene_id"].to_list(), names)
    taken = sheet.filter(pl.col("chosen"))
    found = ncrna_root.results(frame, taken, names, tested="tumour").with_columns(
        call=classify(pl.col("padj"), pl.col("log2FoldChange"))
    )
    _write(root, "de", "tumour_normal_results", found, day)
    spec = tmp / f"{DATASET}.toml"
    spec.write_text(spec_text(root, target_files(root)))
    return root, spec
