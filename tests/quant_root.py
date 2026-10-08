"""A ctx root for `ctx quant` on the quant fixture: the chr22_se samples, gathered.

The five samples' gather outputs come from ``process_root``, their cell type set to the
expression profile and their batch to the read length each was simulated with
(``tests/fixtures/chr22_se/samples.tsv``). The dataset file declares the chr22 slice's GTF with
two multi-isoform genes typed lncRNA, so the isoform figure has two types to compare, and one
contrast, case over control, as its ``[de] featured``. The fixture's outdir is read in place.
"""

import hashlib
from pathlib import Path

import polars as pl
from process_root import DATASET, SE, annotation, plant

ROOT = Path(__file__).resolve().parents[1]
QUANT = ROOT / "tests" / "fixtures" / "quant"
RNASEQ = QUANT / "rnaseq"
SALMON = RNASEQ / "salmon"
CONTRAST = "case_control"
LNCRNA = ("ENSG00000099995.21", "ENSG00000100003.18")  # SF3A1 and SEC14L2, retyped here
BATCHES = {76: "short", 101: "long"}


def _sheet(root: Path) -> None:
    """Each sample's cell type from its profile and batch from its read length."""
    made = pl.read_csv(SE / "samples.tsv", separator="\t")
    cells = dict(zip(made["sample"], (f"type_{p}" for p in made["profile"]), strict=True))
    batches = dict(zip(made["sample"], (BATCHES[n] for n in made["read_length"]), strict=True))
    path = next((root / "outputs" / "01_gather").glob(f"*_{DATASET}_samplesheet_v1.parquet"))
    pl.read_parquet(path).with_columns(
        cell_type=pl.col("sample").replace_strict(cells),
        batch=pl.col("sample").replace_strict(batches)).write_parquet(path)  # fmt: skip


def quant_root(root: Path) -> Path:
    """Plant the gather outputs and write the dataset file; return the dataset file."""
    plant(root)
    _sheet(root)
    gtf = annotation(root)
    lines = []
    for line in gtf.read_text().splitlines(keepends=True):
        retyped = any(f'gene_id "{g}"' in line for g in LNCRNA)
        lines.append(line.replace('gene_type "protein_coding"', 'gene_type "lncRNA"')
                     if retyped else line)  # fmt: skip
    gtf.write_text("".join(lines))
    digest = hashlib.sha256(gtf.read_bytes()).hexdigest()
    spec = root / f"{DATASET}.toml"
    spec.write_text(
        f'[dataset]\nslug = "{DATASET}"\n\n[guard]\nids = []\npatterns = []\n\n[[resources]]\n'
        f'name = "annotation_gtf"\nurls = ["http://127.0.0.1/slice_plus.gtf"]\n'
        f'bytes = {gtf.stat().st_size}\nsha256 = "{digest}"\ndest = "{gtf.relative_to(root)}"\n'
        f'\n[[contrasts]]\nname = "{CONTRAST}"\nvariable = "condition"\nreference = "control"\n'
        f'tested = "case"\nn = [3, 2]\n\n[de]\nfeatured = "{CONTRAST}"\n'
    )  # fmt: skip
    return spec


__all__ = ["DATASET", "RNASEQ", "SALMON", "quant_root"]
