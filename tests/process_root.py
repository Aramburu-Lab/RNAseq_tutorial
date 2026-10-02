"""A ctx root for `ctx process` on the chr22_se fixture: gather outputs, FASTQs, annotation.

The five single-end samples become a gathered dataset: the sheet names their gzipped FASTQ
files in a fetchngs-shaped outdir, and a dataset file declares the fixture slice's GTF as its
`annotation_gtf`, with the SHA-256 the run checks. A second dataset can stand beside it with
Salmon-like counts and ctx qc metrics for the same samples, for `--compare-dataset`.
"""

import gzip
import hashlib
import shutil
from datetime import date
from pathlib import Path

import polars as pl
from synth import Synthetic, design_counts, write_gather_outputs

ROOT = Path(__file__).resolve().parents[1]
SE = ROOT / "tests" / "fixtures" / "chr22_se"
SLICE = ROOT / "tests" / "fixtures" / "chr22" / "slice.fa"
GTF = ROOT / "tests" / "fixtures" / "chr22" / "slice.gtf"
SAMPLES = ("se_a", "se_b", "se_c", "se_d", "se_e")
DATASET, COMPARE = "chrse", "chrse_full"
# The parent a subset's dataset file names (D-042); its page says it is the laptop exercise.
PARENT = "chrse_cohort"
DAY = date(2026, 10, 1)
# A gene on another contig, kept in the compare counts so the expected share is below 1.
OFF_CONTIG = ("ENSG00000000001.1", "elsewhere")


def annotation(root: Path) -> Path:
    """The slice's GTF plus one gene on another contig, under the root's resources."""
    gtf = root / "resources" / "gencode" / "slice_plus.gtf"
    gtf.parent.mkdir(parents=True, exist_ok=True)
    gene, contig = OFF_CONTIG
    extra = (f'{contig}\tctx\tgene\t1\t900\t.\t+\t.\tgene_id "{gene}"; gene_type '
             f'"protein_coding"; gene_name "ELSEWHERE";\n')  # fmt: skip
    gtf.write_text(GTF.read_text() + extra)
    return gtf


def dataset_file(root: Path, folder: Path, slug: str = DATASET, *, subset_of: str = "") -> Path:
    """A dataset file whose annotation_gtf is the root's GTF, with its real size and SHA-256."""
    gtf = annotation(root)
    digest = hashlib.sha256(gtf.read_bytes()).hexdigest()
    head = f'[dataset]\nslug = "{slug}"\n' + (f'subset_of = "{subset_of}"\n' if subset_of else "")
    spec = folder / f"{slug}.toml"
    spec.write_text(
        head + '\n[guard]\nids = []\npatterns = []\n\n[[resources]]\nname = "annotation_gtf"\n'
        f'urls = ["http://127.0.0.1/slice_plus.gtf"]\nbytes = {gtf.stat().st_size}\n'
        f'sha256 = "{digest}"\ndest = "{gtf.relative_to(root)}"\n'
    )  # fmt: skip
    return spec


def plant(root: Path, samples: tuple[str, ...] = SAMPLES, slug: str = DATASET) -> pl.DataFrame:
    """Gather outputs naming the fixture's reads, gzipped into ``data/raw/<slug>/fastq``."""
    fastq = root / "data" / "raw" / slug / "fastq"
    fastq.mkdir(parents=True, exist_ok=True)
    for name in samples:
        with (SE / f"{name}.fq").open("rb") as src, gzip.open(fastq / f"{name}.fastq.gz",
                                                              "wb") as dst:  # fmt: skip
            shutil.copyfileobj(src, dst)
    synthetic = design_counts((("control", "short", 3), ("case", "short", 2)))
    names = dict(zip(synthetic.samples["sample"], samples, strict=False))
    kept = synthetic.samples.filter(pl.col("sample").is_in(list(names)))
    renamed = Synthetic(synthetic.counts.select("gene_id", *names).rename(names),
                        kept.with_columns(pl.col("sample").replace_strict(names)),
                        synthetic.truth)  # fmt: skip
    write_gather_outputs(root, slug, renamed, DAY)
    path = next((root / "outputs" / "01_gather").glob(f"*_{slug}_samplesheet_v1.parquet"))
    sheet = pl.read_parquet(path).with_columns(
        fastq_1=pl.format("fastq/{}.fastq.gz", "sample"), status_1=pl.lit("verified"),
        fastq_2=pl.lit(""), status_2=pl.lit(""))  # fmt: skip
    sheet.write_parquet(path)
    return sheet


def plant_compare(root: Path, slug: str = COMPARE, samples: tuple[str, ...] = SAMPLES,
                  off: dict[str, int] | None = None) -> None:  # fmt: skip
    """Salmon-like counts and qc metrics of the same samples under another dataset's name.

    Each sample holds 1,000 counts: 100 on the off-contig gene, or ``off``'s number, and the
    rest spread evenly over the slice's genes, so the expected share is 0.9 by default; Salmon
    maps 80 %, and se_e alone is unstranded.
    """
    genes = [line.split('gene_id "')[1].split('"')[0] for line in GTF.read_text().splitlines()
             if "\tgene\t" in line]  # fmt: skip
    elsewhere = {s: (off or {}).get(s, 100) for s in samples}
    counts = pl.DataFrame({"gene_id": [*genes, OFF_CONTIG[0]],
                           **{s: [(1000 - n) // len(genes)] * len(genes) + [n]
                              for s, n in elsewhere.items()}})  # fmt: skip
    gathered, qc = root / "outputs" / "01_gather", root / "outputs" / "02_qc"
    gathered.mkdir(parents=True, exist_ok=True)
    qc.mkdir(parents=True, exist_ok=True)
    counts.write_parquet(gathered / f"{DAY}_gather_{slug}_counts_v1.parquet")
    pl.DataFrame({
        "sample": list(samples), "salmon-percent_mapped": [80.0] * len(samples),
        "salmon-strand_mapping_bias": [0.5 if s == "se_e" else 0.07 for s in samples],
    }).write_parquet(qc / f"{DAY}_qc_{slug}_metrics_v1.parquet")  # fmt: skip
