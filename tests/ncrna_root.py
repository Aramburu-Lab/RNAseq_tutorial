"""A made-up ctx root for ``ctx ncrna``: gather's tables, M5's results, an annotation, a file.

Used by tests/test_ncrna.py. Counts are tests/synth.py's design_counts; the annotation gives
each gene a type and a place, so every pair kind occurs, and one lncRNA, ``g00201`` (gene
200, as synth numbers genes from 1), is planted to follow its protein-coding neighbour
``g00001`` sample by sample. M5's results are made here
by M5's own filter and a Welch t-test, not by a fit: ``ctx ncrna`` reads only their columns.
Nothing comes from a real dataset.
"""

import hashlib
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control, ttest_ind
from synth import Synthetic, design_counts, write_gather_outputs

from ctx.de.filters import MIN_COUNT, expressed
from ctx.io.naming import stem
from ctx.io.readme_companion import write_table

DATASET = "ncsynth"
CONTRAST = "a_case"
GROUPS = (("control", "short", 6, "type_a"), ("case", "short", 6, "type_a"))
N_GENES = 400
PLANTED = ("g00201", "g00001")
# Gene index ranges and their types; the last six genes are synth's sex genes.
TYPES = ((0, 200, "protein_coding"), (200, 320, "lncRNA"), (320, 350, "unprocessed_pseudogene"),
         (350, 352, "Mt_rRNA"), (352, 372, "miRNA"), (372, 392, "snRNA"),
         (392, N_GENES, "TEC"))  # fmt: skip
PC_SPACING, PC_LENGTH, LNC_LENGTH = 100_000, 5_000, 2_000
GTF = "resources/annotation.gtf"


def gene_type(i: int) -> str:
    """The made-up type of gene ``i``."""
    return next(t for low, high, t in TYPES if low <= i < high)


def locus(i: int) -> tuple[str, int, int, str]:
    """chrom, start, end, strand of gene ``i``: lncRNA ``200 + j`` sits by protein-coding gene j.

    By j mod 4: overlapping it on the other strand, overlapping on its strand, 5,000 nt past
    its end, or 50,000 nt past it, which leaves the next gene nearer but still far.
    """
    kind = gene_type(i)
    if kind == "protein_coding":
        start = PC_SPACING * i + 1
        # Genes 150 on lie on chr3, beside no lncRNA: the background's other chromosome.
        chrom = "chr1" if i < 150 else "chr3"  # noqa: PLR2004 - see the comment
        return chrom, start, start + PC_LENGTH - 1, "+" if i % 2 == 0 else "-"
    if kind == "Mt_rRNA":
        return "chrM", 1000 * i, 1000 * i + 900, "+"
    if kind != "lncRNA":
        return "chr2", 10_000 * i + 1, 10_000 * i + 500, "+"
    j = i - 200
    pc_start, pc_end, pc_strand = PC_SPACING * j + 1, PC_SPACING * j + PC_LENGTH, "+-"[j % 2]
    other = "-" if pc_strand == "+" else "+"
    placed = {0: (pc_start + 1000, other), 1: (pc_start + 1000, pc_strand),
              2: (pc_end + 5001, pc_strand), 3: (pc_end + 50_001, pc_strand)}  # fmt: skip
    start, strand = placed[j % 4]
    return "chr1", start, start + LNC_LENGTH - 1, strand


def write_gtf(path: Path, ids: list[str], names: list[str]) -> None:
    """A GTF of one gene line and one exon line per gene."""
    lines = ["##description: made-up annotation for tests/test_ncrna.py\n"]
    for i, (gene, name) in enumerate(zip(ids, names, strict=True)):
        chrom, start, end, strand = locus(i)
        attributes = f'gene_id "{gene}"; gene_type "{gene_type(i)}"; gene_name "{name}";'
        for feature in ("gene", "exon"):
            lines.append(f"{chrom}\tHAVANA\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t"
                         f"{attributes}\n")  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))


def counts() -> tuple[pl.DataFrame, Synthetic]:
    """design_counts with the planted lncRNA following its neighbour and near-silent miRNAs."""
    synthetic = design_counts(GROUPS, n_genes=N_GENES, seed=11)
    frame = synthetic.counts
    lead = frame.filter(pl.col("gene_id") == PLANTED[1])
    rng = np.random.default_rng(3)
    planted = {c: int(lead[c][0] // 2 + rng.integers(0, 3)) for c in frame.columns[1:]}
    mirna = [f"g{i + 1:05d}" for i in range(352, 372)]
    return frame.with_columns(
        pl.when(pl.col("gene_id") == PLANTED[0]).then(pl.lit(planted[c]))
        .when(pl.col("gene_id").is_in(mirna)).then(pl.lit(1 if j % 5 == 0 else 0))
        .otherwise(pl.col(c)).cast(pl.Int32).alias(c)
        for j, c in enumerate(frame.columns[1:])
    ), synthetic  # fmt: skip


def results(frame: pl.DataFrame, sheet: pl.DataFrame, names: list[str]) -> pl.DataFrame:
    """M5's columns for the genes M5's filter keeps: Welch t on log2 CPM, BH over them."""
    samples = sheet["sample"].to_list()
    case = (sheet["condition"] == "case").to_numpy()
    matrix = frame.select(samples).to_numpy().astype(float)
    keep = expressed(matrix, [int((~case).sum()), int(case.sum())], MIN_COUNT)
    logs = np.log2(matrix / matrix.sum(axis=0) * 1e6 + 1)[keep]
    p = ttest_ind(logs[:, case], logs[:, ~case], axis=1, equal_var=False).pvalue
    return pl.DataFrame({
        "gene_id": frame["gene_id"].filter(pl.Series(keep)),
        "baseMean": matrix[keep].mean(axis=1),
        "log2FoldChange": logs[:, case].mean(axis=1) - logs[:, ~case].mean(axis=1),
        "pvalue": p, "padj": false_discovery_control(p, method="bh"),
        "gene_name": [n for n, k in zip(names, keep, strict=True) if k],
    })  # fmt: skip


def spec_text(gtf: Path) -> str:
    """The dataset file: the contrast, the annotation with its hash, and [ncrna]."""
    digest = hashlib.sha256(gtf.read_bytes()).hexdigest()
    return f"""[dataset]
slug = "{DATASET}"

[[resources]]
name = "annotation_gtf"
urls = ["http://127.0.0.1/annotation.gtf"]
bytes = {gtf.stat().st_size}
sha256 = "{digest}"
dest = "{GTF}"

[guard]
ids = []
patterns = []

[[contrasts]]
name = "{CONTRAST}"
variable = "condition"
reference = "control"
tested = "case"
n = [6, 6]
filters = {{ cell_type = ["type_a"], batch = ["short"] }}

[ncrna]
featured = "{CONTRAST}"
pairs = 5

[ncrna.sanity]
registered_on = "2026-10-08 00:00"

[[ncrna.sanity.rules]]
id = "N1"
text = "miRNA genes hold under 0.1 % of counted reads in every poly-A library."
source = "test"
conditions = [{{ metric = "mirna.max_percent", op = "<", value = 0.1 }}]
"""


def make(tmp: Path, day: date) -> tuple[Path, Path]:
    """Write the root and the dataset file under ``tmp``.

    Returns:
        The root and the dataset file.
    """
    root = tmp / "root"
    frame, synthetic = counts()
    sheet = synthetic.samples.with_columns(library_selection=pl.lit("Oligo-dT"))
    synthetic = type(synthetic)(frame, sheet, synthetic.truth)
    write_gather_outputs(root, DATASET, synthetic, day)
    genes = pl.read_parquet(next((root / "outputs" / "01_gather").glob("*_genes_v1.parquet")))
    ids, names = genes["gene_id"].to_list(), genes["gene_name"].to_list()
    write_gtf(root / GTF, ids, names)
    found = results(frame, sheet, names)
    out = root / "outputs" / "05_de"
    write_table(found, out, stem("de", f"{DATASET}_{CONTRAST}_results", 1, day=day),
                columns={c: ("none", c) for c in found.columns}, data_releases={})  # fmt: skip
    spec = tmp / f"{DATASET}.toml"
    spec.write_text(spec_text(root / GTF))
    return root, spec
