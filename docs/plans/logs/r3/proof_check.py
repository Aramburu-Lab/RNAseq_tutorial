"""R3-1, registered in docs/plans/2026-10-18_r3_pxd038437_registered_check.md before this ran.

The protein table through ctx's generic reader, read as the proof's dataset file declares it,
against the paper's Table S2 RNA on its use_for_correlation columns, matched by TMT channel;
per-gene Pearson on genes with at least ten paired values; the median beside 0.29.
Run through env.sh: python docs/plans/logs/r3/proof_check.py <root> <dataset file>
The supplement is the dataset file's supplement_table_s2, read from the root, where ctx setup
placed it (R3.4); no path outside the root is read for it.
"""

import sys
from pathlib import Path

import numpy as np
import openpyxl
import polars as pl

from ctx import concordance
from ctx.dataset import load
from ctx.proteo import quantities
from ctx.proteo.identifiers import rna_names
from ctx.proteo.study import id_map, protein_path, text_of

PUBLISHED, BAND, MIN_PAIRS = 0.29, 0.05, 10


class Root:
    """The settings the reader's helpers take: the root alone."""

    def __init__(self, root: Path) -> None:
        """Hold the root."""
        self.root = root


root, spec_path = (Path(a) for a in sys.argv[1:3])
spec = load(spec_path)
supplement = root / spec.resource("supplement_table_s2").dest
p = spec.proteo
assert p is not None and p.table is not None
book = openpyxl.load_workbook(supplement, read_only=True)
rows = book["Table S2"].iter_rows(values_only=True)
head = []
for row in rows:
    head.append(row)
    if row[0] == "Gene_id":
        break
use = next(r for r in head if r[4] == "use_for_correlation")
ids = head[-1]
keep = [j for j in range(5, len(ids)) if str(use[j]).upper() == "TRUE"]
genes, values = [], []
for row in rows:
    genes.append((str(row[0]), str(row[1])))
    values.append([np.nan if row[j] in (None, "NA") else float(row[j]) for j in keep])
rna_genes = pl.DataFrame(genes, schema=["gene_id", "gene_name"], orient="row")
rna = np.array(values)
rna_tmt = [ids[j] for j in keep]
print(f"Table S2: {len(genes):,} genes, {len(ids) - 5} columns, {len(keep)} use_for_correlation")

problems: list[str] = []
samples = quantities.sample_table(spec_path.parent / p.table.samples, p.table, problems)
assert samples is not None and not problems, problems
settings = Root(root)
path = protein_path(settings, spec)  # type: ignore[arg-type]
text, id_col = text_of(path, p.member, p.table, samples)
read = quantities.read(text, samples, p.table, rna_genes, id_col, id_map(settings, spec, p.table))  # type: ignore[arg-type]
layer = read.layer
print(f"MaxQuant drops {read.dropped}; mapping {read.mapped}; tokens {read.tokens}")
print(f"layer: {len(layer.genes):,} genes by {len(layer.samples)} channels (references out)")
assert read.samples is not None
column = {r: j for j, r in enumerate(read.samples.rna) if r}
both = [t for t in rna_tmt if t in column]
print(f"channels with both layers: {len(both)} of {len(rna_tmt)} use_for_correlation columns")
shown = {n: i for i, n in enumerate(rna_names(rna_genes))}
r_values = []
for g, i in layer.row().items():
    x = rna[shown[g], [rna_tmt.index(t) for t in both]]
    y = layer.values[i, [column[t] for t in both]]
    ok = ~(np.isnan(x) | np.isnan(y))
    if ok.sum() >= MIN_PAIRS and np.ptp(x[ok]) > 0 and np.ptp(y[ok]) > 0:
        r_values.append(float(concordance.correlation(x[ok], y[ok])))
median = float(np.median(r_values))
verdict = "PASS" if abs(median - PUBLISHED) <= BAND else "FAIL"
print(
    f"R3-1: median per-gene Pearson {median:.4f} over {len(r_values):,} genes with at least "
    f"{MIN_PAIRS} pairs; Ramberger 2024: {PUBLISHED}; band {BAND}: {verdict}"
)
