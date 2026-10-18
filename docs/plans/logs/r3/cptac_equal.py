"""R3 item 2: could ccRCC's read_tsv be an adapter of the generic reader with every number equal?

Reads the ccRCC table both ways from Rytis's root, read-only, and compares (a) the genes and
values that reach M7's per-gene step, bit for bit, and (b) the protein layer the page's F1 and
the paired fold changes read, which the generic reader limits to genes the RNA side holds.
Run through env.sh: python docs/plans/logs/r3/cptac_equal.py ~/ctx_kidney
"""

import sys
from pathlib import Path

import numpy as np
import polars as pl

from ctx.dataset import ProteoTable, load
from ctx.de.run import newest
from ctx.proteo import adapters, quantities, tables
from ctx.proteo.run import unique_symbols

root = Path(sys.argv[1]).expanduser()
spec = load(Path(__file__).resolve().parents[4] / "conf" / "datasets" / "cptac_ccrcc.toml")
p = spec.proteo
assert p is not None and p.member
archive = root / spec.resource(p.protein).dest
genes = pl.read_parquet(newest(root / "outputs" / "01_gather", spec.slug, "genes"))

old = tables.read_tsv(archive, p.member, p.descriptors)
names = genes["gene_name"].to_list()
single, _ = unique_symbols(names)
rna_row = {n: i for i, n in enumerate(names) if single[i]}
old_joined = [g for g in old.genes if g in rna_row]

data = quantities.raw_bytes(archive, p.member)
header = adapters.header(data, quoted=True)
samples = tuple(header[len(p.descriptors) :])
declared = ProteoTable("table", "x.csv", "log", ("NA", "NaN", ""), "symbol", "tmt", "plex")
sheet = quantities.SampleTable(samples, samples, ("p",) * len(samples), pl.DataFrame(), "")
text = adapters.table(data, "\t", p.descriptors[0], samples)
new = quantities.read(text, sheet, declared, genes, p.descriptors[0]).layer

print(
    f"read_tsv layer: {len(old.genes):,} genes; joined with unique RNA symbols {len(old_joined):,}"
)
print(f"generic layer:  {len(new.genes):,} genes (only genes the RNA side holds)")
same_genes = old_joined == list(new.genes)
print(f"genes reaching per_gene equal, in order: {same_genes}")
if same_genes:
    a = old.values[[old.row()[g] for g in old_joined]]
    b = new.values
    bits = np.array_equal(a.view(np.uint64), b.view(np.uint64))
    print(f"values equal to the bit: {bits}")
print(f"collapse report read_tsv: {old.collapse}")
print(f"collapse report generic:  {new.collapse}")
