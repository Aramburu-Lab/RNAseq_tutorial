"""A made-up root for ``ctx proteo`` on a ``[proteo.table]``: gather's outputs written as any
route writes them, a MaxQuant or DIA-NN export, and a sample table, for tests/test_proteo_study.py.

The RNA side is three gather tables (sheet, counts, genes) of a study with no GDC column; the
protein side follows RNA for the planted metabolic genes and G9, not for the complex subunits.
Every row the reader must leave out is planted once: a reverse decoy and a contaminant that
would win G9's and G1's collapse if kept, a group of two genes, a symbol the RNA lacks, a
symbol two RNA genes share. Nothing here comes from a real dataset.
"""

from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
from clinical_root import ctx
from proteo_root import COMPLEX, METABOLIC, RULES, resource_toml, resources

from ctx.de.filters import log2cpm

SAMPLES = 24
# The planted RNA genes: G1 to G12 and a symbol two gene IDs share.
GENES = [f"G{i}" for i in range(1, 13)]
SLUG = "fake_study"


def gather(root: Path, rng: np.random.Generator) -> tuple[list[str], np.ndarray]:
    """Gather's three tables, as any route writes them; returns the samples and log2 CPM."""
    out = root / "outputs" / "01_gather"
    out.mkdir(parents=True)
    samples = [f"rna_{i:02d}" for i in range(SAMPLES)]
    ids = [f"ENSG{i:011d}.{i % 3 + 1}" for i in range(len(GENES) + 2)]
    names = [*GENES, "DUP", "DUP"]
    means = rng.uniform(4, 10, len(names))
    counts = rng.poisson(np.exp(means[:, None] + rng.normal(0, 0.8, (len(names), SAMPLES))))
    day = date.today().isoformat()
    pl.DataFrame({"sample": samples, "condition": ["a", "b"] * (SAMPLES // 2)}).write_parquet(
        out / f"{day}_gather_{SLUG}_samplesheet_v1.parquet"
    )
    pl.DataFrame(
        {"gene_id": ids, **{s: counts[:, j] for j, s in enumerate(samples)}}
    ).write_parquet(out / f"{day}_gather_{SLUG}_counts_v1.parquet")
    pl.DataFrame({"gene_id": ids, "gene_name": names}).write_parquet(
        out / f"{day}_gather_{SLUG}_genes_v1.parquet"
    )
    return samples, log2cpm(counts[: len(GENES)])


def protein(rna: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """log2 protein: metabolic genes and G9 follow RNA, the rest do not; G10 missing in half."""
    follows = np.array([g in METABOLIC or g == "G9" for g in GENES])
    values = 20 + np.where(follows[:, None], rna - rna.mean(axis=1, keepdims=True), 0.0)
    values += rng.normal(0, np.where(follows, 0.02, 1.0)[:, None], rna.shape)
    values[GENES.index("G10"), ::2] = np.nan
    values[GENES.index("G9"), 0] = np.nan  # the real G9 row is one sample short of the decoy's
    return values


def maxquant_text(values: np.ndarray, columns: list[str], rng: np.random.Generator) -> str:
    """proteinGroups.txt: LFQ columns, 0 for missing, the planted rows each flagged once."""
    head = [
        "Protein IDs",
        "Majority protein IDs",
        "Gene names",
        "Intensity",
        *[f"LFQ intensity {c}" for c in columns],
        "Only identified by site",
        "Reverse",
        "Potential contaminant",
        "id",
    ]
    noise = 2 ** rng.normal(20, 2, len(columns))

    def row(ids: str, gene: str, cells: np.ndarray, flags: tuple[str, str, str]) -> str:
        """One line, NaN written as MaxQuant's 0."""
        lfq = ["0" if np.isnan(v) else f"{2**v:.1f}" for v in cells]
        return "\t".join([ids, ids, gene, "1", *lfq, *flags, "0"])

    lines = ["\t".join(head)]
    lines += [row(f"P{i}", g, values[i], ("", "", "")) for i, g in enumerate(GENES)]
    lines += [
        row("REV__P9", "G9", np.log2(noise), ("", "+", "")),
        row("CON__P1", "G1", np.log2(noise), ("", "", "+")),
        row("P90", "G2", np.log2(noise), ("+", "+", "")),  # site and reverse: once
        row("P91;P92", "G2;G3", np.log2(noise), ("", "", "")),
        row("P93", "NOTINRNA", np.log2(noise), ("", "", "")),
        row("P94", "DUP", np.log2(noise), ("", "", "")),
        row("P95", "", np.log2(noise), ("", "", "")),
    ]
    return "\n".join(lines) + "\n"


def diann_text(values: np.ndarray, columns: list[str]) -> str:
    """pg_matrix.tsv: each run named by a Windows path, missing empty, no decoy rows."""
    runs = [f"D:\\runs\\{c}.raw" for c in columns]
    head = [
        "Protein.Group",
        "Protein.Ids",
        "Protein.Names",
        "Genes",
        "First.Protein.Description",
        *runs,
    ]
    lines = ["\t".join(head)]
    for i, g in enumerate(GENES):
        cells = ["" if np.isnan(v) else f"{2**v:.1f}" for v in values[i]]
        lines.append("\t".join([f"P{i}", f"P{i}", f"N{i}", g, "made up", *cells]))
    return "\n".join(lines) + "\n"


def build(tmp: Path, fmt: str) -> tuple[Path, Path]:
    """The root and dataset file; ctx proteo not yet run."""
    rng = np.random.default_rng(11)
    root = tmp / "root"
    samples, rna = gather(root, rng)
    values = protein(rna, rng)
    columns = [f"ms_{j:02d}" for j in range(SAMPLES)] + ["ms_extra"]
    values = np.hstack([values, values[:, :1]])
    body = maxquant_text(values, columns, rng) if fmt == "maxquant" else diann_text(values, columns)
    name = "proteinGroups.txt" if fmt == "maxquant" else "report.pg_matrix.tsv"
    (tmp / name).write_text(body)
    # The extra protein sample has no RNA: listed, drawn in the missingness figures, not paired.
    rows = [
        f"{c},{samples[j] if j < SAMPLES else ''},plex{j % 3},{('a', 'b')[j % 2]}"
        for j, c in enumerate(columns)
    ]
    (tmp / "fake_study_proteo.csv").write_text(
        "protein_sample,rna_sample,plex,condition\n" + "\n".join(rows) + "\n"
    )
    resources(root)
    gmt = root / "resources" / "msigdb" / "classes.gmt"
    gmt.write_text(
        f"SET_COMPLEX\tx\t{chr(9).join(COMPLEX)}\nSET_METABOLIC\tx\t{chr(9).join(METABOLIC)}\n"
    )
    spec = tmp / f"{SLUG}.toml"
    spec.write_text(
        f'[dataset]\nslug = "{SLUG}"\n[guard]\nids = []\npatterns = []\n'
        + resource_toml("classes", "resources/msigdb/classes.gmt")
        + '\n[proteo]\ndesign = "across_samples"\ngenes = ["G9"]\ntracks = ["condition"]\n'
        + '[proteo.classes.complex]\nresource = "classes"\nsets = ["SET_COMPLEX"]\n'
        + '[proteo.classes.metabolic]\nresource = "classes"\nsets = ["SET_METABOLIC"]\n'
        + f'[proteo.table]\nformat = "{fmt}"\nfile = "{name}"\nsamples = "fake_study_proteo.csv"\n'
        + 'scale = "linear"\nid_kind = "symbol"\nacquisition = "tmt"\nbatch = "plex"\n'
        + RULES.format(extra="")
    )
    return root, spec


def run(tmp: Path, fmt: str) -> Path:
    """``build``, then ctx proteo on it."""
    root, spec = build(tmp, fmt)
    ctx("--root", str(root), "proteo", "--dataset", SLUG, "--dataset-file", str(spec))
    return root
