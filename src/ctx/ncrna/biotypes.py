"""M8's one grouping of GENCODE gene types, its colours, and each sample's reads by group.

Six groups, so six colours: Okabe-Ito without its yellow, which a projector loses (1.32:1 on
white). The rRNA classes are QC's (``ctx.qc.checkpoint``): mitochondrial rRNA is a group of its
own, the reads poly-A capture leaves behind most variably; nuclear rRNA and its pseudogenes,
5S and 5.8S copies, pool with the other small RNAs. Order matters: the rRNA and small-RNA
types are decided before the pseudogene rule, which would otherwise take rRNA_pseudogene.
"""

from collections.abc import Sequence

import numpy as np
import polars as pl

from ctx.ncrna import NcrnaError
from ctx.qc.checkpoint import MT_RRNA, RRNA
from ctx.viz.style import DATA_COLOURS, PROJECTOR_UNSAFE

CODING, LNCRNA, PSEUDOGENE, MITO, SMALL, OTHER = (
    "protein_coding", "lncRNA", "pseudogene", "Mt_rRNA", "small_RNA", "other",
)  # fmt: skip
GROUPS = (CODING, LNCRNA, PSEUDOGENE, MITO, SMALL, OTHER)
# The groups M8 tests by default, each a stratum of its own BH.
TESTED = (CODING, LNCRNA)
MIRNA = "miRNA"
# GENCODE 50's small RNA gene types, miRNA hairpins among them, and nuclear rRNA.
SMALL_TYPES = ("miRNA", "snRNA", "snoRNA", "misc_RNA", "scaRNA", "scRNA", "sRNA", "vault_RNA",
               "ribozyme", "Y_RNA", *RRNA)  # fmt: skip
COLOURS = dict(zip(GROUPS, [c for c in DATA_COLOURS if c not in PROJECTOR_UNSAFE], strict=True))
# library_selection values of poly-A capture; Step-1's 210 runs all read Oligo-dT.
POLY_A = ("Oligo-dT", "PolyA")


def group(gene_type: pl.Expr) -> pl.Expr:
    """The group of each gene type, rRNA and small RNA decided before pseudogene."""
    return (
        pl.when(gene_type.is_in(list(MT_RRNA))).then(pl.lit(MITO))
        .when(gene_type.is_in(list(SMALL_TYPES))).then(pl.lit(SMALL))
        .when(gene_type.str.contains("pseudogene")).then(pl.lit(PSEUDOGENE))
        .when(gene_type == CODING).then(pl.lit(CODING))
        .when(gene_type == LNCRNA).then(pl.lit(LNCRNA))
        .otherwise(pl.lit(OTHER))
    )  # fmt: skip


def annotate(gene_ids: pl.Series, types: pl.DataFrame) -> pl.DataFrame:
    """The annotation's record of each gene, in the given order, with its ``biotype`` group.

    Args:
        gene_ids: Genes, such as the rows of the count matrix.
        types: ``ctx.qc.checkpoint.gene_types`` of the annotation the counts were made with.

    Raises:
        NcrnaError: If a gene is not in the annotation; the counts then came from another one,
            and every share would be wrong.
    """
    typed = gene_ids.to_frame("gene_id").join(
        types, on="gene_id", how="left", maintain_order="left"
    )
    missing = typed.filter(pl.col("gene_type").is_null())["gene_id"]
    if missing.len():
        raise NcrnaError(f"{missing.len():,} genes are not in the annotation, for example "
                         f"{missing.head(5).to_list()}")  # fmt: skip
    return typed.with_columns(biotype=group(pl.col("gene_type")))


def parse(value: Sequence[str]) -> tuple[str, ...]:
    """The groups an option names, in ``GROUPS`` order.

    Raises:
        NcrnaError: If one is no group, or none is named.
    """
    if bad := [v for v in value if v not in GROUPS]:
        raise NcrnaError(f"--biotypes {bad}: each one of {GROUPS}")
    if not value:
        raise NcrnaError(f"--biotypes: name at least one of {GROUPS}")
    return tuple(g for g in GROUPS if g in value)


def shares(matrix: np.ndarray, rows: np.ndarray, levels: Sequence[str]) -> np.ndarray:
    """Each level's percent of each sample's counted reads.

    Args:
        matrix: Counts, genes by samples.
        rows: The level of each gene.
        levels: The levels, in the order of the result's rows.

    Returns:
        Levels by samples, in percent; each column sums to 100 over every level the genes hold.
    """
    total = matrix.sum(axis=0)
    return np.stack([100 * matrix[rows == level].sum(axis=0) / total for level in levels])
