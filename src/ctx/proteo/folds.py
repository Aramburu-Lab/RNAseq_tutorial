"""Protein fold changes and their calls, and F4: one contrast's RNA and protein calls compared.

``protein_calls`` is the one protein test of M7: Welch's t on present values, BH over the
proteins quantified in enough samples, then D-012's rule through ``ctx.de.stats.classify``. The
cross-study design and F4 both call it. F4's RNA calls are M5's DESeq2 results for the same
contrast, called by the same ``classify``; its protein samples are the tumours M5 compared, taken
by M5's own ``select_samples``, so both layers are tested on the same tumours.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from ctx import concordance
from ctx.dataset import Contrast, DatasetError
from ctx.de.deseq import select_samples
from ctx.de.run import newest
from ctx.de.stats import adjust, classify
from ctx.proteo import ProteoError, measure

# Where a gene was called, as F4 groups it: at both layers, at one, or at neither.
KINDS = ("both", "rna", "protein", "neither")


def protein_calls(genes: Sequence[str], values: np.ndarray, tested: np.ndarray,
                  minimum: float, lfc_min: float) -> pl.DataFrame:  # fmt: skip
    """Per protein: the share quantified, Welch's fold change and p, BH, and the call.

    Order: the quantified filter in exact fractions over every sample of the comparison, then
    Welch's t, then BH over the filtered proteins with a p, then ``classify``.

    Args:
        genes: One symbol per row of ``values``.
        values: Genes by samples, log2, NaN where missing.
        tested: True for each column in the tested group, False for the reference.
        minimum: Share of the comparison's samples a protein needs a value in to be tested.
        lfc_min: Absolute log2 fold change a call needs (D-012).

    Returns:
        One row per gene: gene, n_samples, quantified_share, tested, protein_lfc, protein_p,
        protein_padj and protein_call.
    """
    present = (~np.isnan(values)).sum(axis=1)
    kept = measure.enough(present, values.shape[1], minimum)
    lfc, pvalue = measure.welch(values, tested)
    padj = np.full(lfc.size, np.nan)
    ok = kept & ~np.isnan(pvalue)
    if ok.any():
        padj[ok] = adjust(pvalue[ok], "padj")
    frame = pl.DataFrame({
        "gene": list(genes), "n_samples": present.astype(np.int64),
        "quantified_share": present / values.shape[1], "tested": kept,
        "protein_lfc": lfc, "protein_p": pvalue, "protein_padj": padj,
    }).with_columns(pl.col("protein_p", "protein_padj").fill_nan(None))  # fmt: skip
    return frame.with_columns(
        protein_call=classify(pl.col("protein_padj"), pl.col("protein_lfc"), lfc_min=lfc_min)
    )


def kind(rna_call: pl.Expr, protein_call: pl.Expr) -> pl.Expr:
    """Where a gene was called: at both layers, RNA only, protein only, or neither."""
    rna, protein = rna_call != "ns", protein_call != "ns"
    return (pl.when(rna & protein).then(pl.lit(KINDS[0])).when(rna).then(pl.lit(KINDS[1]))
            .when(protein).then(pl.lit(KINDS[2])).otherwise(pl.lit(KINDS[3])))  # fmt: skip


@dataclass(frozen=True)
class Folds:
    """F4 for one contrast.

    Attributes:
        contrast: The M5 contrast.
        genes: One row per gene tested at both layers: both fold changes, adjusted p, calls,
            and ``kind``.
        n_protein: (reference, tested) tumours with a protein value table column.
        rho: Spearman rho of RNA against protein fold change, and its bootstrap interval.
    """

    contrast: Contrast
    genes: pl.DataFrame
    n_protein: tuple[int, int]
    rho: tuple[float | None, float | None, float | None]

    def counts(self) -> dict[str, int]:
        """Calls per layer and direction, and genes called at both layers in opposite ways."""
        g = self.genes
        out = {f"{layer}_{d}": int((g[f"{layer}_call"] == d).sum())
               for layer in ("rna", "protein") for d in ("up", "down")}  # fmt: skip
        out |= {k: int((g["kind"] == k).sum()) for k in KINDS}
        out["opposite"] = int(((g["kind"] == KINDS[0])
                               & (g["rna_call"] != g["protein_call"])).sum())  # fmt: skip
        return out

    def metrics(self) -> list[tuple[str, float, str]]:
        """(metric, value, detail) rows for the metrics table."""
        name, detail = self.contrast.name, f"{self.genes.height} genes tested at both layers"
        rows = [(f"f4.{name}.{k}", float(v), detail) for k, v in self.counts().items()]
        rows += [(f"f4.{name}.spearman{k}", np.nan if v is None else v,
                  f"{detail}, bootstrap 1000, seed 42")
                 for k, v in zip(("", ".low", ".high"), self.rho, strict=True)]  # fmt: skip
        return rows


def rna_calls(results: pl.DataFrame, lfc_min: float) -> pl.DataFrame:
    """M5's tested genes, one row per symbol, called by D-012's rule."""
    tested = results.filter(pl.col("padj").is_not_null() & pl.col("gene_name").is_not_null())
    shared = tested["gene_name"].is_duplicated()
    return tested.filter(~shared).select(
        gene="gene_name", rna_lfc="log2FoldChange", rna_padj="padj",
        rna_call=classify(pl.col("padj"), pl.col("log2FoldChange"), lfc_min=lfc_min),
    )  # fmt: skip


def contrast_folds(  # noqa: PLR0913 - each a separate input of one comparison
    results: pl.DataFrame, sheet: pl.DataFrame, pairs: Sequence[tuple[int, str]],
    protein: tuple[Sequence[str], np.ndarray], contrast: Contrast, *,
    order: dict[str, tuple[str, ...]], limits: tuple[float, float],
) -> Folds:  # fmt: skip
    """F4 for one contrast on the tumours M5 compared that have a protein column.

    Args:
        results: M5's results table for the contrast, with ``gene_name``.
        sheet: Gather's sheet as M5 reads it (``gdc.chosen``), with the contrast's variable.
        pairs: (protein column, sheet sample) of every paired tumour.
        protein: The protein genes and their genes-by-columns values.
        contrast: The contrast.
        order: The dataset's level orders, as M5 passes them.
        limits: (minimum quantified share, lfc_min).

    Raises:
        ProteoError: If the sheet lacks the contrast's variable, or a group has no protein.
    """
    if contrast.variable not in sheet.columns:
        raise ProteoError(f"{contrast.name}: gather's sheet has no column {contrast.variable!r}")
    used, _ = select_samples(sheet, contrast, order)
    group = dict(used.select("sample", "group").rows())
    picked = [(j, group[s] == "tested") for j, s in pairs if s in group]
    columns = np.array([j for j, _ in picked], dtype=np.int64)
    is_tested = np.array([t for _, t in picked], dtype=bool)
    n = (int((~is_tested).sum()), int(is_tested.sum()))
    if min(n) < 2:  # noqa: PLR2004 - Welch's t needs two values per group
        raise ProteoError(f"{contrast.name}: protein for {n[0]} and {n[1]} tumours")
    genes, values = protein
    layer = protein_calls(genes, values[:, columns], is_tested, *limits)
    joined = rna_calls(results, limits[1]).join(
        layer.filter(pl.col("tested") & pl.col("protein_lfc").is_not_nan()),
        on="gene", how="inner", maintain_order="right",
    ).with_columns(kind=kind(pl.col("rna_call"), pl.col("protein_call")))  # fmt: skip
    rho, low, high, _ = concordance.spearman(joined["rna_lfc"], joined["protein_lfc"])
    return Folds(contrast, joined, n, (rho, low, high))


def results_file(de_dir: Path, slug: str, contrast: Contrast) -> Path | None:
    """M5's newest results file for the contrast, or None when ``ctx de`` has not written it."""
    try:
        return newest(de_dir, slug, f"{contrast.name}_results")
    except DatasetError:
        return None
