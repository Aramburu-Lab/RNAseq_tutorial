"""M8's miRNA half on a GDC dataset: the paired fit, target shifts, pairs and agreement.

Mature miRNA quantities come from the project's miRNA libraries alone (CONVENTIONS.md, Statistics
rules); the STAR counts' miRNA genes are hairpin loci counted from long fragments and serve only
to show that they disagree. A miRNA is paired with mRNA only within one sample barcode, the
join gather made (``rna_sample``), and only within one condition: pooled, the tumour-normal
gap drives rho (plan decision 4).
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

from ctx.dataset import Contrast
from ctx.de.deseq import Fit, fit
from ctx.de.stats import adjust
from ctx.function.lite import ora
from ctx.ncrna import neighbours, targets

CONDITIONS = ("tumour", "normal")
# Lower bound of a sample stratum the pairs are correlated in: rho of fewer is not reported.
MIN_SAMPLES = 10


@dataclass(frozen=True)
class Layers:
    """The samples with both layers, on one sample barcode, and each layer's values.

    Attributes:
        samples: mirna_sample, rna_sample, case, condition; tumours first, then normals.
        mirna: log2 of 1 plus reads per million miRNA-mapped reads, MIMATs by samples.
        mimats: The rows of ``mirna``.
        rna: log2 of 1 plus counts per million counted reads, genes by samples.
        genes: The rows of ``rna``, unversioned Ensembl IDs.
    """

    samples: pl.DataFrame
    mirna: np.ndarray
    mimats: tuple[str, ...]
    rna: np.ndarray
    genes: tuple[str, ...]


def names(mimats: Sequence[str], named: Mapping[str, str]) -> pl.DataFrame:
    """gene_id (the MIMAT) and gene_name (its miRBase name, else the MIMAT) per mature miRNA."""
    return pl.DataFrame({"gene_id": list(mimats),
                         "gene_name": [named.get(m, m) for m in mimats]})  # fmt: skip


def fit_mature(mature: pl.DataFrame, sheet: pl.DataFrame, contrast: Contrast,
               threads: int) -> Fit:  # fmt: skip
    """M5's fit, unchanged, on the mature counts of the chosen miRNA libraries."""
    counts = mature.rename({"mimat": "gene_id"})
    return fit(counts, sheet, contrast, {}, threads=threads)


def layers(sheet: pl.DataFrame, mature: pl.DataFrame, counts: pl.DataFrame) -> Layers:
    """Both layers of every chosen miRNA library whose sample has a chosen STAR library.

    Args:
        sheet: Gather's miRNA sheet, through ``gdc.chosen``, with ``rna_sample``.
        mature: Mature counts, ``mimat`` then one column per miRNA aliquot.
        counts: STAR counts, ``gene_id`` then one column per STAR aliquot.
    """
    both = (sheet.filter(pl.col("rna_sample").is_not_null()
                         & pl.col("rna_sample").is_in(counts.columns[1:])
                         & pl.col("condition").is_in(list(CONDITIONS)))
            .sort(pl.col("condition").replace_strict(list(CONDITIONS), [0, 1]),
                  "sample"))  # fmt: skip
    reads = both["mirna_reads"].to_numpy().astype(float)
    mir = mature.select(both["sample"].to_list()).to_numpy().astype(float)
    rna = counts.select(both["rna_sample"].to_list()).to_numpy().astype(float)
    return Layers(
        samples=both.select(mirna_sample="sample", rna_sample="rna_sample", case="case",
                            condition="condition"),
        mirna=np.log2(mir / reads * 1e6 + 1), mimats=tuple(mature["mimat"].to_list()),
        rna=np.log2(rna / rna.sum(axis=0) * 1e6 + 1), genes=tuple(counts["gene_id"].to_list()),
    )  # fmt: skip


def pair_table(both: Layers, pairs: pl.DataFrame) -> pl.DataFrame:
    """Spearman rho of each miRNA and target across the samples of each condition apart.

    Args:
        both: The two layers.
        pairs: mimat, mirna, gene_id, gene_name and source, one row per pair.

    Returns:
        The pairs with n, rho, p and BH q per condition; BH over the pairs within it. A
        condition with fewer than ``MIN_SAMPLES`` samples gets no rho.
    """
    mir_row = {m: i for i, m in enumerate(both.mimats)}
    gene_row = {g: i for i, g in enumerate(both.genes)}
    kept = pairs.filter(pl.col("mimat").is_in(list(mir_row)) & pl.col("gene_id").is_in(
        list(gene_row)))  # fmt: skip
    a = both.mirna[[mir_row[m] for m in kept["mimat"]]]
    b = both.rna[[gene_row[g] for g in kept["gene_id"]]]
    columns: dict[str, Any] = {}
    for condition in CONDITIONS:
        cols = (both.samples["condition"] == condition).to_numpy()
        n = int(cols.sum())
        rho = np.full(kept.height, np.nan)
        if n >= MIN_SAMPLES and kept.height:
            rho = neighbours.spearman(a[:, cols], b[:, cols])
        p = neighbours.p_value(rho, n) if n > 2 else np.full(kept.height, np.nan)  # noqa: PLR2004 - t needs n - 2 > 0
        q = np.full(kept.height, np.nan)
        held = np.isfinite(p)
        q[held] = adjust(p[held], "padj")
        columns |= {f"n_{condition}": [n] * kept.height, f"rho_{condition}": rho,
                    f"p_{condition}": p, f"q_{condition}": q}  # fmt: skip
    return kept.with_columns(**{k: pl.Series(k, v) for k, v in columns.items()})


def target_pairs(found: targets.Targets, mirnas: Sequence[str], mimat: Mapping[str, str],
                 tested: pl.DataFrame) -> pl.DataFrame:  # fmt: skip
    """The featured miRNAs' validated and predicted-only targets among the tested genes.

    Predicted at the default preset; the two sources stay apart in ``source``.
    """
    parts = []
    for mirna in mirnas:
        if mirna not in mimat:
            continue
        classed = targets.classes(found, mirna, tested, targets.DEFAULT_PRESET)
        parts.append(classed.filter(pl.col("target").is_in([targets.VALIDATED, targets.PREDICTED]))
                     .select(mimat=pl.lit(mimat[mirna]), mirna=pl.lit(mirna), gene_id="gene_id",
                             gene_name="gene_name", source="target"))  # fmt: skip
    schema = {k: pl.String for k in ("mimat", "mirna", "gene_id", "gene_name", "source")}
    return pl.concat(parts) if parts else pl.DataFrame(schema=schema)


def shifts(found: targets.Targets, mirnas: Sequence[str], results: pl.DataFrame) -> pl.DataFrame:
    """Each target class's fold changes against the non-targets', per miRNA and preset.

    The validated class does not depend on the preset and is tested once; the predicted-only
    class once per preset. BH runs over every row, KS and Mann-Whitney each on its own.

    Args:
        found: The target files.
        mirnas: The featured miRNAs.
        results: M5's tested genes: gene_id, gene_name and log2FoldChange.
    """
    rows = []
    for mirna in mirnas:
        for preset in targets.PRESETS:
            classed = targets.classes(found, mirna, results, preset)
            lfc = {c: classed.filter(pl.col("target") == c)["log2FoldChange"].to_numpy()
                   for c in targets.CLASSES}  # fmt: skip
            tested = [targets.PREDICTED]
            tested += [targets.VALIDATED] if preset == targets.DEFAULT_PRESET else []
            for target in tested:
                found_shift = targets.shift(lfc[target], lfc[targets.NON_TARGET])
                rows.append({"mirna": mirna, "target": target,
                             "preset": preset if target == targets.PREDICTED else None,
                             "n_non": lfc[targets.NON_TARGET].size, **found_shift})  # fmt: skip
    frame = pl.DataFrame(rows, schema_overrides={"preset": pl.Float64, "n": pl.Int64})
    for test in ("ks", "mw"):
        p = frame[f"p_{test}"].to_numpy()
        q = np.full(p.size, np.nan)
        held = np.isfinite(p)
        q[held] = adjust(p[held], "padj")
        frame = frame.with_columns(pl.Series(f"q_{test}", q))
    return frame


def curves(found: targets.Targets, mirna: str, results: pl.DataFrame,
           preset: float = targets.DEFAULT_PRESET) -> pl.DataFrame:  # fmt: skip
    """The cumulative share of each class's genes at or below each fold change, at a preset."""
    classed = targets.classes(found, mirna, results, preset)
    parts = []
    for target in targets.CLASSES:
        x = np.sort(classed.filter(pl.col("target") == target)["log2FoldChange"].to_numpy())
        if x.size:
            parts.append(pl.DataFrame({"target": target, "log2FoldChange": x,
                                       "share": np.arange(1, x.size + 1) / x.size}))  # fmt: skip
    return pl.concat(parts)


def enrichment(found: targets.Targets, lfc: Mapping[str, float],
               results: pl.DataFrame) -> pl.DataFrame:  # fmt: skip
    """Each featured miRNA's targets among the genes called against its own change.

    A miRNA that rises should leave its targets among the down genes, one that falls among the
    up genes. Hypergeometric (``ctx.function.lite.ora``) against the tested genes in
    TargetScan's universe; validated and predicted-only targets are two sets; BH over all.

    Args:
        found: The target files.
        lfc: Each featured miRNA's log2 fold change from the miRNA fit.
        results: M5's tested genes: gene_id, gene_name, log2FoldChange and call.
    """
    background = results.filter(pl.col("gene_id").is_in(list(found.universe)))
    rows = []
    for mirna, change in lfc.items():
        if not np.isfinite(change):
            continue
        direction = "down" if change > 0 else "up"
        hits = background.filter(pl.col("call") == direction)["gene_id"].to_list()
        classed = targets.classes(found, mirna, background, targets.DEFAULT_PRESET)
        sets = {t: classed.filter(pl.col("target") == t)["gene_id"].to_list()
                for t in (targets.VALIDATED, targets.PREDICTED)}  # fmt: skip
        for result in ora(hits, background["gene_id"].to_list(), sets, min_size=1,
                          max_size=background.height):  # fmt: skip
            rows.append({"mirna": mirna, "target": result.name, "genes": direction,
                         "size": result.size, "hits": result.hits, "expected": result.expected,
                         "p": result.p})  # fmt: skip
    frame = pl.DataFrame(rows, schema={"mirna": pl.String, "target": pl.String,
                                       "genes": pl.String, "size": pl.Int64, "hits": pl.Int64,
                                       "expected": pl.Float64, "p": pl.Float64})  # fmt: skip
    q = adjust(frame["p"].to_numpy(), "padj") if frame.height else np.empty(0)
    return frame.with_columns(q=pl.Series(q, dtype=pl.Float64))


def placed(loci: pl.DataFrame, types: pl.DataFrame) -> pl.DataFrame:
    """Hairpins matched one to one with an annotation miRNA gene they overlap on its strand.

    Args:
        loci: Gather's hairpin loci: miRNA_ID, chrom, start, end, strand.
        types: The annotation's gene records, unversioned, with gene_type.
    """
    genes = types.filter(pl.col("gene_type") == "miRNA").select(
        "gene_id", "chrom", "strand", g_start="start", g_end="end"
    )
    hits = (loci.drop_nulls("chrom").join(genes, on=["chrom", "strand"])
            .filter((pl.col("start") <= pl.col("g_end"))
                    & (pl.col("g_start") <= pl.col("end"))))  # fmt: skip
    once = hits.filter((pl.col("miRNA_ID").count().over("miRNA_ID") == 1)
                       & (pl.col("gene_id").count().over("gene_id") == 1))  # fmt: skip
    return once.select("miRNA_ID", "gene_id")


def agreement(both: Layers, hairpins: pl.DataFrame, matched: pl.DataFrame,
              counts: pl.DataFrame) -> pl.DataFrame:  # fmt: skip
    """Per sample, Spearman rho of the hairpins' miRNA-library reads and their genes' STAR counts.

    Args:
        both: The samples with both layers.
        hairpins: ``hairpin`` then one column per miRNA aliquot.
        matched: ``placed``'s one-to-one hairpins and genes.
        counts: STAR counts, ``gene_id`` then one column per STAR aliquot.

    Returns:
        mirna_sample, rna_sample, condition and rho, one row per sample.
    """
    h_row = {h: i for i, h in enumerate(hairpins["hairpin"].to_list())}
    g_row = {g: i for i, g in enumerate(counts["gene_id"].to_list())}
    keep = matched.filter(pl.col("miRNA_ID").is_in(list(h_row)) & pl.col("gene_id").is_in(
        list(g_row)))  # fmt: skip
    small = hairpins.select(both.samples["mirna_sample"].to_list()).to_numpy()
    long = counts.select(both.samples["rna_sample"].to_list()).to_numpy()
    a = small[[h_row[h] for h in keep["miRNA_ID"]]].T.astype(float)
    b = long[[g_row[g] for g in keep["gene_id"]]].T.astype(float)
    return both.samples.select("mirna_sample", "rna_sample", "condition").with_columns(
        rho=pl.Series(neighbours.spearman(a, b)), hairpins=pl.lit(keep.height)
    )
