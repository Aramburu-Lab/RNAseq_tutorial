"""Synthetic count matrix with planted truth, for the tests and the Pyodide spike.

Not part of the ctx wheel: its consumers are tests/ and scripts/python/spike_pyodide.py,
which loads it by path. Nothing here comes from a real dataset.
"""

from dataclasses import dataclass

import numpy as np
import polars as pl

CELL_TYPES = ("type_a", "type_b")
CONDITIONS = ("control", "case")
# Fractions of n_genes, assigned in this order. Cell type is planted as the dominant axis
# (markers at |LFC| 4) so PC1 splits cell types before condition, as sorted immune cells do.
CLASSES = (
    ("marker", 0.10),
    ("condition", 0.05),
    ("interaction", 0.05),
    ("low", 0.10),
    ("null", 0.70),
)
LFC_MARKER = 4.0
LFC_CONDITION = 2.5
DISPERSION = 0.1
DONOR_SD_LOG2 = 0.25
SIZE_FACTOR_SD = 0.2
INT32_MAX = 2**31 - 1


@dataclass(frozen=True)
class Synthetic:
    """One synthetic experiment: counts, the sample sheet and the planted truth."""

    counts: pl.DataFrame
    samples: pl.DataFrame
    truth: pl.DataFrame


def sample_sheet(donors_per_condition: int) -> pl.DataFrame:
    """Donors nested in condition; every donor gives both cell types."""
    rows = []
    for d in range(2 * donors_per_condition):
        condition = CONDITIONS[d // donors_per_condition]
        # Batch alternates by donor, so it is balanced against condition and carries no effect.
        rows += [
            (f"d{d + 1:03d}_{cell}", f"d{d + 1:03d}", cell, condition, f"batch{d % 2 + 1}")
            for cell in CELL_TYPES
        ]
    return pl.DataFrame(
        rows, schema=["sample", "donor", "cell_type", "condition", "batch"], orient="row"
    )


def planted_truth(rng: np.random.Generator, n_genes: int) -> pl.DataFrame:
    """Gene classes, base means and planted log2 fold changes (type_b/type_a, case/control)."""
    sizes = [round(frac * n_genes) for _, frac in CLASSES[:-1]]
    sizes.append(n_genes - sum(sizes))
    gene_class = np.repeat([name for name, _ in CLASSES], sizes)
    sign = np.where(np.arange(n_genes) % 2 == 0, 1.0, -1.0)
    base = np.exp(rng.normal(np.log(200.0), 1.0, n_genes))
    low = gene_class == "low"
    base[low] = rng.uniform(0.5, 5.0, low.sum())
    lfc_cell = np.where(gene_class == "marker", sign * LFC_MARKER, 0.0)
    affected = np.isin(gene_class, ["condition", "interaction"])
    lfc_a = np.where(affected, sign * LFC_CONDITION, 0.0)
    lfc_b = np.where(gene_class == "condition", sign * LFC_CONDITION, 0.0)
    return pl.DataFrame(
        {
            "gene_id": [f"g{i + 1:05d}" for i in range(n_genes)],
            "gene_class": gene_class,
            "base_mean": base,
            "lfc_cell_type": lfc_cell,
            "lfc_condition_type_a": lfc_a,
            "lfc_condition_type_b": lfc_b,
        }
    )


def counts(seed: int = 42, n_genes: int = 1000, donors_per_condition: int = 6) -> Synthetic:
    """Negative-binomial counts with per-sample size factors and a per-donor effect."""
    rng = np.random.default_rng(seed)
    sheet = sample_sheet(donors_per_condition)
    truth = planted_truth(rng, n_genes)
    is_b = (sheet["cell_type"] == "type_b").to_numpy()
    is_case = (sheet["condition"] == "case").to_numpy()
    donor_idx = sheet["donor"].rank("dense").cast(pl.Int64).to_numpy() - 1

    donor_effect = rng.normal(0.0, DONOR_SD_LOG2, (n_genes, 2 * donors_per_condition))
    size_factor = np.exp(rng.normal(0.0, SIZE_FACTOR_SD, sheet.height))
    lfc_a = truth["lfc_condition_type_a"].to_numpy()[:, None]
    lfc_b = truth["lfc_condition_type_b"].to_numpy()[:, None]
    lfc_cond = np.where(is_b, lfc_b, lfc_a)
    log2_mu = (
        truth["lfc_cell_type"].to_numpy()[:, None] * is_b
        + lfc_cond * is_case
        + donor_effect[:, donor_idx]
    )
    mu = size_factor * truth["base_mean"].to_numpy()[:, None] * np.exp2(log2_mu)
    n = 1.0 / DISPERSION
    draws = rng.negative_binomial(n, n / (n + mu))
    assert draws.max() <= INT32_MAX, "counts exceed Int32; do not downcast"
    matrix = pl.DataFrame(draws.astype(np.int32), schema=sheet["sample"].to_list())
    return Synthetic(matrix.insert_column(0, truth["gene_id"]), sheet, truth)


def demo_box(synthetic: Synthetic, planted: int = 3, null: int = 5):  # returns ctx BoxData
    """A box-plot contract from the type_a samples: planted condition genes, then null genes."""
    from ctx.viz.box import box_data  # noqa: PLC0415 - keeps synth importable without ctx

    matrix = synthetic.counts.drop("gene_id").to_numpy().astype(float)
    log_cpm = np.log2(matrix / matrix.sum(axis=0) * 1e6 + 1)
    truth, sheet = synthetic.truth, synthetic.samples
    affected = truth.filter(pl.col("lfc_condition_type_a") != 0)["gene_id"].head(planted)
    nulls = truth.filter(pl.col("gene_class") == "null")["gene_id"].head(null)
    row = {g: i for i, g in enumerate(truth["gene_id"])}
    keep = [j for j in range(sheet.height) if sheet["cell_type"][j] == "type_a"]
    frame = pl.DataFrame(
        [
            (g, sheet["condition"][j], sheet["sample"][j], float(log_cpm[row[g], j]))
            for g in [*affected, *nulls]
            for j in keep
        ],
        schema=["feature", "group", "sample", "value"],
        orient="row",
    )
    columns = {
        "feature": ("none", "Synthetic gene identifier."),
        "group": ("none", "Condition of the donor."),
        "sample": ("none", "Sample identifier: donor and cell type."),
        "value": ("log2 CPM", "log2 of counts per million plus 1."),
    }
    return box_data(
        frame, columns=columns, feature_label="Gene", value_label="Expression",
        groups=("control", "case"),
    )  # fmt: skip
