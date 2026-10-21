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
# design_counts draws each gene's dispersion around DESeq2's own trend shape, a0 + a1 / mean,
# scattered log-normally. With one flat value there is no trend: pydeseq2's parametric fit
# did not converge and fell back to the mean on 313 of 480 fits of the test designs, and the
# prior's variance sat at its floor. These three, the best of 8 settings swept, gave 0 such
# fits in test_de and test_stats over seeds 0-29 and 100-119, and 1 in 240 fits of the
# CLI's contrasts (knowledge/decisions.md D-017).
DISPERSION_ASYMPTOTE = 0.05
DISPERSION_EXTRA_POISSON = 10.0
DISPERSION_SCATTER_SD = 0.3
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


def _log_cpm(synthetic: Synthetic) -> np.ndarray:
    """log2 CPM + 1, genes by samples."""
    matrix = synthetic.counts.drop("gene_id").to_numpy().astype(float)
    return np.log2(matrix / matrix.sum(axis=0) * 1e6 + 1)


def demo_volcano(synthetic: Synthetic, width: str = "double", labels: int = 6):  # ScatterData
    """A volcano contract from the type_a samples: Welch t-test per gene, BH, live calls."""
    from scipy.stats import false_discovery_control, ttest_ind  # noqa: PLC0415

    from ctx.de.stats import ALPHA, LFC_MIN, classify, p_line  # noqa: PLC0415
    from ctx.viz.scatter import CUTOFF_PRESETS, Calls, ScatterData  # noqa: PLC0415
    from ctx.viz.style import ROLES  # noqa: PLC0415

    sheet = synthetic.samples
    values = _log_cpm(synthetic)
    a = (sheet["cell_type"] == "type_a").to_numpy()
    case, control = (
        a & (sheet["condition"] == "case").to_numpy(),
        a & ~(sheet["condition"] == "case").to_numpy(),
    )
    p = ttest_ind(values[:, case], values[:, control], axis=1).pvalue
    lfc = values[:, case].mean(axis=1) - values[:, control].mean(axis=1)
    frame = pl.DataFrame({"gene_id": synthetic.truth["gene_id"], "lfc": lfc, "p": p}).with_columns(
        padj=pl.Series(false_discovery_control(np.nan_to_num(p, nan=1.0), method="bh")),
        neg_log10_p=-pl.col("p").log10(),
    )
    frame = frame.with_columns(
        call=classify(pl.col("padj"), pl.col("lfc")),
        keep=pl.col("p") < max(CUTOFF_PRESETS),
    )
    top = set(frame.filter(pl.col("call") != "ns").sort("p").head(labels)["gene_id"])
    frame = frame.with_columns(
        label=pl.when(pl.col("gene_id").is_in(top)).then(pl.col("gene_id")).otherwise(pl.lit(""))
    )
    columns = {
        "gene_id": ("none", "Synthetic gene identifier."),
        "lfc": ("log2 fold change", "Mean log2 CPM, case minus control, type_a samples."),
        "p": ("dimensionless", "Welch two-sample t-test p-value."),
        "padj": ("dimensionless", "Benjamini-Hochberg adjusted p over all genes."),
        "neg_log10_p": ("-log10 p", "Minus log10 of p."),
        "call": ("none", "up, down or ns at padj < 0.05 and |lfc| >= 1."),
        "keep": ("none", "p below the loosest cutoff offered; always drawn on the page."),
        "label": ("none", "Gene printed beside its point; empty for none."),
    }
    passing = frame.filter(pl.col("padj") < ALPHA)["p"]
    calls = Calls(
        lfc="lfc", p="p", adjusted={"padj": "BH-adjusted p"}, raw_words="raw p",
        tested="case", reference="control", alpha=ALPHA, lfc_min=LFC_MIN,
        max_cutoff=max(CUTOFF_PRESETS), n_tested=frame.height, p_floor=1e-300,
        p_line=p_line("padj", ALPHA, frame.height, passing.len()),
        reading=("Each dot is one of {n_tested} genes.", "Up in {tested}: {line_rule}."),
        summary="Up in {tested}: {n_up}. Down: {n_down}.",
        line_rule="above -log10 p = {y_line}", no_line_rule="none passes {basis} < {cutoff}",
        unit="genes",
    )  # fmt: skip
    return ScatterData(
        frame=frame, columns=columns, x="lfc", y="neg_log10_p", x_label="Fold change",
        y_label="Significance", colour="call",
        colours={"up": ROLES["comparison"], "down": ROLES["reference"], "ns": None},
        caption=f"Welch t-test, BH over {frame.height} genes; n = 6 control, 6 case.",
        label="label", hover=("gene_id", "padj"), width=width, keep="keep", calls=calls,
    )  # fmt: skip


def demo_ma(synthetic: Synthetic):  # ScatterData
    """An MA contract: mean counts on a log x axis, the way M5's MA plot draws it."""
    import dataclasses  # noqa: PLC0415

    volcano = demo_volcano(synthetic)
    mean = synthetic.counts.select(pl.exclude("gene_id")).mean_horizontal() + 1
    frame = volcano.frame.with_columns(mean_count=mean)
    columns = {**volcano.columns, "mean_count": ("reads", "Mean count over samples, plus 1.")}
    return dataclasses.replace(
        volcano, frame=frame, columns=columns, x="mean_count", y="lfc", x_label="Mean count",
        y_label="Fold change", log_x=True, calls=None, keep=None, label=None,
        caption="Mean count against log2 fold change; no test.",
    )  # fmt: skip


def demo_pca(synthetic: Synthetic):  # ScatterData
    """A PCA contract: all samples, coloured by cell type, shaped by condition."""
    from ctx.viz.annotate import permanova  # noqa: PLC0415
    from ctx.viz.scatter import ScatterData  # noqa: PLC0415
    from ctx.viz.style import category_colours  # noqa: PLC0415

    values = _log_cpm(synthetic)
    centred = values - values.mean(axis=1, keepdims=True)
    u, s, _ = np.linalg.svd(centred.T, full_matrices=False)
    scores = u[:, :2] * s[:2]
    sheet = synthetic.samples
    test = permanova(scores, sheet["cell_type"].to_list(), strata=sheet["donor"].to_list())
    frame = sheet.with_columns(pc1=pl.Series(scores[:, 0]), pc2=pl.Series(scores[:, 1]))
    columns = {
        "sample": ("none", "Sample identifier."),
        "donor": ("none", "Donor identifier."),
        "cell_type": ("none", "Cell type."),
        "condition": ("none", "Condition."),
        "batch": ("none", "Batch."),
        "pc1": ("log2 CPM", "Score on principal component 1."),
        "pc2": ("log2 CPM", "Score on principal component 2."),
    }
    return ScatterData(
        frame=frame, columns=columns, x="pc1", y="pc2", x_label="PC1", y_label="PC2",
        colour="cell_type", colours=category_colours(list(CELL_TYPES)),
        caption=test.text("cell type"), symbol="condition", hover=("sample", "donor"),
    )  # fmt: skip


def demo_bar(synthetic: Synthetic, horizontal: bool = False):  # BarData
    """A bar contract: samples per cell type and condition, stacked by batch."""
    from ctx.viz.bar import BarData  # noqa: PLC0415
    from ctx.viz.style import category_colours  # noqa: PLC0415

    frame = (
        synthetic.samples.with_columns(group=pl.format("{} {}", "cell_type", "condition"))
        .group_by("group", "batch", maintain_order=True)
        .agg(n=pl.len().cast(pl.Int64))
    )
    columns = {
        "group": ("none", "Cell type and condition."),
        "batch": ("none", "Batch."),
        "n": ("samples", "Samples in this group and batch."),
    }
    return BarData(
        frame=frame, columns=columns, category="group", value="n", category_label="Group",
        value_label="Samples", series="batch", colours=category_colours(["batch1", "batch2"]),
        caption=f"{synthetic.samples.height} samples; counts, no test.", horizontal=horizontal,
        stacked=True, hline=(3.0, "half a group"), hover=("batch",),
    )  # fmt: skip


def demo_bars_at_scale(horizontal: bool):  # BarData
    """Bars at Step-1 scale: 20 labelled groups stacked by 2 batches, or 210 unlabelled runs."""
    from ctx.viz.bar import BarData  # noqa: PLC0415
    from ctx.viz.style import category_colours  # noqa: PLC0415

    rng = np.random.default_rng(3)
    if horizontal:
        groups = [f"cell type {i // 5} condition {i % 5}" for i in range(20)]
        frame = pl.DataFrame(
            {"group": groups * 2, "batch": ["short"] * 20 + ["long"] * 20,
             "n": rng.integers(0, 15, 40).astype(np.int64)}
        )  # fmt: skip
        category, value, label = "group", "n", "Samples"
        columns = {"group": ("none", "Group."), "batch": ("none", "Batch."),
                   "n": ("samples", "Samples.")}  # fmt: skip
    else:
        frame = pl.DataFrame(
            {"run": [f"r{i:03d}" for i in range(210)], "batch": ["short"] * 111 + ["long"] * 99,
             "million_reads": rng.uniform(5.5, 232.0, 210)}
        )  # fmt: skip
        category, value, label = "run", "million_reads", "Reads"
        columns = {"run": ("none", "Run."), "batch": ("none", "Batch."),
                   "million_reads": ("million reads", "Reads / 1e6.")}  # fmt: skip
    return BarData(
        frame=frame, columns=columns, category=category, value=value, category_label="Group",
        value_label=label, series="batch", colours=category_colours(["short", "long"]),
        caption="Counts at the Step-1 scale; no test.", horizontal=horizontal,
        stacked=True, width="double",
    )  # fmt: skip


def design_counts(  # noqa: PLR0913 - each argument is one knob of the planted design
    groups: tuple[tuple, ...] = (("control", "short", 6), ("case", "short", 6)),
    *,
    seed: int = 7,
    n_genes: int = 1000,
    paired: bool = False,
    donor_sd: float = DONOR_SD_LOG2,
    batch_lfc: float = 0.0,
    condition_lfc: float = LFC_CONDITION,
) -> Synthetic:
    """Counts for DE tests: declared (condition, batch, n) groups; planted condition, batch genes.

    Separate from ``counts()``, whose arguments the P1.5 spike record pins. Genes 0-9 % carry
    a condition effect
    (case over control, alternating sign, |LFC| ``condition_lfc``); genes 10-19 % carry
    a batch effect (long over short, ``batch_lfc``); the rest are null. With ``paired`` the
    i-th sample of every group comes from donor i.

    Returns:
        Counts (``gene_id`` then one column per sample), a curated-style sheet, and the truth.
    """
    rng = np.random.default_rng(seed)
    rows, donor_of = [], []
    for g, (condition, batch, n, *extra) in enumerate(groups):
        cell = extra[0] if extra else "type_a"
        # A fifth element names the individuals: groups sharing a prefix share donors.
        prefix = extra[1] if len(extra) > 1 else None
        for i in range(n):
            donor = f"{prefix}{i:03d}" if prefix else (f"d{i:03d}" if paired else f"d{g}{i:02d}")
            rows.append((f"XRR{len(rows) + 1:06d}", f"S{len(rows) + 1:03d}", donor, cell,
                         condition, batch, "", "M" if i % 2 else "F", "PB"))  # fmt: skip
            donor_of.append(donor)
    sheet = pl.DataFrame(
        rows, orient="row",
        schema=["run_accession", "sample", "donor", "cell_type", "condition", "batch",
                "time_point", "sex", "tissue"],
    )  # fmt: skip
    sign = np.where(np.arange(n_genes) % 2 == 0, 1.0, -1.0)
    index = np.arange(n_genes)
    lfc_condition = np.where(index < 0.1 * n_genes, condition_lfc * sign, 0.0)
    lfc_batch = np.where((index >= 0.1 * n_genes) & (index < 0.2 * n_genes), batch_lfc * sign, 0.0)
    base = np.exp(rng.normal(np.log(200.0), 1.0, n_genes))
    donors = sorted(set(donor_of))
    effect = rng.normal(0.0, donor_sd, (n_genes, len(donors)))[
        :, [donors.index(d) for d in donor_of]
    ]
    is_case = (sheet["condition"] == "case").to_numpy()
    is_long = (sheet["batch"] == "long").to_numpy()
    is_b = (sheet["cell_type"] == "type_b").to_numpy()
    # Genes 20-30 % mark type_b (|LFC| 4), so cell type leads the PCA as in sorted cells.
    lfc_cell = np.where((index >= 0.2 * n_genes) & (index < 0.3 * n_genes), 4.0 * sign, 0.0)
    log2_mu = (
        lfc_condition[:, None] * is_case
        + lfc_batch[:, None] * is_long
        + lfc_cell[:, None] * is_b
        + effect
    )
    mu = np.exp(rng.normal(0.0, SIZE_FACTOR_SD, sheet.height)) * base[:, None] * np.exp2(log2_mu)
    alpha = (DISPERSION_ASYMPTOTE + DISPERSION_EXTRA_POISSON / base) * np.exp(
        rng.normal(0.0, DISPERSION_SCATTER_SD, n_genes)
    )
    n = 1.0 / alpha[:, None]
    draws = rng.negative_binomial(n, n / (n + mu)).astype(np.int32)
    truth = pl.DataFrame({
        "gene_id": [f"g{i + 1:05d}" for i in range(n_genes)],
        "lfc_condition": lfc_condition,
        "lfc_batch": lfc_batch,
    })  # fmt: skip
    matrix = pl.DataFrame(draws, schema=sheet["sample"].to_list(), orient="row")
    return Synthetic(matrix.insert_column(0, truth["gene_id"]), sheet, truth)


# Null genes at mean 20, dispersed as design_counts disperses them. 1,500 of them make
# independent filtering pay on every seed measured (at least 141 genes removed over seeds
# 0-29 and 100-119); 600 left it removing nothing on 7 of 30 seeds, and a flat dispersion on
# them brought back the trend warning.
LOW_NULLS = 1500
LOW_NULL_MEAN = 20.0


def filtered_counts(synthetic: Synthetic, outlier: str = "g00500") -> pl.DataFrame:
    """The counts plus ``LOW_NULLS`` low null genes and one Cook's outlier in ``outlier``.

    Both of pydeseq2's filters then remove genes, which a test needs whenever the tested set
    must differ from every gene: the reasons, BH over the tested set, the volcano's m.
    """
    rng = np.random.default_rng(1)
    names = synthetic.samples["sample"].to_list()
    alpha = (DISPERSION_ASYMPTOTE + DISPERSION_EXTRA_POISSON / LOW_NULL_MEAN) * np.exp(
        rng.normal(0.0, DISPERSION_SCATTER_SD, LOW_NULLS)
    )
    n = 1.0 / alpha
    low = pl.DataFrame(
        {"gene_id": [f"low{i:04d}" for i in range(LOW_NULLS)],
         **{c: rng.negative_binomial(n, n / (n + LOW_NULL_MEAN)).astype(np.int32)
            for c in names}}
    )  # fmt: skip
    first = names[0]
    return pl.concat([synthetic.counts, low]).with_columns(
        pl.when(pl.col("gene_id") == outlier).then(pl.col(first) * 100 + 5000)
        .otherwise(pl.col(first)).alias(first)
    )  # fmt: skip


SEX_GENES = ("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")


def write_gather_outputs(root, dataset: str, synthetic: Synthetic, day, ensembl=False) -> None:
    """Sheet, counts and genes as gather writes them, from a design_counts experiment.

    The last six genes become the sex genes, expressed by the sheet's sex: XIST in F, the
    chrY genes in M. With ``ensembl`` each gene's ID is a made-up Ensembl ID, ENSG00090000001
    on, under the same symbol, so M5's picker can link it (``ctx.de.page.picking``).
    """
    from ctx.gather.harmonise import MODULE  # noqa: PLC0415 - keeps synth importable without ctx
    from ctx.io.naming import stem  # noqa: PLC0415
    from ctx.io.readme_companion import write_table  # noqa: PLC0415

    out = root / "outputs" / "01_gather"
    counts = synthetic.counts
    ids = counts["gene_id"].to_list()
    names = [g.upper() for g in ids[: -len(SEX_GENES)]] + list(SEX_GENES)
    female = set(synthetic.samples.filter(pl.col("sex") == "F")["sample"])
    sex_ids = ids[-len(SEX_GENES) :]
    counts = counts.with_columns([
        pl.when(pl.col("gene_id") == sex_ids[0]).then(3000 if c in female else 0)
        .when(pl.col("gene_id").is_in(sex_ids[1:])).then(0 if c in female else 800)
        .otherwise(pl.col(c)).cast(pl.Int32).alias(c)
        for c in counts.columns[1:]
    ])  # fmt: skip
    if ensembl:
        made_up = {g: f"ENSG{90_000_000 + i + 1:011d}" for i, g in enumerate(ids)}
        counts = counts.with_columns(pl.col("gene_id").replace_strict(made_up))
        ids = [made_up[g] for g in ids]
    genes = pl.DataFrame({"gene_id": ids, "gene_name": names})
    sheet = synthetic.samples.with_columns(donor_sex=pl.lit(""), replicate=pl.lit(1, pl.Int32))
    described = {c: ("none", c) for c in {*sheet.columns, *counts.columns, *genes.columns}}
    for what, frame in (("samplesheet", sheet), ("counts", counts), ("genes", genes)):
        table_stem = stem(MODULE, f"{dataset}_{what}", 1, day=day)
        write_table(frame, out, table_stem, columns={c: described[c] for c in frame.columns},
                    data_releases={})  # fmt: skip


# The groups tests/fixtures/datasets/de_synthetic.toml declares its contrasts on.
DE_GROUPS = (
    ("control", "short", 6, "type_a"),
    ("case", "short", 6, "type_a"),
    ("control", "long", 3, "type_a"),
    ("case", "long", 6, "type_a"),
    ("control", "short", 6, "type_b"),
    ("case", "short", 6, "type_b"),
)


# Mean reads under which a type_b contrast's expression filter (10 reads at the median library)
# drops a gene, with room below it, and over which every type_a contrast keeps one.
TYPE_B_SILENT, TYPE_A_EXPRESSED = 8.0, 20.0


def type_a_only(synthetic: Synthetic) -> list[str]:
    """Genes expressed in type_a and silent in type_b: a set of them is tested by the type_a
    contrasts alone, so a contrast's set list differs from another's (TODO.md V2.1)."""
    cell = synthetic.samples["cell_type"]
    a, b = (synthetic.samples.filter(cell == c)["sample"].to_list() for c in CELL_TYPES)
    mean_a, mean_b = (synthetic.counts.select(s).mean_horizontal().to_numpy() for s in (a, b))
    keep = (mean_b < TYPE_B_SILENT) & (mean_a >= TYPE_A_EXPRESSED)
    return [g.upper() for g, k in zip(synthetic.counts["gene_id"], keep, strict=True) if k]


def write_gene_sets(
    path, n_genes: int = 1000, seed: int = 5, n_random: int = 20, synthetic: Synthetic | None = None
) -> None:
    """A GMT for design_counts genes: UP_SET and DOWN_SET hold the planted condition genes
    (uppercase symbols, as write_gather_outputs names them); twenty random sets hold others.
    With ``synthetic``, TYPE_A_SET holds its ``type_a_only`` genes."""
    rng = np.random.default_rng(seed)
    planted = [f"G{i + 1:05d}" for i in range(int(0.1 * n_genes))]
    others = [f"G{i + 1:05d}" for i in range(int(0.3 * n_genes), n_genes - len(SEX_GENES))]
    # At most 100 genes each, near the Step-1 IFN-alpha set's 97 and under GSEA's 500.
    sets = {"UP_SET": planted[0::2][:100], "DOWN_SET": planted[1::2][:100]}
    sets |= {
        f"RANDOM_{j:02d}": list(rng.choice(others, int(rng.integers(15, 200)), replace=False))
        for j in range(n_random)
    }
    if synthetic is not None:
        sets["TYPE_A_SET"] = type_a_only(synthetic)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f"{name}\tsynthetic\t" + "\t".join(genes) + "\n" for name, genes in sets.items())
    )


def demo_strip(log: bool = False, width: str = "double"):  # returns ctx StripData
    """A QC strip: one metric in three strata of eight, one planted low sample, two bands.

    The grey bars are each stratum's median ± 3 scaled MAD, computed here with scipy, so the
    figure's numbers are checked against a second computation.
    """
    from scipy.stats import median_abs_deviation  # noqa: PLC0415

    from ctx.viz.strip import Range, Span, StripData  # noqa: PLC0415
    from ctx.viz.style import category_colours  # noqa: PLC0415

    rng = np.random.default_rng(3)
    centres = {"PB / NDN": 9000.0, "PB / Mono": 12000.0, "BM / NDN": 10000.0}
    scale = 1000.0 if log else 1.0
    rows = [
        (f"S{8 * k + i:03d}", stratum, "case" if i % 2 else "control",
         float(centre * scale * np.exp(rng.normal(0.0, 0.05))))
        for k, (stratum, centre) in enumerate(centres.items()) for i in range(8)
    ]  # fmt: skip
    rows[3] = (*rows[3][:3], 4000.0 * scale)
    frame = pl.DataFrame(rows, schema=["sample", "stratum", "condition", "value"], orient="row")
    ranges = []
    for stratum in centres:
        values = frame.filter(pl.col("stratum") == stratum)["value"].to_numpy()
        x = np.log10(values) if log else values
        median, mad = np.median(x), median_abs_deviation(x, scale="normal")
        ends = np.array([median, median - 3 * mad, median + 3 * mad])
        centre, low, high = (float(v) for v in (10**ends if log else ends))
        ranges.append(Range(stratum, centre, low, high))
    low = {r.level: r.low for r in ranges}
    frame = frame.with_columns(
        flagged=pl.col("value") < pl.col("stratum").replace_strict(low, return_dtype=pl.Float64)
    ).with_columns(
        reason=pl.when("flagged").then(pl.lit("below its stratum")).otherwise(pl.lit(""))
    )
    columns = {
        "sample": ("none", "Synthetic sample."), "stratum": ("none", "Tissue / cell type."),
        "condition": ("none", "Condition."),
        "value": ("reads" if log else "genes", "Library size." if log else "Genes detected."),
        "flagged": ("none", "Below its stratum's median minus 3 scaled MAD."),
        "reason": ("none", "Why the sample is flagged."),
    }  # fmt: skip
    return StripData(
        frame=frame, columns=columns, category="stratum", value="value",
        category_label="Stratum", value_label="Library size" if log else "Genes detected",
        colour="condition", colours=category_colours(["control", "case"]),
        caption="Grey bars: median ± 3 scaled MAD per stratum; triangles are flagged. "
        "Shaded: the metric's absolute bands. Synthetic data; no test.",
        flagged="flagged", flag_label="flagged", hover=("sample", "reason"),
        spans=(Span("bad", "bad", None, 5000.0 * scale),
               Span("good", "good", 8000.0 * scale, None)),
        ranges=tuple(ranges), range_label="median ± 3 MAD", log=log, width=width,
    )  # fmt: skip


def demo_bubble(width: str = "double"):  # returns ctx ScatterData
    """Read fates as bubbles: raw against useful fragments, size the MT + rRNA share.

    Log axes, the identity line, and the median ± 3 scaled MAD of the useful fraction as
    three lines y = f x; one planted sample keeps a third of what the others keep.
    """
    from scipy.stats import median_abs_deviation  # noqa: PLC0415

    from ctx.viz.scatter import Line, ScatterData  # noqa: PLC0415
    from ctx.viz.style import category_colours  # noqa: PLC0415

    rng = np.random.default_rng(11)
    cells = ("NDN", "LDN", "Monocyte")
    raw = np.exp(rng.normal(np.log(2.5e7), 0.6, 24))
    useful_share = np.clip(rng.normal(0.8, 0.03, 24), 0, 1)
    useful_share[5] = 0.27
    frame = pl.DataFrame({
        "sample": [f"S{i:03d}" for i in range(24)], "cell_type": [cells[i % 3] for i in range(24)],
        "raw": raw, "useful": raw * useful_share, "useful_share": useful_share,
        "mt_rrna_share": np.clip(rng.normal(0.08, 0.03, 24), 0.0, 1.0),
    })  # fmt: skip
    median = float(np.median(useful_share))
    mad = float(median_abs_deviation(useful_share, scale="normal"))
    ends = (float(raw.min()) * 0.8, float(raw.max()) * 1.2)
    lines = (Line("y = x", ends, ends),) + tuple(
        Line(name, ends, (f * ends[0], f * ends[1]), dash=True)
        for name, f in (("median ± 3 MAD", median - 3 * mad), ("median ± 3 MAD", median + 3 * mad))
    )  # fmt: skip
    columns = {
        "sample": ("none", "Synthetic sample."), "cell_type": ("none", "Cell type."),
        "raw": ("fragments", "Fragments sequenced."),
        "useful": ("fragments", "Fragments Salmon assigned to a transcript."),
        "useful_share": ("fraction", "useful / raw."),
        "mt_rrna_share": ("fraction", "Share of counted reads on MT and rRNA genes."),
    }  # fmt: skip
    return ScatterData(
        frame=frame, columns=columns, x="raw", y="useful", x_label="Raw", y_label="Useful",
        colour="cell_type", colours=category_colours(list(cells)),
        caption="Bubble area: MT + rRNA share. Dashed: the useful fraction's median ± 3 scaled "
        "MAD. Synthetic data; no test.",
        hover=("sample", "useful_share", "mt_rrna_share"), lines=lines, log_x=True, log_y=True,
        width=width, size="mt_rrna_share", size_label="MT + rRNA",
    )  # fmt: skip


def demo_curves(width: str = "double"):  # returns ctx CurveData
    """Mean quality per position for eight samples, one with a degraded tail.

    FastQC's three shaded regions and its median thresholds (limits.txt, 0.12.1) over them.
    """
    from ctx.viz.curve import CurveData  # noqa: PLC0415
    from ctx.viz.scatter import Line  # noqa: PLC0415
    from ctx.viz.strip import Span  # noqa: PLC0415
    from ctx.viz.style import category_colours  # noqa: PLC0415

    rng = np.random.default_rng(5)
    positions = np.arange(1, 102)
    rows = []
    for i in range(8):
        quality = 36.0 - 0.02 * positions + rng.normal(0.0, 0.3, positions.size)
        if i == 3:
            quality[-30:] = np.linspace(30.0, 2.0, 30)
        rows += [(f"S{i:03d}", "degraded" if i == 3 else "typical", int(x), float(q))
                 for x, q in zip(positions, quality, strict=True)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["sample", "group", "position", "quality"], orient="row")
    columns = {
        "sample": ("none", "Synthetic sample."), "group": ("none", "Planted degraded or typical."),
        "position": ("nt", "Position in the read."), "quality": ("Phred", "Mean base quality."),
    }  # fmt: skip
    ends = (1.0, 101.0)
    return CurveData(
        frame=frame, columns=columns, series="sample", x="position", y="quality",
        x_label="Position in read", y_label="Mean quality", colour="group",
        colours=category_colours(["typical", "degraded"]),
        caption="Shaded as FastQC shades it: good from Phred 28, reasonable from 20. Dashed: "
        "FastQC's median warn (25) and fail (20). Synthetic data; no test.",
        spans=(Span("good", "good", 28.0, None), Span("usable", "reasonable", 20.0, 28.0),
               Span("bad", "poor", None, 20.0)),
        thresholds=(Line("warn, fail", ends, (25.0, 25.0), dash=True),
                    Line("warn, fail", ends, (20.0, 20.0), dash=True)),
        width=width,
    )  # fmt: skip


def demo_heatmap(n_genes: int = 30):  # returns ctx HeatmapData
    """z-scored log CPM of the planted genes, clustered both ways, two tracks, a p per row."""
    from scipy.stats import false_discovery_control, ttest_ind  # noqa: PLC0415

    from ctx.viz.heatmap import HeatmapData, Track, cluster_order, diverging_levels  # noqa: PLC0415
    from ctx.viz.style import DATA_COLOURS  # noqa: PLC0415

    synthetic = design_counts(
        (("control", "short", 6, "type_a"), ("case", "short", 6, "type_a"),
         ("control", "short", 6, "type_b"), ("case", "short", 6, "type_b")), seed=3)  # fmt: skip
    sheet = synthetic.samples
    log_cpm = _log_cpm(synthetic)[:n_genes]
    genes = synthetic.counts["gene_id"].to_list()[:n_genes]
    z = (log_cpm - log_cpm.mean(axis=1, keepdims=True)) / log_cpm.std(axis=1, keepdims=True)
    samples = sheet["sample"].to_list()
    case = (sheet["condition"] == "case").to_numpy()
    p = ttest_ind(log_cpm[:, case], log_cpm[:, ~case], axis=1).pvalue
    stats = dict(zip(genes, (float(q) for q in false_discovery_control(p)), strict=True))
    values = [float(v) for row in z for v in row]
    levels, swatches = diverging_levels(values)
    frame = pl.DataFrame({
        "gene": [g for g in genes for _ in samples], "sample": samples * len(genes),
        "z": values, "level": levels,
    })  # fmt: skip
    tracks = (
        Track("cell type", dict(zip(samples, sheet["cell_type"], strict=True)),
              {"type_a": DATA_COLOURS[2], "type_b": DATA_COLOURS[6]}),
        Track("condition", dict(zip(samples, sheet["condition"], strict=True)),
              {"control": DATA_COLOURS[1], "case": DATA_COLOURS[0]}),
    )  # fmt: skip
    return HeatmapData(
        frame=frame, row="gene", column="sample", value="z", level="level",
        columns={"gene": ("none", "Synthetic gene."), "sample": ("none", "Synthetic sample."),
                 "z": ("SD", "log2 CPM, centred and scaled per gene."),
                 "level": ("none", "The z step the cell is drawn in.")},
        row_label="Gene", column_label="Sample", swatches=swatches,
        caption="z-scored log2 CPM; rows and columns ordered by average-linkage clustering on "
        "1 − Pearson. Right: BH-adjusted Welch t-test p, case against control. Synthetic data.",
        row_order=cluster_order(z, genes), column_order=cluster_order(z.T, samples),
        tracks=tracks, row_stats=stats, row_stats_label="BH-adjusted p",
        row_stats_format=".1g",
    )  # fmt: skip


def demo_band_map():  # returns ctx HeatmapData
    """Samples × metrics coloured by QC band, one metric with no band."""
    from ctx.viz.heatmap import HeatmapData, Swatch  # noqa: PLC0415
    from ctx.viz.style import BAND_COLOURS  # noqa: PLC0415

    metrics = ("Reads", "Q30 bases", "Adapter content", "Genes detected")
    bands = ("good", "usable", "bad", "no band")
    rows = [(f"S{i:03d}", m, bands[(i * 7 + j * 3) % 4] if m != "Genes detected" else "no band")
            for i in range(12) for j, m in enumerate(metrics)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["sample", "metric", "band"], orient="row", strict=False)
    frame = frame.with_columns(value=pl.lit(None, pl.Float64))
    return HeatmapData(
        frame=frame, row="sample", column="metric", value="value", level="band",
        columns={"sample": ("none", "Synthetic sample."), "metric": ("none", "QC metric."),
                 "band": ("none", "The metric's band; empty when it has none."),
                 "value": ("none", "Not used by this map.")},
        row_label="Sample", column_label="Metric",
        swatches=(*(Swatch(b, BAND_COLOURS[b]) for b in ("good", "usable", "bad")),
                  Swatch("no band", None)),
        caption="Each metric's absolute band per sample; grey: no band. Synthetic data; no test.",
    )  # fmt: skip


def demo_flow(width: str = "double"):  # returns ctx FlowData
    """Where the reads of a cohort go: sequenced, Trim Galore's fates, Salmon's fates."""
    from ctx.viz.flow import FlowData, Node  # noqa: PLC0415
    from ctx.viz.style import DATA_COLOURS  # noqa: PLC0415

    nodes = (
        Node("raw", "sequenced", 0, DATA_COLOURS[4]),
        Node("kept", "written", 1, DATA_COLOURS[2]), Node("short", "too short", 1),
        Node("n", "too many N", 1),
        Node("assigned", "assigned", 2, DATA_COLOURS[2]), Node("unassigned", "not assigned", 2),
        Node("decoy", "decoy", 2, DATA_COLOURS[6]),
    )  # fmt: skip
    links = pl.DataFrame(
        [("raw", "kept", 961_204.0), ("raw", "short", 36_512.0), ("raw", "n", 2_284.0),
         ("kept", "assigned", 781_530.0), ("kept", "unassigned", 151_020.0),
         ("kept", "decoy", 28_654.0)],
        schema=["source", "target", "value"], orient="row",
    )  # fmt: skip
    columns = {"source": ("none", "Node the reads leave."), "target": ("none", "Node they reach."),
               "value": ("reads", "Reads moving along the link.")}  # fmt: skip
    return FlowData(
        nodes=nodes, links=links, columns=columns,
        stage_labels=("Sequenced", "Trim Galore", "Salmon"),
        caption="Every read of the cohort, summed over samples: Trim Galore's fates, then "
        "Salmon's for the reads it wrote. Synthetic data; no test.",
        width=width,
    )  # fmt: skip


def demo_surface(width: str = "double"):  # returns ctx SurfaceData
    """Reads kept across fastp -q and -l, a smooth synthetic surface; fastp's default marked."""
    from ctx.viz.contour import Mark, SurfaceData  # noqa: PLC0415

    rows = [
        (q, n, float(100.0 / (1.0 + np.exp((q - 28.0) / 2.5)) * (1.0 - (n / 150.0) ** 2)))
        for q in range(5, 40, 2) for n in range(15, 106, 5)
    ]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["q", "l", "kept"], orient="row")
    return SurfaceData(
        frame=frame, x="q", y="l", z="kept", x_label="-q", y_label="-l", z_label="Reads kept",
        columns={"q": ("Phred", "fastp -q: a base below this is unqualified."),
                 "l": ("nt", "fastp -l: a read shorter than this is dropped."),
                 "kept": ("%", "Reads that pass, of the reads in.")},
        caption="Reads fastp keeps at each -q and -l, the other settings at fastp's defaults; "
        "eight equal bands; the line marks 80 %. Synthetic surface, 18 x 19 grid.",
        isolines=(80.0,), isoline_format=".0f", mark=Mark(15.0, 15.0, "fastp default"),
        width=width,
    )  # fmt: skip


def demo_density(width: str = "double"):  # returns ctx DensityData
    """Library size against genes detected: two groups with densities, one of 12 as points."""
    import dataclasses  # noqa: PLC0415

    from ctx.viz.density import DensityData  # noqa: PLC0415
    from ctx.viz.style import category_colours  # noqa: PLC0415

    rng = np.random.default_rng(21)
    rows = []
    for group, n, depth, genes in (("NDN", 60, 2.0e7, 9000.0), ("Monocyte", 45, 3.0e7, 12500.0),
                                   ("LDN", 12, 2.5e7, 10000.0)):  # fmt: skip
        size = np.exp(rng.normal(np.log(depth), 0.35, n))
        detected = genes + 1800.0 * np.log(size / depth) + rng.normal(0.0, 500.0, n)
        rows += [(f"{group}_{i:02d}", group, float(s), float(g)) for i, (s, g) in
                 enumerate(zip(size, detected, strict=True))]  # fmt: skip
    schema = ["sample", "cell_type", "library_size", "genes"]
    frame = pl.DataFrame(rows, schema=schema, orient="row")
    frame = frame.with_columns(flagged=pl.col("sample").is_in(["NDN_03", "NDN_17", "Monocyte_05"]))
    data = DensityData(
        frame=frame, x="library_size", y="genes", group="cell_type",
        x_label="Library size", y_label="Genes detected",
        columns={"sample": ("none", "Synthetic sample."), "cell_type": ("none", "Cell type."),
                 "library_size": ("reads", "Reads counted on genes."),
                 "genes": ("genes", "Genes with at least one read."),
                 "flagged": ("none", "Flagged by the QC checkpoint.")},
        colours=category_colours(["NDN", "Monocyte", "LDN"]), caption="", highlight="flagged",
        highlight_label="flagged", hover=("sample", "flagged"), log_x=True, width=width,
    )  # fmt: skip
    n = {g: d.n for g, d in data.densities.items()}
    caption = (f"Solid: 50 %, dashed: 90 % highest-density regions of NDN (n = {n['NDN']}) and "
               f"Monocyte (n = {n['Monocyte']}); LDN (n = 12) is drawn as points, too few for a "
               f"density. {data.bandwidth()}. Synthetic data; no test.")  # fmt: skip
    return dataclasses.replace(data, caption=caption)


# Okabe and Ito's own names for the seven data colours, in DATA_COLOURS order.
COLOUR_NAMES = ("orange", "sky blue", "bluish green", "yellow", "blue", "vermillion",
                "reddish purple")  # fmt: skip


def demo_linked(n: int = 48):  # returns ctx LinkedData
    """Samples on two PCs beside their tissue, cell type, condition, batch and QC status."""
    from ctx.viz.linked import LinkedData  # noqa: PLC0415
    from ctx.viz.style import BAND_COLOURS  # noqa: PLC0415

    rng = np.random.default_rng(4)
    cells = rng.choice(["NDN", "Monocyte", "LDN"], n, p=[0.5, 0.35, 0.15])
    shift = {"NDN": (-6.0, 0.0), "Monocyte": (6.0, 1.0), "LDN": (-2.0, 5.0)}
    frame = pl.DataFrame({
        "sample": [f"S{i:03d}" for i in range(n)],
        "pc1": [shift[c][0] + rng.normal() for c in cells],
        "pc2": [shift[c][1] + rng.normal() for c in cells],
        "tissue": rng.choice(["PB", "BM"], n, p=[0.8, 0.2]).tolist(),
        "cell_type": cells.tolist(),
        "condition": rng.choice(["steady_state", "gcsf", "pdac"], n).tolist(),
        "batch": rng.choice(["short", "long"], n).tolist(),
        "qc_status": rng.choice(["pass", "flag"], n, p=[0.85, 0.15]).tolist(),
    })  # fmt: skip
    return LinkedData(
        frame=frame, x="pc1", y="pc2", x_label="PC1", y_label="PC2",
        columns={"sample": ("none", "Synthetic sample."), "pc1": ("VST units", "PC1 score."),
                 "pc2": ("VST units", "PC2 score."), "tissue": ("none", "Tissue."),
                 "cell_type": ("none", "Cell type."), "condition": ("none", "Condition."),
                 "batch": ("none", "Read-length batch."),
                 "qc_status": ("none", "QC status: pass or flag.")},
        dimensions=("tissue", "cell_type", "condition", "batch", "qc_status"),
        dimension_labels=("Tissue", "Cell type", "Condition", "Batch", "QC"),
        colour_dimension="qc_status",
        colours={"pass": BAND_COLOURS["good"], "flag": BAND_COLOURS["bad"]},
        caption="Samples by tissue, cell type, condition, batch and QC status; bands into QC "
        "are coloured by status. Synthetic data; no test.",
        names=COLOUR_NAMES, summary="{n} painted {colour}: {detail}",
        empty="Nothing painted yet.", detail="cell_type", hover=("sample", "qc_status"),
    )  # fmt: skip
