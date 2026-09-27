"""M5: ``ctx de``. Fit every declared contrast, then write the tables, the figures and the page.

Reads the newest gather outputs of one dataset (sheet, counts, genes) and the dataset file's
``[[contrasts]]`` and ``[de]``. Writes under ``outputs/05_de``: one contrasts table, one
results table per contrast, a metrics table the sanity rules read, the figures and the page.
"""

import hashlib
import math
import re
from collections.abc import Mapping
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from ctx import runtime
from ctx.dataset import Dataset, DatasetError, PcaSpec, load
from ctx.de.deseq import RESULT_COLUMNS, TOP_GENES, Fit, fit, pca
from ctx.de.stats import ALPHA, classify
from ctx.gather.curation import rank_replicates
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table
from ctx.report.builder import (
    FigureBlock,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    write,
)
from ctx.viz.annotate import permanova
from ctx.viz.bar import BarData
from ctx.viz.box import box_data
from ctx.viz.save import Figure, figure_block, write_figure
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import ROLES, category_colours

MODULE = "de"
VOLCANO_STEPS = (0.01, 0.05, 0.1)
P_BINS = 20
FLAT_FROM = 0.5
MAX_VOLCANO_LABELS = 10
SEX_X = "XIST"
SEX_Y = ("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")
CALL_COLOURS = {"up": ROLES["comparison"], "down": ROLES["reference"], "ns": None}
METRIC_COLUMNS = {
    "metric": ("none", "Name the sanity rules read, e.g. permanova.<pca>.<factor>.r2."),
    "value": ("varies", "The measured value; the unit follows from the metric's name."),
    "detail": ("none", "What the value was computed on, or why it is missing."),
}
CONTRAST_COLUMNS = {
    "contrast": ("none", "Contrast name from the dataset file."),
    "design": ("none", "pydeseq2 design formula."),
    "reference": ("none", "Baseline level; LFC is tested over reference."),
    "tested": ("none", "Compared level."),
    "n_reference": ("samples", "Samples in the reference group after every step."),
    "n_tested": ("samples", "Samples in the tested group after every step."),
    "paired": ("none", "Only donors with both levels kept."),
    "dropped_runs": ("none", "Runs removed as later replicates or unpaired; ; separated."),
    "groups": ("none", "Samples by group and batch (and isolation where curated)."),
    "genes_tested": ("genes", "Genes passing the expression filter and tested."),
    "up": ("genes", f"padj < {ALPHA} and raw LFC > 0."),
    "down": ("genes", f"padj < {ALPHA} and raw LFC < 0."),
    "cooks_outlier": ("genes", "No p-value: a Cook's distance outlier."),
    "independent_filter": ("genes", "A p-value but no padj: removed by independent filtering."),
}

QUESTION = "Which genes change between two groups, and how sure can we be?"
CONCEPT = (
    "DESeq2 models each gene's counts as negative binomial, with a size factor per sample for "
    "depth and a dispersion per gene for extra variability, shared across genes of similar "
    "expression. A Wald test asks whether the fold change between the groups is zero; "
    "Benjamini-Hochberg turns thousands of p-values into a false discovery rate. The design "
    "says what else differs between samples: an individual, or a sequencing batch."
)
PRACTICE = (
    "Read the p-value histogram before the volcano: a peak near 0 on a flat floor is signal, "
    "anything else is a model problem. Compare contrasts at equal n, and put a known batch in "
    "the design instead of dropping samples. Shrunk fold changes rank genes; raw ones carry "
    "the test."
)
THINK = (
    SelfCheck(
        "The naive and the adjusted PDAC contrasts use the same samples. Why do they call "
        "different genes?",
        "PDAC samples are all in the long-read batch and most healthy ones in the short: the "
        "naive design reads the batch as disease. The adjusted design estimates the batch "
        "from the healthy samples sequenced in both, and removes it.",
    ),
    SelfCheck(
        "The paired design finds more genes than the unpaired one on the same 18 samples. "
        "Where does the power come from?",
        "Each individual is compared with themself, so the variation between individuals "
        "leaves the error term.",
    ),
    SelfCheck(
        "Why can a gene with a large fold change have padj = NA?",
        "Either a single sample is an outlier (Cook's distance), or the gene is so weakly "
        "expressed that independent filtering removed it before BH.",
    ),
)
MISTAKES = (
    "Comparing a within-batch contrast of 4 samples with a pooled one of 17 and calling the "
    "difference batch: the smaller one also has less power.",
    "Counting one individual twice, e.g. several follow-ups, as independent samples.",
    "Ranking genes for GSEA by padj: it has no sign and saturates.",
)


def newest(directory: Path, dataset: str, what: str) -> Path:
    """The newest ``<date>_<module>_<dataset>_<what>_v<N>.parquet`` by date, then version.

    Raises:
        DatasetError: If there is none.
    """
    pattern = re.compile(
        rf"(\d{{4}}-\d{{2}}-\d{{2}})_[a-z]+_{re.escape(dataset)}_{what}_v(\d+)\.parquet"
    )
    found = [
        (m.group(1), int(m.group(2)), p)
        for p in directory.glob("*.parquet")
        if (m := pattern.fullmatch(p.name))
    ]
    if not found:
        raise DatasetError(
            f"no {dataset} {what} table in {directory}; run the earlier module first"
        )
    return max(found, key=lambda x: (x[0], x[1]))[2]


def _sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _subset(
    sheet: pl.DataFrame, filters: Mapping[str, tuple[str, ...]], one_per: str | None
) -> pl.DataFrame:
    """Samples matching every filter, optionally one per value of a column (its first run)."""
    for column, allowed in filters.items():
        sheet = sheet.filter(pl.col(column).is_in(list(allowed)))
    sheet = sheet.unique("sample", keep="first", maintain_order=True)
    if one_per:
        sheet = sheet.filter(rank_replicates(sheet, (one_per,), {}) == 1)
    return sheet


def _pca_figure(
    spec: PcaSpec, sheet: pl.DataFrame, counts: pl.DataFrame, threads: int
) -> tuple[ScatterData, list[tuple[str, float, str]]]:
    """One PCA scatter, and a PERMANOVA per declared factor as metrics."""
    used = _subset(sheet, spec.filters, spec.one_per)
    result = pca(counts, used["sample"].to_list(), threads=threads)
    frame = result.scores.join(
        used.select(
            "sample",
            "donor",
            *dict.fromkeys([*spec.factors, *([spec.symbol] if spec.symbol else [])]),
        ),
        on="sample",
    )
    metrics: list[tuple[str, float, str]] = []
    first = None
    for factor in spec.factors:
        labels = frame[factor].to_list()
        if not 1 < len(set(labels)) < len(labels):
            metrics.append(
                (f"permanova.{spec.name}.{factor}.r2", math.nan, "fewer than two groups")
            )
            continue
        # Strata keep labels within an individual, so only the first factor, the one that
        # varies within individuals, gets a p-value; R squared is descriptive for every one.
        strata = frame[spec.strata].to_list() if spec.strata and first is None else None
        test = permanova(result.coords, labels, strata=strata)
        where = f"{test.n} samples, {test.groups} groups, top {TOP_GENES} genes by VST variance"
        metrics.append((f"permanova.{spec.name}.{factor}.r2", test.r2, where))
        if first is None:
            metrics.append((f"permanova.{spec.name}.{factor}.p", test.p, where))
            first = test.text(factor.replace("_", " "))
    share = ", ".join(f"PC{i + 1} {s:.0%}" for i, s in enumerate(result.share))
    caption = (
        f"{first} Blind VST, the {TOP_GENES} most variable genes after the transform; "
        f"{share} of their variance; n = {frame.height} samples."
    )
    colour = spec.factors[0]
    columns = {
        "sample": ("none", "Sample (experiment accession)."),
        "donor": ("none", "Individual."),
        "pc1": ("VST units", "Score on principal component 1."),
        "pc2": ("VST units", "Score on principal component 2."),
        **{
            c: ("none", f"{c}, curated.")
            for c in frame.columns
            if c not in {"sample", "donor", "pc1", "pc2"}
        },
    }
    data = ScatterData(
        frame=frame,
        columns=columns,
        x="pc1",
        y="pc2",
        x_label="PC1",
        y_label="PC2",
        colour=colour,
        colours=category_colours(sorted(set(frame[colour]))),
        caption=caption,
        symbol=spec.symbol,
        hover=("sample", "donor"),
        width="double",
    )
    return data, metrics


def _volcano(result: Fit, genes: pl.DataFrame) -> ScatterData:
    """Volcano of one contrast, with the padj slider."""
    frame = result.results.join(genes, on="gene_id", how="left").filter(
        pl.col("pvalue").is_not_null()
    )
    frame = frame.with_columns(
        neg_log10_p=-pl.col("pvalue").log10(),
        **{
            f"call_{a}": classify(pl.col("padj"), pl.col("log2FoldChange"), a)
            for a in VOLCANO_STEPS
        },
    )
    frame = frame.with_columns(keep=pl.col(f"call_{max(VOLCANO_STEPS)}") != "ns")
    top = set(frame.filter(pl.col("call") != "ns").sort("padj").head(MAX_VOLCANO_LABELS)["gene_id"])
    frame = frame.with_columns(
        label=pl.when(pl.col("gene_id").is_in(top))
        .then(pl.coalesce("gene_name", "gene_id"))
        .otherwise(pl.lit(""))
    )
    frame = frame.select(
        "gene_id",
        "gene_name",
        "log2FoldChange",
        "neg_log10_p",
        "padj",
        *(f"call_{a}" for a in VOLCANO_STEPS),
        "keep",
        "label",
    )
    c = result.contrast
    columns = {
        "gene_id": RESULT_COLUMNS["gene_id"],
        "gene_name": ("none", "Gene symbol."),
        "log2FoldChange": RESULT_COLUMNS["log2FoldChange"],
        "padj": RESULT_COLUMNS["padj"],
        "neg_log10_p": ("-log10 p", "Minus log10 of the Wald p-value."),
        **{f"call_{a}": ("none", f"up, down or ns at padj < {a}.") for a in VOLCANO_STEPS},
        "keep": ("none", "Called at the loosest threshold; always drawn on the page."),
        "label": ("none", "Symbol printed beside the point; empty for none."),
    }
    dropped = result.results.filter(pl.col("pvalue").is_null()).height
    caption = (
        f"{c.name}: DESeq2 Wald test of {c.tested} over {c.reference}, design {c.design}; "
        f"n = {result.samples.filter(pl.col('group') == 'reference').height} {c.reference}, "
        f"{result.samples.filter(pl.col('group') == 'tested').height} {c.tested}; BH over "
        f"{result.results['padj'].is_not_null().sum():,} tested genes; colours at padj < {ALPHA}. "
        f"{dropped} Cook's outliers have no p-value and are not drawn."
    )
    return ScatterData(
        frame=frame,
        columns=columns,
        x="log2FoldChange",
        y="neg_log10_p",
        x_label=f"Fold change, {c.tested} over {c.reference}",
        y_label="Significance",
        colour=f"call_{ALPHA}",
        colours=CALL_COLOURS,
        caption=caption,
        label="label",
        hover=("gene_name", "padj"),
        width="double",
        steps=tuple((f"padj < {a}", f"call_{a}") for a in VOLCANO_STEPS),
        keep="keep",
    )


def _ma(result: Fit, genes: pl.DataFrame) -> ScatterData:
    """Mean against fold change, raw and shrunk as two series."""
    base = result.results.join(genes, on="gene_id", how="left").filter(pl.col("baseMean") > 0)
    frame = pl.concat(
        [
            base.select(
                "gene_id",
                "gene_name",
                "baseMean",
                "padj",
                "call",
                pl.col(column).alias("lfc"),
                pl.lit(label).alias("estimate"),
            )
            for column, label in (("log2FoldChange", "raw"), ("log2FoldChange_shrunk", "shrunk"))
        ]
    ).with_columns(keep=pl.col("call") != "ns")
    columns = {
        "gene_id": RESULT_COLUMNS["gene_id"],
        "gene_name": ("none", "Gene symbol."),
        "baseMean": RESULT_COLUMNS["baseMean"],
        "padj": RESULT_COLUMNS["padj"],
        "call": RESULT_COLUMNS["call"],
        "lfc": ("log2 fold change", "Raw or shrunk LFC."),
        "estimate": ("none", "raw or shrunk."),
        "keep": ("none", "Called; always drawn."),
    }
    caption = (
        f"{result.contrast.name}: raw and shrunk LFC against mean expression for "
        f"{base.height:,} genes. Shrinkage pulls noisy, weakly expressed genes toward "
        "0; the test uses the raw estimate."
    )
    return ScatterData(
        frame=frame,
        columns=columns,
        x="baseMean",
        y="lfc",
        x_label="Mean expression",
        y_label="Fold change",
        colour="estimate",
        colours=category_colours(["raw", "shrunk"]),
        caption=caption,
        hover=("gene_name", "padj"),
        log_x=True,
        width="double",
        keep="keep",
    )


def p_histogram(result: Fit) -> tuple[BarData, float, int]:
    """Twenty p-value bins; the flat level is the mean count per bin over [0.5, 1]."""
    p = result.results["pvalue"].drop_nulls().to_numpy()
    counts, edges = np.histogram(p, bins=P_BINS, range=(0.0, 1.0))
    flat = float(counts[int(FLAT_FROM * P_BINS) :].mean())
    low = int((p < ALPHA).sum())
    frame = pl.DataFrame(
        {
            "bin": [f"{a:.2f}-{b:.2f}" for a, b in zip(edges[:-1], edges[1:], strict=True)],
            "genes": counts.astype(np.int64),
            "series": ["p-values"] * P_BINS,
        }
    )
    columns = {
        "bin": ("dimensionless", "Range of Wald p-values."),
        "genes": ("genes", "Genes in the bin."),
        "series": ("none", "One series."),
    }
    caption = (
        f"{result.contrast.name}: Wald p-values of {len(p):,} tested genes. {low:,} have "
        f"p < {ALPHA} against a flat level of {flat:,.1f} per bin, the mean over "
        f"[{FLAT_FROM}, 1]. A peak at 0 over a flat floor is signal; a hill elsewhere "
        "is a model problem."
    )
    data = BarData(
        frame=frame,
        columns=columns,
        category="bin",
        value="genes",
        category_label="p-value",
        value_label="Genes",
        series="series",
        colours=category_colours(["p-values"]),
        caption=caption,
        hline=(flat, "flat level"),
        width="double",
    )
    return data, flat, low


def _dispersion(result: Fit) -> ScatterData:
    """Gene-wise and final dispersions against mean expression, with the fitted trend."""
    d = result.dispersions.filter(pl.col("baseMean") > 0).drop_nans()
    frame = (
        pl.concat(
            [
                d.select(
                    "gene_id",
                    "baseMean",
                    pl.col(column).alias("dispersion"),
                    pl.lit(label).alias("estimate"),
                )
                for column, label in (("genewise", "gene-wise"), ("final", "final"))
            ]
        )
        .filter(pl.col("dispersion") > 0)
        .with_columns(keep=pl.lit(False))
    )
    trend = d.sort("baseMean")
    line = Line(
        "fitted trend", tuple(trend["baseMean"].to_list()), tuple(trend["fitted"].to_list())
    )
    columns = {
        "gene_id": RESULT_COLUMNS["gene_id"],
        "baseMean": RESULT_COLUMNS["baseMean"],
        "dispersion": ("dimensionless", "Negative binomial dispersion."),
        "estimate": ("none", "gene-wise or final (shrunk toward the trend)."),
        "keep": ("none", "Never kept; the page thins these points."),
    }
    caption = (
        f"{result.contrast.name}: dispersion estimates of {d.height:,} genes; the final "
        "values are pulled toward the fitted trend, most for weakly expressed genes."
    )
    return ScatterData(
        frame=frame,
        columns=columns,
        x="baseMean",
        y="dispersion",
        x_label="Mean expression",
        y_label="Dispersion",
        colour="estimate",
        colours=category_colours(["gene-wise", "final"]),
        caption=caption,
        lines=(line,),
        hover=("gene_id",),
        log_x=True,
        log_y=True,
        width="double",
        keep="keep",
    )


def _calls(title: str, names: tuple[str, ...], fits: Mapping[str, Fit]) -> BarData:
    """Up and down calls per contrast, for one lesson."""
    rows = [
        (n, call, int((fits[n].results["call"] == call).sum()))
        for n in names
        for call in ("up", "down")
    ]
    frame = pl.DataFrame(rows, schema=["contrast", "call", "genes"], orient="row")
    sizes = "; ".join(
        f"{n}: n = {fits[n].contrast.expected_n[0]} vs {fits[n].contrast.expected_n[1]}, "
        f"{fits[n].contrast.design}"
        for n in names
    )
    columns = {
        "contrast": ("none", "Contrast."),
        "call": ("none", "up or down."),
        "genes": ("genes", f"Genes at padj < {ALPHA}."),
    }
    return BarData(
        frame=frame,
        columns=columns,
        category="contrast",
        value="genes",
        category_label=title,
        value_label="Genes called",
        series="call",
        colours={"up": ROLES["comparison"], "down": ROLES["reference"]},
        caption=f"{title}: genes at padj < {ALPHA}. {sizes}.",
        horizontal=True,
        stacked=True,
        width="double",
    )


def _lfc_scatter(x: Fit, y: Fit, genes: pl.DataFrame) -> tuple[ScatterData, float]:
    """One contrast's raw LFC against another's, on genes tested in both."""
    frame = (
        x.results.select(
            "gene_id", pl.col("log2FoldChange").alias("lfc_x"), pl.col("call").alias("call_x")
        )
        .join(
            y.results.select(
                "gene_id", pl.col("log2FoldChange").alias("lfc_y"), pl.col("call").alias("call_y")
            ),
            on="gene_id",
        )
        .drop_nulls()
        .join(genes, on="gene_id", how="left")
        .with_columns(
            called=pl.when((pl.col("call_x") != "ns") | (pl.col("call_y") != "ns"))
            .then(pl.lit("called in either"))
            .otherwise(pl.lit("neither"))
        )
        .with_columns(keep=pl.col("called") == "called in either")
    )
    rho = float(spearmanr(frame["lfc_x"], frame["lfc_y"]).statistic)
    shared = set(x.samples["donor"]) & set(y.samples["donor"])
    note = (
        f" {len(shared)} individuals are in both contrasts, so the two estimates are not "
        "independent."
        if shared
        else " No individual is in both contrasts."
    )
    ends = (
        float(frame.select(pl.min_horizontal("lfc_x", "lfc_y").min()).item()),
        float(frame.select(pl.max_horizontal("lfc_x", "lfc_y").max()).item()),
    )
    columns = {
        "gene_id": RESULT_COLUMNS["gene_id"],
        "gene_name": ("none", "Gene symbol."),
        "lfc_x": ("log2 FC", f"Raw log2 fold change in {x.contrast.name}."),
        "lfc_y": ("log2 FC", f"Raw log2 fold change in {y.contrast.name}."),
        "call_x": ("none", f"Call in {x.contrast.name}."),
        "call_y": ("none", f"Call in {y.contrast.name}."),
        "called": ("none", "Called in either contrast, or neither."),
        "keep": ("none", "Called in either; always drawn on the page."),
    }
    caption = (
        f"Raw LFC of {frame.height:,} genes tested in both: Spearman rho = {rho:.2f}; the "
        f"line is x = y.{note}"
    )
    data = ScatterData(
        frame=frame,
        columns=columns,
        x="lfc_x",
        y="lfc_y",
        x_label=x.contrast.name,
        y_label=y.contrast.name,
        colour="called",
        colours={"called in either": ROLES["comparison"], "neither": None},
        caption=caption,
        hover=("gene_name",),
        lines=(Line("x = y", ends, ends),),
        width="double",
        keep="keep",
    )
    return data, rho


def _box(
    result: Fit, counts: pl.DataFrame, genes: pl.DataFrame, symbols: tuple[str, ...]
) -> tuple[Figure | None, list[str]]:
    """Normalised expression of named genes in one contrast, MWU in the caption, padj on hover."""
    named = genes.filter(pl.col("gene_name").is_in(list(symbols))).unique("gene_name", keep="first")
    missing = [s for s in symbols if s not in set(named["gene_name"])]
    if named.is_empty():
        return None, missing
    c = result.contrast
    labels = {"reference": c.reference, "tested": c.tested}
    padj = dict(result.results.select("gene_id", "padj").rows())
    sf = dict(result.size_factors.rows())
    rows = []
    for gene_id, symbol in named.select("gene_id", "gene_name").rows():
        values = counts.filter(pl.col("gene_id") == gene_id)
        for sample, group, donor in result.samples.select("sample", "group", "donor").rows():
            level = float(values[sample][0]) / sf[sample]
            rows.append(
                (symbol, labels[group], sample, donor, math.log2(level + 1), padj.get(gene_id))
            )
    frame = pl.DataFrame(
        rows, schema=["feature", "group", "sample", "donor", "value", "deseq2_padj"], orient="row"
    )
    columns = {
        "feature": ("none", "Gene symbol."),
        "group": ("none", "Compared level."),
        "sample": ("none", "Sample."),
        "donor": ("none", "Individual."),
        "value": ("log2 reads", "log2(count / DESeq2 size factor + 1): normalised."),
        "deseq2_padj": ("dimensionless", "DESeq2 padj of the gene in this contrast."),
    }
    data = box_data(
        frame,
        columns=columns,
        feature_label="Gene",
        value_label="Expression",
        groups=(c.reference, c.tested),
        hover=("sample", "donor", "deseq2_padj"),
    )
    return data, missing


def sex_check(
    sheet: pl.DataFrame, counts: pl.DataFrame, genes: pl.DataFrame
) -> tuple[pl.DataFrame | None, list[tuple[str, float, str]]]:
    """Sex from XIST against five chrY genes, per sample, against the recorded sex.

    A transplanted sample is compared with its donor's sex, where one is recorded: myeloid
    cells after HSCT come from the donor.
    """
    if "sex" not in sheet.columns:
        return None, [("sex.agree", math.nan, "not computed: the sheet has no sex column")]
    ids = dict(genes.select("gene_name", "gene_id").unique("gene_name", keep="first").rows())
    if missing := [g for g in (SEX_X, *SEX_Y) if g not in ids]:
        return None, [("sex.agree", math.nan, f"not computed: no gene named {missing}")]
    samples = sheet.unique("sample", keep="first", maintain_order=True)
    names = samples["sample"].to_list()
    matrix = counts.select(names).to_numpy().astype(float)
    cpm = matrix / matrix.sum(axis=0) * 1e6
    row = {g: i for i, g in enumerate(counts["gene_id"].to_list())}
    xist = np.log2(cpm[row[ids[SEX_X]]] + 1)
    chr_y = np.log2(cpm[[row[ids[g]] for g in SEX_Y]] + 1).mean(axis=0)
    donor = samples["donor_sex"] if "donor_sex" in samples.columns else pl.Series([""] * len(names))
    expected = [
        d if d else s
        for s, d in zip(samples["sex"].fill_null(""), donor.fill_null(""), strict=True)
    ]
    frame = pl.DataFrame(
        {
            "sample": names,
            "expected": expected,
            "xist_log2cpm": xist,
            "chry_log2cpm": chr_y,
            "expression_sex": np.where(chr_y > xist, "M", "F"),
        }
    ).with_columns(
        agrees=pl.when(pl.col("expected") == "")
        .then(None)
        .otherwise(pl.col("expected") == pl.col("expression_sex"))
    )
    compared = frame.filter(pl.col("agrees").is_not_null())
    wrong = compared.filter(~pl.col("agrees"))["sample"].to_list()
    unrecorded = frame.filter(pl.col("agrees").is_null())["sample"].to_list()
    return frame, [
        (
            "sex.compared",
            float(compared.height),
            f"{len(unrecorded)} with no recorded sex: {', '.join(unrecorded)}",
        ),
        ("sex.disagree", float(len(wrong)), ", ".join(wrong) or "none"),
    ]


def _contrast_row(result: Fit, sheet_columns: list[str]) -> dict[str, object]:
    """One row of the contrasts table."""
    c, r = result.contrast, result.results
    keys = [k for k in ("group", "batch", "isolation") if k in sheet_columns or k == "group"]
    groups = "; ".join(
        f"{'/'.join(str(v) for v in key)}: {n}"
        for *key, n in result.samples.group_by(keys).len().sort(keys).rows()
    )
    reasons = dict(r["na_reason"].value_counts().rows())
    return {
        "contrast": c.name,
        "design": c.design,
        "reference": c.reference,
        "tested": c.tested,
        "n_reference": c.expected_n[0],
        "n_tested": c.expected_n[1],
        "paired": c.paired,
        "dropped_runs": ";".join(result.dropped),
        "groups": groups,
        "genes_tested": r.height,
        "up": int((r["call"] == "up").sum()),
        "down": int((r["call"] == "down").sum()),
        "cooks_outlier": reasons.get("cooks_outlier", 0),
        "independent_filter": reasons.get("independent_filter", 0),
    }


GENE_NAME = ("none", "Gene symbol.")
SEX_COLUMNS = {
    "sample": ("none", "Sample."),
    "expected": ("none", "Recorded sex; the donor's after HSCT."),
    "xist_log2cpm": ("log2 CPM", "XIST expression."),
    "chry_log2cpm": ("log2 CPM", f"Mean of {', '.join(SEX_Y)}."),
    "expression_sex": ("none", "M when chrY genes exceed XIST, else F."),
    "agrees": ("none", "Expression sex equals the expected; empty when none recorded."),
}


class Writer:
    """Where one module run writes: its tables, figures and page blocks, under one stem rule."""

    def __init__(self, out: Path, dataset: str, day: date, releases: dict[str, str]) -> None:
        """Remember the output directory, naming and provenance every file shares."""
        self.out, self.dataset, self.day, self.releases = out, dataset, day, releases
        self.files: list[Path] = []
        self.blocks: list[FigureBlock | str] = []

    def stem(self, what: str) -> str:
        """The stem of one output of this run."""
        return stem(self.out.name.split("_", 1)[1], f"{self.dataset}_{what}", 1, day=self.day)

    def table(self, frame: pl.DataFrame, what: str, columns: Mapping[str, tuple[str, str]]) -> None:
        """Write one table with its README."""
        self.files += write_table(
            frame, self.out, self.stem(what), columns=columns, data_releases=self.releases
        )

    def figure(self, data: Figure, what: str, title: str) -> None:
        """Write one figure and put it on the page."""
        figure_stem = self.stem(what)
        self.files += write_figure(
            data, self.out, figure_stem, title=title, data_releases=self.releases
        )
        self.blocks.append(figure_block(data, figure_stem))


def _fit_metrics(fits: Mapping[str, Fit]) -> list[tuple[str, float, str]]:
    """p-value histogram and call counts of every contrast, for the sanity rules."""
    metrics = []
    for name, result in fits.items():
        _, flat, low = p_histogram(result)
        tested = float(result.results["padj"].is_not_null().sum())
        metrics += [
            (f"phist.{name}.low", float(low), f"p < {ALPHA}"),
            (f"phist.{name}.flat", flat, f"mean count per bin over [{FLAT_FROM}, 1]"),
            (
                f"calls.{name}.called",
                float((result.results["call"] != "ns").sum()),
                f"padj < {ALPHA}",
            ),
            (f"calls.{name}.tested", tested, "padj present"),
        ]
    return metrics


def _figures(
    writer: Writer,
    spec: Dataset,
    fits: Mapping[str, Fit],
    data: tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame],
    threads: int,
) -> list[tuple[str, float, str]]:
    """Every M5 figure, in page order; returns the metrics they measured."""
    sheet, counts, genes = data
    assert spec.de is not None  # noqa: S101 - run() refuses a dataset file without [de]
    metrics = []
    for pca_spec in spec.de.pca:
        scatter, found = _pca_figure(pca_spec, sheet, counts, threads)
        metrics += found
        writer.figure(scatter, f"pca_{pca_spec.name}", f"M5 PCA {pca_spec.name}")
    featured = fits[spec.de.featured]
    name = featured.contrast.name
    writer.figure(_volcano(featured, genes), f"volcano_{name}", "M5 volcano")
    writer.figure(_ma(featured, genes), f"ma_{name}", "M5 MA plot")
    writer.figure(p_histogram(featured)[0], f"phist_{name}", "M5 p-value histogram")
    writer.figure(_dispersion(featured), f"dispersion_{name}", "M5 dispersion")
    for title, names in spec.de.lessons.items():
        slug = re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_")
        writer.figure(_calls(title, names, fits), f"calls_{slug}", f"M5 {title}")
    for x, y in spec.de.scatters:
        scatter, rho = _lfc_scatter(fits[x], fits[y], genes)
        metrics.append((f"spearman.{x}.{y}", rho, "raw LFC, genes tested in both"))
        writer.figure(scatter, f"lfc_{x}_vs_{y}", f"M5 {x} against {y}")
    if spec.de.box_contrast:
        box, missing = _box(fits[spec.de.box_contrast], counts, genes, spec.de.box_genes)
        if box is not None:
            writer.figure(box, f"box_{spec.de.box_contrast}", "M5 named genes")
        if missing:
            writer.blocks.append(f"Named genes not in the annotation: {', '.join(missing)}.")
    return metrics + _fit_metrics(fits)


def run(  # noqa: PLR0913, PLR0915 - one call from the CLI; the page is one long sequence
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset: str,
    *,
    dataset_file: Path,
    min_count: int,
) -> list[Path]:
    """Fit the dataset's contrasts and write the M5 tables, figures and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: Slug naming the dataset.
        dataset_file: The dataset file, with ``[[contrasts]]`` and ``[de]``.
        min_count: The expression filter's reads at the median library size.

    Returns:
        The page, then every table and figure file.

    Raises:
        DatasetError: If the dataset file has no contrasts, or gather has not been run.
    """
    spec = load(dataset_file)
    if spec.slug != dataset or spec.de is None or not spec.contrasts:
        raise DatasetError(f"{dataset_file}: needs slug {dataset!r}, [[contrasts]] and [de]")
    day = date.today()
    gathered = settings.root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    inputs = {what: newest(gathered, dataset, what) for what in ("samplesheet", "counts", "genes")}
    sheet, counts, genes = (pl.read_parquet(inputs[w]) for w in ("samplesheet", "counts", "genes"))
    releases = {
        "dataset file": str(spec.path),
        f"sha256:{spec.path.name}": _sha256(spec.path),
        **{f"sha256:{p.name}": _sha256(p) for p in inputs.values()},
    }
    fits = {
        c.name: fit(counts, sheet, c, spec.order, threads=settings.threads, min_count=min_count)
        for c in spec.contrasts
    }
    writer = Writer(
        settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset, day, releases
    )
    writer.blocks.append(
        f"Inputs: {', '.join(p.name for p in inputs.values())} (newest gather outputs). "
        f"{len(fits)} contrasts fitted with DESeq2 (pydeseq2), BH at padj < {ALPHA}; the "
        "contrasts table gives each design, n, dropped runs and calls."
    )
    rows = [_contrast_row(f, sheet.columns) for f in fits.values()]
    writer.table(pl.DataFrame(rows), "contrasts", CONTRAST_COLUMNS)
    for name, result in fits.items():
        frame = result.results.join(genes, on="gene_id", how="left")
        writer.table(frame, f"{name}_results", {**RESULT_COLUMNS, "gene_name": GENE_NAME})
    metrics = _figures(writer, spec, fits, (sheet, counts, genes), settings.threads)
    sexes, found = sex_check(sheet, counts, genes)
    if sexes is not None:
        writer.table(sexes, "sex_check", SEX_COLUMNS)
    frame = pl.DataFrame(metrics + found, schema=["metric", "value", "detail"], orient="row")
    writer.table(frame, "metrics", METRIC_COLUMNS)
    out, files, blocks = writer.out, writer.files, writer.blocks
    usage = runtime.measure(settings.started)
    budget = runtime.budget(MODULE)
    bound = usage.peak_rss_self + settings.threads * usage.peak_rss_children
    page = ModulePage(
        title=f"M5 Differential expression: {dataset}",
        question=QUESTION,
        concept=CONCEPT,
        parameters=cards,
        results=tuple(blocks),
        practice=PRACTICE,
        think=THINK,
        mistakes=MISTAKES,
        reproduce=Reproduction(
            command=command,
            wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=(
                f"Track B: under {budget.minutes} min and {budget.ram_gb} GB; RAM bound "
                f"self + {settings.threads} x largest child = {bound:,} bytes"
            ),
            versions=runtime.versions(),
        ),
    )
    return [write(page, out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"), *files]
