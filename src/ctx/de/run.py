"""M5: ``ctx de``. Fit every declared contrast, then write the tables, the figures and the page.

Reads the newest gather outputs of one dataset (sheet, counts, genes) and the dataset file's
``[[contrasts]]`` and ``[de]``. Writes under ``outputs/05_de``: one contrasts table, one
results table per contrast, a metrics table the sanity rules read, the figures and the page.
"""

import hashlib
import math
import re
from collections.abc import Mapping
from dataclasses import replace
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from ctx import runtime
from ctx.dataset import Dataset, DatasetError, PcaSpec, load
from ctx.de.deseq import RESULT_COLUMNS, TOP_GENES, Fit, fit, pca
from ctx.de.page import (
    BASES,
    CONCEPT,
    LINE_RULE,
    MISTAKES,
    NO_LINE_RULE,
    PRACTICE,
    QUESTION,
    RAW_WORDS,
    TERMS,
    VOLCANO_READING,
    VOLCANO_SUMMARY,
    box_text,
    calls_text,
    dispersion_text,
    lfc_text,
    ma_text,
    pca_text,
    phist_text,
    think,
    volcano_text,
)
from ctx.de.stats import ADJUSTMENTS, ALPHA, LFC_MIN, P_FLOOR, classify, p_line
from ctx.gather.curation import rank_replicates
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table
from ctx.report.builder import (
    FigureBlock,
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    write,
)
from ctx.report.text import fill, word
from ctx.viz.annotate import Permanova, permanova
from ctx.viz.bar import BarData
from ctx.viz.box import box_data
from ctx.viz.save import Figure, figure_block, write_figure
from ctx.viz.scatter import CUTOFF_PRESETS, Calls, Line, ScatterData
from ctx.viz.style import ROLES, category_colours

MODULE = "de"
# The content file holding the page's sentences moved so far (P2S.4 moves the rest).
PAGE = "de"
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
    "genes_fitted": ("genes", "Genes past the expression filter, fitted by DESeq2."),
    "genes_tested": (
        "genes",
        "Genes with an adjusted p: fitted, less Cook's outliers and independent filtering.",
    ),
    "up": ("genes", f"padj < {ALPHA} and raw LFC >= {LFC_MIN}."),
    "down": ("genes", f"padj < {ALPHA} and raw LFC <= -{LFC_MIN}."),
    "cooks_outlier": ("genes", "No p-value: a Cook's distance outlier."),
    "independent_filter": ("genes", "A p-value but no padj: removed by independent filtering."),
}


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


def n_tested(results: pl.DataFrame) -> int:
    """Genes tested: those with an adjusted p, after Cook's and independent filtering.

    The one definition every page count, metric and caption uses.
    """
    return int(results["padj"].is_not_null().sum())


def _sizes(result: Fit) -> str:
    """Group sizes as the page prints them everywhere: tested group first, each named."""
    c = result.contrast
    n = {g: result.samples.filter(pl.col("group") == g).height for g in ("tested", "reference")}
    return f"{n['tested']} {c.tested} vs {n['reference']} {c.reference}"


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
) -> tuple[ScatterData, list[tuple[str, float, str]], FigureText]:
    """One PCA scatter, a PERMANOVA per declared factor as metrics, and its page words."""
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
    first: str | None = None
    shown: tuple[str, Permanova] | None = None
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
            first, shown = test.text(factor.replace("_", " ")), (factor, test)
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
    # The page describes the first factor's test, the one with a p-value.
    text = pca_text(spec.name, spec.filters, colour, spec.symbol, result.share, shown)
    return data, metrics, text


def _volcano(result: Fit, genes: pl.DataFrame) -> ScatterData:
    """Volcano of one contrast: the default call drawn, and every cutoff callable live.

    Every gene some setting of the page could call is drawn; the rest, above the loosest
    cutoff offered or never tested, are thinned. All four adjustments count over the same
    tested genes, so the live view and the table agree.
    """
    # Only tested genes: the others have no adjusted p, so no setting of the page can call
    # them, and a dot the rule cannot reach would contradict the reading guide.
    frame = result.results.join(genes, on="gene_id", how="left").filter(
        pl.col("padj").is_not_null()
    )
    line = p_line("padj", ALPHA, frame.height, int((frame["padj"] < ALPHA).sum()))
    frame = frame.with_columns(
        gene=pl.coalesce("gene_name", "gene_id"),
        neg_log10_p=-pl.max_horizontal(pl.col("pvalue"), pl.lit(P_FLOOR)).log10(),
        keep=pl.col("pvalue") < max(CUTOFF_PRESETS),
    )
    top = set(frame.filter(pl.col("call") != "ns").sort("padj").head(MAX_VOLCANO_LABELS)["gene_id"])
    frame = frame.with_columns(
        label=pl.when(pl.col("gene_id").is_in(top)).then(pl.col("gene")).otherwise(pl.lit(""))
    ).select(
        "gene", "gene_id", "log2FoldChange", "pvalue", *ADJUSTMENTS, "neg_log10_p", "call",
        "keep", "label",
    )  # fmt: skip
    c = result.contrast
    columns = {
        "gene": ("none", "Gene symbol, or the identifier where there is none."),
        "gene_id": RESULT_COLUMNS["gene_id"],
        "log2FoldChange": RESULT_COLUMNS["log2FoldChange"],
        "pvalue": RESULT_COLUMNS["pvalue"],
        **{column: RESULT_COLUMNS[column] for column in ADJUSTMENTS},
        "neg_log10_p": ("-log10 p", f"Minus log10 of the raw p-value, drawn at {P_FLOOR:g} "
                                    "at least."),
        "call": RESULT_COLUMNS["call"],
        "keep": ("none", f"Raw p below {max(CUTOFF_PRESETS)}: callable at some setting, so "
                         "always drawn."),
        "label": ("none", "Symbol printed beside the point; empty for none."),
    }  # fmt: skip
    untested = result.results.height - frame.height
    floored = int((frame["pvalue"] < P_FLOOR).sum())
    caption = (
        f"{c.name}: DESeq2 Wald test of {c.tested} over {c.reference}, design {c.design}; "
        f"n = {_sizes(result)}. {fill(PAGE, 'figures.volcano.tested', tested=frame.height)} "
        f"Default call: BH-adjusted p below {ALPHA} and log2 FC at or beyond ±{LFC_MIN:g}."
        + (f" {untested:,} genes with no adjusted p (Cook's outliers or independent "
           "filtering) are not drawn." if untested else "")
        + (f" {floored} p-values below {P_FLOOR:g} are drawn at it." if floored else "")
    )  # fmt: skip
    calls = Calls(
        lfc="log2FoldChange", p="pvalue", adjusted=BASES, raw_words=RAW_WORDS, tested=c.tested,
        reference=c.reference, alpha=ALPHA, lfc_min=LFC_MIN, max_cutoff=max(CUTOFF_PRESETS),
        n_tested=frame.height, p_floor=P_FLOOR, p_line=line, reading=VOLCANO_READING,
        summary=VOLCANO_SUMMARY, line_rule=LINE_RULE, no_line_rule=NO_LINE_RULE,
    )  # fmt: skip
    return ScatterData(
        frame=frame,
        columns=columns,
        x="log2FoldChange",
        y="neg_log10_p",
        x_label=f"log2 fold change, {c.tested} over {c.reference}",
        y_label="-log10 raw p",
        colour="call",
        colours=CALL_COLOURS,
        caption=caption,
        label="label",
        hover=("gene", "gene_id", "padj"),
        width="double",
        keep="keep",
        names={"up": f"up in {c.tested}", "down": f"down in {c.tested}", "ns": "not called"},
        calls=calls,
    )


def _ma(result: Fit, genes: pl.DataFrame) -> ScatterData:
    """Mean against fold change, raw and shrunk as two series; raw alone without shrinkage.

    A page that thins the genes not called draws both estimates of a gene or neither, so the
    pull of shrinkage stays visible gene by gene; genes in order of their raw estimate.
    """
    base = result.results.join(genes, on="gene_id", how="left").filter(pl.col("baseMean") > 0)
    shrunk = base["log2FoldChange_shrunk"].is_not_null().any()
    estimates = (("log2FoldChange", "raw"), ("log2FoldChange_shrunk", "shrunk"))[: 1 + shrunk]
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
            for column, label in estimates
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
        "keep": ("none", "Called at the default rule; always drawn."),
    }
    caption = (
        f"{result.contrast.name}: raw and shrunk LFC against mean expression for "
        f"{base.height:,} genes. Shrinkage pulls noisy, weakly expressed genes toward "
        "0; the test uses the raw estimate."
        if shrunk
        else f"{result.contrast.name}: raw LFC against mean expression for {base.height:,} "
        "genes; shrinkage was off (--no-lfc-shrink). The test uses the raw estimate."
    )
    return ScatterData(
        frame=frame,
        columns=columns,
        x="baseMean",
        y="lfc",
        x_label="Mean expression",
        y_label="log2 fold change",
        colour="estimate",
        colours=category_colours([label for _, label in estimates]),
        names={"raw": "raw estimate", "shrunk": "shrunk estimate"},
        caption=caption,
        hover=("gene_name", "padj"),
        log_x=True,
        width="double",
        keep="keep",
        unit="gene_id",
        thinned=word(PAGE, "figures.ma.thinned"),
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
        f"{result.contrast.name}: Wald p-values of the {len(p):,} genes with one, in "
        f"{P_BINS} bins; {low:,} have p < {ALPHA}. The flat level, {flat:,.1f}, is the mean "
        f"count per bin over [{FLAT_FROM}, 1]."
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
    """Gene-wise and final dispersions against mean expression, with the fitted trend.

    A page that thins them draws both of a gene's estimates or neither, as the MA plot does.
    """
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
        f"{result.contrast.name}: gene-wise and final dispersion estimates of {d.height:,} "
        "genes, with the trend DESeq2 fitted across them."
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
        unit="gene_id",
        thinned=word(PAGE, "figures.dispersion.thinned"),
    )


def _calls(title: str, names: tuple[str, ...], fits: Mapping[str, Fit]) -> BarData:
    """Up and down calls per contrast, for one lesson."""
    rows = [
        (n, call, int((fits[n].results["call"] == call).sum()))
        for n in names
        for call in ("up", "down")
    ]
    frame = pl.DataFrame(rows, schema=["contrast", "call", "genes"], orient="row")
    sizes = "; ".join(f"{n}: n = {_sizes(fits[n])}, {fits[n].contrast.design}" for n in names)
    columns = {
        "contrast": ("none", "Contrast."),
        "call": ("none", "up or down."),
        "genes": ("genes", f"Genes called: padj < {ALPHA} and log2 FC at or beyond ±{LFC_MIN:g}."),
    }
    return BarData(
        frame=frame,
        columns=columns,
        category="contrast",
        value="genes",
        category_label="Contrast",
        value_label="Genes called",
        series="call",
        colours={"up": ROLES["comparison"], "down": ROLES["reference"]},
        caption=(
            f"{title}: genes called per contrast (BH-adjusted p < {ALPHA}, log2 FC at or beyond ±"
            f"{LFC_MIN:g}). {sizes}."
        ),
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
        .join(
            x.results.filter(pl.col("padj").is_not_null())
            .select("gene_id")
            .join(y.results.filter(pl.col("padj").is_not_null()).select("gene_id"), on="gene_id"),
            on="gene_id",
        )
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
        x_label=f"log2 FC, {x.contrast.name}",
        y_label=f"log2 FC, {y.contrast.name}",
        colour="called",
        colours={"called in either": ROLES["comparison"], "neither": None},
        caption=caption,
        hover=("gene_name",),
        lines=(Line("x = y", ends, ends),),
        width="double",
        keep="keep",
        thinned=word(PAGE, "figures.lfc.thinned"),
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
        "genes_fitted": r.height,
        "genes_tested": n_tested(r),
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

    def stem(self, what: str) -> str:
        """The stem of one output of this run."""
        return stem(self.out.name.split("_", 1)[1], f"{self.dataset}_{what}", 1, day=self.day)

    def table(self, frame: pl.DataFrame, what: str, columns: Mapping[str, tuple[str, str]]) -> None:
        """Write one table with its README."""
        self.files += write_table(
            frame, self.out, self.stem(what), columns=columns, data_releases=self.releases
        )

    def figure(self, data: Figure, what: str, text: FigureText) -> FigureBlock:
        """Write one figure, its page linked back to the module page; return its page block."""
        figure_stem = self.stem(what)
        text = replace(text, back=f"{self.stem('overview')}.html")
        self.files += write_figure(
            data, self.out, figure_stem, title=text.title, data_releases=self.releases, text=text
        )
        return figure_block(data, figure_stem, text)


def _fit_metrics(fits: Mapping[str, Fit]) -> list[tuple[str, float, str]]:
    """p-value histogram and call counts of every contrast, for the sanity rules."""
    metrics = []
    for name, result in fits.items():
        _, flat, low = p_histogram(result)
        r = result.results
        # The rules on the model's false positives and power read padj alone: a fold-change
        # filter would flatter both, so lfc_min 0 counts every significant gene.
        any_fold = classify(pl.col("padj"), pl.col("log2FoldChange"), lfc_min=0.0)
        metrics += [
            (f"phist.{name}.low", float(low), f"p < {ALPHA}"),
            (f"phist.{name}.flat", flat, f"mean count per bin over [{FLAT_FROM}, 1]"),
            (f"calls.{name}.significant", float(r.select(any_fold != "ns").sum().item()),
             f"padj < {ALPHA}, any fold change"),
            (f"calls.{name}.up", float((r["call"] == "up").sum()),
             f"padj < {ALPHA} and raw LFC >= {LFC_MIN}"),
            (f"calls.{name}.down", float((r["call"] == "down").sum()),
             f"padj < {ALPHA} and raw LFC <= -{LFC_MIN}"),
            (f"calls.{name}.tested", float(n_tested(r)), "genes with an adjusted p"),
        ]  # fmt: skip
    return metrics


def _lessons(
    writer: Writer, spec: Dataset, fits: Mapping[str, Fit], genes: pl.DataFrame
) -> tuple[list[Result], list[tuple[str, float, str]]]:
    """One group per lesson and one for the scatters no lesson holds; rho as metrics.

    A lesson shows its call counts, then the fold-change scatters of its own contrasts.
    """
    assert spec.de is not None  # noqa: S101 - run() refuses a dataset file without [de]
    groups: list[Result] = []
    metrics: list[tuple[str, float, str]] = []
    placed: dict[tuple[str, str], FigureBlock] = {}
    for x, y in spec.de.scatters:
        scatter, rho = _lfc_scatter(fits[x], fits[y], genes)
        metrics.append((f"spearman.{x}.{y}", rho, "raw LFC, genes tested in both"))
        placed[(x, y)] = writer.figure(scatter, f"lfc_{x}_vs_{y}", lfc_text(x, y, rho))
    for title, names in spec.de.lessons.items():
        slug = re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_")
        items = [writer.figure(_calls(title, names, fits), f"calls_{slug}",
                               calls_text(title, names, fits))]  # fmt: skip
        items += [placed.pop(pair) for pair in list(placed) if set(pair) <= set(names)]
        groups.append(Group(f"Lesson: {title}", tuple(items)))
    if placed:
        groups.append(Group("Across contrasts", tuple(placed.values())))
    return groups, metrics


def _figures(
    writer: Writer,
    spec: Dataset,
    fits: Mapping[str, Fit],
    data: tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame],
    threads: int,
) -> tuple[list[Result], list[tuple[str, float, str]]]:
    """Every M5 figure, grouped as the page shows them; and the metrics they measured."""
    sheet, counts, genes = data
    assert spec.de is not None  # noqa: S101 - run() refuses a dataset file without [de]
    groups: list[Result] = []
    metrics = []
    pcas = []
    for pca_spec in spec.de.pca:
        scatter, found, text = _pca_figure(pca_spec, sheet, counts, threads)
        metrics += found
        pcas.append(writer.figure(scatter, f"pca_{pca_spec.name}", text))
    if pcas:
        groups.append(Group("Samples", tuple(pcas), lead="Do samples group by biology, or by "
                            "how they were processed?", tabs=len(pcas) > 1))  # fmt: skip
    featured = fits[spec.de.featured]
    c, name = featured.contrast, featured.contrast.name
    groups.append(Group(
        f"Featured contrast: {name}",
        (writer.figure(_volcano(featured, genes), f"volcano_{name}", volcano_text(c)),
         writer.figure(_ma(featured, genes), f"ma_{name}",
                       ma_text(c, featured.results["log2FoldChange_shrunk"].is_not_null().any()))),
        lead="Which genes changed, and by how much?", tabs=True,
    ))  # fmt: skip
    phist, flat, low = p_histogram(featured)
    with_p = featured.results["pvalue"].drop_nulls().len()
    groups.append(Group(
        "Model health",
        (writer.figure(phist, f"phist_{name}", phist_text(
            name, with_p, low, flat, int(phist.frame["genes"][0]), 1 / P_BINS)),
         writer.figure(_dispersion(featured), f"dispersion_{name}", dispersion_text(name))),
        lead="Check these before trusting any call.", tabs=True,
    ))  # fmt: skip
    lessons, found = _lessons(writer, spec, fits, genes)
    groups += lessons
    metrics += found
    if spec.de.box_contrast:
        box, missing = _box(fits[spec.de.box_contrast], counts, genes, spec.de.box_genes)
        named: list[FigureBlock | str] = []
        if box is not None:
            n_genes = len(spec.de.box_genes) - len(missing)
            named.append(
                writer.figure(
                    box,
                    f"box_{spec.de.box_contrast}",
                    box_text(fits[spec.de.box_contrast].contrast, n_genes),
                )
            )
        if missing:
            named.append(f"Named genes not in the annotation: {', '.join(missing)}.")
        groups.append(Group("Named genes", tuple(named)))
    return groups, metrics + _fit_metrics(fits)


def _contrasts_table(fits: Mapping[str, Fit]) -> ResultTable:
    """The contrasts, as the page shows them; the CSV beside it holds every column."""
    rows = []
    for name, f in fits.items():
        c, r = f.contrast, f.results
        rows.append((
            name, f"{c.tested} vs {c.reference}", c.design,
            f"{c.expected_n[1]} vs {c.expected_n[0]}", f"{n_tested(r):,}",
            f"{(r['call'] == 'up').sum():,}", f"{(r['call'] == 'down').sum():,}",
        ))  # fmt: skip
    columns = ("Contrast", "Compares", "Design", "n", "Tested", "Up", "Down")
    return ResultTable("Contrasts; tested means with an adjusted p", columns, tuple(rows),
                       numeric=(4, 5, 6))  # fmt: skip


def _keys(sheet: pl.DataFrame, fits: Mapping[str, Fit], featured: str) -> tuple[KeyNumber, ...]:
    """Samples, contrasts, and the featured contrast's tested and called genes."""
    r = fits[featured].results
    up, down = int((r["call"] == "up").sum()), int((r["call"] == "down").sum())
    return (
        KeyNumber("Samples", f"{sheet['sample'].n_unique():,}", "in the sample sheet"),
        KeyNumber("Contrasts", f"{len(fits)}"),
        KeyNumber("Genes tested", f"{n_tested(r):,}", f"in {featured}"),
        KeyNumber("Genes called", f"{up + down:,}", f"{up:,} up, {down:,} down in {featured}"),
    )


def run(  # noqa: PLR0913, PLR0915 - one call from the CLI; the page is one long sequence
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset: str,
    *,
    dataset_file: Path,
    min_count: int,
    lfc_shrink: bool,
) -> list[Path]:
    """Fit the dataset's contrasts and write the M5 tables, figures and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: Slug naming the dataset.
        dataset_file: The dataset file, with ``[[contrasts]]`` and ``[de]``.
        min_count: The expression filter's reads at the median library size.
        lfc_shrink: Add the apeGLM-shrunk fold changes; off leaves them empty.

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
        c.name: fit(
            counts,
            sheet,
            c,
            spec.order,
            threads=settings.threads,
            min_count=min_count,
            shrink=lfc_shrink,
        )  # fmt: skip
        for c in spec.contrasts
    }
    writer = Writer(
        settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset, day, releases
    )
    rows = [_contrast_row(f, sheet.columns) for f in fits.values()]
    writer.table(pl.DataFrame(rows), "contrasts", CONTRAST_COLUMNS)
    for name, result in fits.items():
        frame = result.results.join(genes, on="gene_id", how="left")
        writer.table(frame, f"{name}_results", {**RESULT_COLUMNS, "gene_name": GENE_NAME})
    groups, metrics = _figures(writer, spec, fits, (sheet, counts, genes), settings.threads)
    sexes, found = sex_check(sheet, counts, genes)
    if sexes is not None:
        writer.table(sexes, "sex_check", SEX_COLUMNS)
    frame = pl.DataFrame(metrics + found, schema=["metric", "value", "detail"], orient="row")
    writer.table(frame, "metrics", METRIC_COLUMNS)
    out, files = writer.out, writer.files
    usage = runtime.measure(settings.started)
    budget = runtime.budget(MODULE)
    bound = usage.peak_rss_self + settings.threads * usage.peak_rss_children
    page = ModulePage(
        title=f"M5 Differential expression: {spec.name}",
        question=QUESTION,
        concept=CONCEPT,
        parameters=cards,
        results=(_contrasts_table(fits), *groups),
        practice=PRACTICE,
        think=think(spec, fits),
        mistakes=MISTAKES,
        key_numbers=_keys(sheet, fits, spec.de.featured),
        terms=TERMS,
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
            inputs=tuple(f"{p.name} (newest gather output)" for p in inputs.values()),
        ),
    )
    return [write(page, out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"), *files]
