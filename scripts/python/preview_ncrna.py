#!/usr/bin/env python3
# /// codebase-meta
# exempt: preview builder; draws the M8 Non-coding RNA page on made-up data until ctx ncrna is built
# ///
"""The M8 Non-coding RNA preview page (``ctx ncrna``), drawn on made-up counts.

Counts come from a seeded negative binomial for the samples of the made-up dataset's sheet, so
the page names the same samples as the other preview pages. The figures use the existing
contracts of ``ctx.viz`` and the expression filter of M5 (``ctx.de.filters``); every sentence
is in ``content/ncrna.toml``. Called by scripts/python/build_preview.py.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control, ks_2samp, mannwhitneyu, spearmanr, ttest_ind

from ctx import runtime
from ctx.de.filters import MIN_COUNT, TOLERANCE, expressed, min_sample_size
from ctx.de.run import Writer, newest
from ctx.de.stats import ALPHA, LFC_MIN, classify
from ctx.io.naming import MODULES, stem
from ctx.report import text
from ctx.report.builder import (
    FigureBlock,
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    write,
)
from ctx.viz.bar import BarData
from ctx.viz.curve import CurveData
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import DATA_COLOURS, PROJECTOR_UNSAFE, ROLES

PAGE = "ncrna"
SEED = 8
LNCRNA, CODING, MT_RRNA = "lncRNA", "protein_coding", "Mt_rRNA"
# Made-up genes per biotype: (genes, mean and SD of the log2 expression weight). Tuned so the
# read shares sit near a poly-A library's, about nine in ten on protein-coding genes, and the
# filter keeps most protein-coding genes and a minority of lncRNAs.
GENES = {"protein_coding": (2000, 2.5, 3.0), "lncRNA": (1500, -1.0, 2.5),
         "pseudogene": (800, -1.0, 2.5), "Mt_rRNA": (2, 10.0, 0.5), "small RNA": (300, 0.0, 2.5),
         "other": (200, -1.0, 2.5)}  # fmt: skip
# Yellow is left out: at 1.32:1 on white it does not survive a projector.
COLOURS = dict(zip(GENES, [c for c in DATA_COLOURS if c not in PROJECTOR_UNSAFE], strict=True))
# log2 SDs: per-sample noise, cell-type offset, mitochondrial rRNA per sample, and the spread of
# every gene's true change between conditions; NB dispersion.
NOISE_SD, CELL_SD, MT_SD, CHANGE_SD, DISPERSION = 0.3, 0.4, 0.8, 0.25, 0.05
# Reads per library, scaled to the made-up gene count rather than to a real transcriptome.
DEPTH, DEPTH_SD = 5e5, 0.35
# Planted in the tested condition: a share of each tested biotype up or down, and target shifts.
DE_SHARE, DE_LFC = 0.06, (1.0, 2.5)
TARGETS = {"validated": (40, -0.3), "predicted": (200, -0.1)}
# The lncRNA and its neighbour: log2 weights, a shared condition and cell-type effect, and a
# shared per-sample term, so most of their correlation is the design, not the pair.
PAIR, PAIR_MU, PAIR_LFC, PAIR_CELL, PAIR_SD = (LNCRNA, CODING), (5.0, 6.0), 1.0, 0.8, 0.15


# One meaning per column, for every figure's CSV README.
COLUMNS = {
    "gene_id": ("none", "Made-up gene."),
    "biotype": ("none", "GENCODE biotype; small RNA and pseudogene pool several types."),
    "sample": ("none", "Sample, from the gather sample sheet."),
    "cell_type": ("none", "Cell type, from the sample sheet."),
    "condition": ("none", "Condition, from the sample sheet."),
    "percent": ("%", "Share of the sample's counted reads."),
    "value": (
        "log2 reads",
        "log2 of 1 plus the reads, at the median library size, that the "
        "filter's sample count reaches.",
    ),
    "share": ("fraction", "Share of the line's genes at or below this value."),
    "target": ("none", "Made-up target class of the made-up miRNA: validated, predicted, non."),
    "target_class": ("none", "Made-up target class, as the legend names it."),
    "pair": ("none", "True for the lncRNA and neighbour of the co-expression figure."),
    "lfc": ("log2 fold change", "Mean log2 CPM, tested minus reference condition."),
    "p": ("dimensionless", "Welch two-sample t-test p-value."),
    "padj": ("dimensionless", "Benjamini-Hochberg adjusted p within the biotype."),
    "neg_log10_p": ("-log10 p", "Minus log10 of p."),
    "call": ("none", "up, down or ns: padj below alpha and |lfc| at least the threshold."),
    "level": ("none", "The biotype of a called gene; ns otherwise."),
    "lncrna": ("log2 CPM", "The made-up lncRNA, log2 of CPM plus 1."),
    "neighbour": ("log2 CPM", "Its made-up neighbour mRNA, log2 of CPM plus 1."),
    "lncrna_residual": ("log2 CPM", "lncrna minus its cell type and condition mean."),
    "neighbour_residual": ("log2 CPM", "neighbour minus its cell type and condition mean."),
}


def _columns(frame: pl.DataFrame) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of a figure's table."""
    return {c: COLUMNS[c] for c in frame.columns}


def _log_cpm(counts: Any) -> Any:
    """log2 of 1 plus counts per million, genes by samples."""
    return np.log2(counts / counts.sum(axis=0) * 1e6 + 1)


def tested_biotypes() -> tuple[str, ...]:
    """The biotypes the page tests: the planned ``--biotypes`` default, one definition."""
    return tuple(text.word(PAGE, "planned.--biotypes.default").split(","))


def _text(key: str, **values: Any) -> FigureText:
    """Title, lead and reading guide of ``figures.<key>``, filled from the same values."""
    reading = tuple(line.format_map(values) for line in text.raw(PAGE, f"figures.{key}.reading"))
    title, lead = (text.fill(PAGE, f"figures.{key}.{part}", **values) for part in ("title", "lead"))
    return FigureText(title=title, lead=lead, reading=reading)


def _caption(key: str, **values: Any) -> str:
    """The caption of ``figures.<key>``, filled."""
    return text.fill(PAGE, f"figures.{key}.caption", **values)


def made_up(sheet: pl.DataFrame, tested: str, rng: np.random.Generator) -> tuple[pl.DataFrame, Any]:
    """Made-up genes, their planted effects and their counts in the sheet's samples.

    Returns:
        One row per gene (``gene_id``, ``biotype``, ``target``, ``pair``) and the counts,
        genes by samples in sheet order.
    """
    biotype = np.concatenate([[b] * n for b, (n, _, _) in GENES.items()])
    mu = np.concatenate([rng.normal(m, sd, n) for n, m, sd in GENES.values()])
    effect, cell = rng.normal(0.0, CHANGE_SD, mu.size), rng.normal(0.0, CELL_SD, mu.size)
    target = np.full(mu.size, "non", dtype=object)
    pair = [int(np.flatnonzero(biotype == b)[0]) for b in PAIR]
    mu[pair], effect[pair], cell[pair] = PAIR_MU, PAIR_LFC, PAIR_CELL
    free = np.setdiff1d(np.flatnonzero(np.isin(biotype, tested_biotypes())), pair)
    de = rng.choice(free, int(DE_SHARE * free.size), replace=False)
    effect[de] = rng.choice([-1.0, 1.0], de.size) * rng.uniform(*DE_LFC, de.size)
    coding = np.setdiff1d(np.flatnonzero(biotype == CODING), [*pair, *de])
    picked = rng.choice(coding, sum(n for n, _ in TARGETS.values()), replace=False)
    split = np.split(picked, [TARGETS["validated"][0]])
    for (name, (_, shift)), chosen in zip(TARGETS.items(), split, strict=True):
        target[chosen] = name
        effect[chosen] += shift
    is_tested = (sheet["condition"] == tested).to_numpy()
    other_cell = (sheet["cell_type"] != sheet["cell_type"][0]).to_numpy()
    log2 = (mu[:, None] + effect[:, None] * is_tested + cell[:, None] * other_cell
            + rng.normal(0.0, NOISE_SD, (mu.size, sheet.height)))  # fmt: skip
    log2[pair] += rng.normal(0.0, PAIR_SD, sheet.height)
    log2[biotype == MT_RRNA] += rng.normal(0.0, MT_SD, sheet.height)
    weight = np.exp2(log2)
    mean = weight / weight.sum(axis=0) * DEPTH * np.exp(rng.normal(0.0, DEPTH_SD, sheet.height))
    counts = rng.poisson(rng.gamma(1 / DISPERSION, mean * DISPERSION))
    genes = pl.DataFrame({"gene_id": [f"MADEUP{i + 1:05d}" for i in range(mu.size)],
                          "biotype": biotype.astype(str), "target": target.astype(str),
                          "pair": np.isin(np.arange(mu.size), pair)})  # fmt: skip
    return genes, counts


def composition(
    writer: Writer, sheet: pl.DataFrame, genes: pl.DataFrame, counts: Any
) -> FigureBlock:
    """Stacked bars: each sample's counted reads by biotype."""
    shares = np.stack([counts[(genes["biotype"] == b).to_numpy()].sum(axis=0) for b in GENES])
    shares = 100 * shares / counts.sum(axis=0)
    frame = pl.DataFrame(
        [(s, b, float(shares[i, j]), c, k) for j, (s, c, k) in
         enumerate(sheet.select("sample", "cell_type", "condition").rows())
         for i, b in enumerate(GENES)],
        schema=["sample", "biotype", "percent", "cell_type", "condition"], orient="row",
    )  # fmt: skip
    mt, lncrna = (shares[list(GENES).index(b)] for b in (MT_RRNA, LNCRNA))
    data = BarData(
        frame=frame, columns=_columns(frame), category="sample", value="percent",
        category_label=text.word(PAGE, "figures.composition.x"),
        value_label=text.word(PAGE, "figures.composition.y"), series="biotype", colours=COLOURS,
        caption=_caption("composition", n_samples=sheet.height, n_genes=genes.height),
        stacked=True, hover=("cell_type", "condition"), width="double",
    )  # fmt: skip
    words = _text("composition", lncrna=float(np.median(lncrna)), mt_low=float(mt.min()),
                  mt_high=float(mt.max()))  # fmt: skip
    return writer.figure(data, "biotype_shares", words)


@dataclass(frozen=True)
class Contrast:
    """The made-up contrast: the sheet's first cell type, its first condition as reference."""

    cell_type: str
    reference: str
    tested: str
    columns: Any  # bool over the sheet's samples: the contrast's own
    is_tested: Any  # bool over the contrast's samples

    @classmethod
    def of(cls, sheet: pl.DataFrame) -> Contrast:
        """The contrast of a sample sheet, one row per sample."""
        cell = str(sheet["cell_type"][0])
        within = sheet.filter(pl.col("cell_type") == cell)
        reference, tested = within["condition"].unique(maintain_order=True).to_list()[:2]
        columns = (sheet["cell_type"] == cell) & sheet["condition"].is_in([reference, tested])
        is_tested = sheet.filter(columns)["condition"] == tested
        return cls(cell, reference, tested, columns.to_numpy(), is_tested.to_numpy())

    def words(self) -> dict[str, Any]:
        """The placeholders every caption on the contrast fills."""
        n_test = int(self.is_tested.sum())
        return {"cell_type": self.cell_type, "reference": self.reference, "tested": self.tested,
                "n_ref": self.is_tested.size - n_test, "n_test": n_test}  # fmt: skip


def expression_filter(writer: Writer, genes: pl.DataFrame, counts: Any,
                      con: Contrast) -> tuple[FigureBlock, Any]:  # fmt: skip
    """Cumulative curves of the value M5's filter tests, per biotype; and which genes pass."""
    sub = counts[:, con.columns]
    sizes = [int((~con.is_tested).sum()), int(con.is_tested.sum())]
    keep = expressed(sub, sizes, MIN_COUNT)
    k = math.ceil(min_sample_size(sizes) - TOLERANCE)
    library = sub.sum(axis=0)
    # The filter's CPM cutoff is MIN_COUNT reads at the median library size, so on that scale
    # the cutoff is MIN_COUNT itself, whatever the depth.
    at_median = np.sort(sub / library * float(np.median(library)))[:, -k]
    reached = genes.with_columns(value=pl.Series(np.log2(at_median + 1)))
    cutoff = math.log2(MIN_COUNT + 1)
    frame = (reached.filter(pl.col("biotype").is_in(tested_biotypes())).sort("biotype", "value")
             .with_columns(share=((pl.int_range(pl.len()) + 1) / pl.len()).over("biotype"))
             .select("gene_id", "biotype", "value", "share"))  # fmt: skip
    value = {b: frame.filter(pl.col("biotype") == b)["value"].to_numpy() for b in PAIR}
    u, p = mannwhitneyu(value[LNCRNA], value[CODING])
    kept = {b: 100 * float(keep[(genes["biotype"] == b).to_numpy()].mean()) for b in PAIR}
    data = CurveData(
        frame=frame, columns=_columns(frame), series="biotype", x="value", y="share",
        x_label=text.word(PAGE, "figures.filter.x"), y_label=text.word(PAGE, "figures.filter.y"),
        colour="biotype", colours={b: COLOURS[b] for b in tested_biotypes()},
        caption=_caption("filter", **con.words(), min_count=MIN_COUNT, p=p,
                         r=2 * u / (value[LNCRNA].size * value[CODING].size) - 1),
        hover=("gene_id",),
        thresholds=(Line(text.word(PAGE, "figures.filter.cutoff"), (cutoff, cutoff), (0.0, 1.0),
                         dash=True),),
    )  # fmt: skip
    words = _text("filter", protein_coding=kept[CODING], lncrna=kept[LNCRNA], k=k)
    return writer.figure(data, "filter_by_biotype", words), keep


def volcano(writer: Writer, genes: pl.DataFrame, logs: Any, keep: Any,
            con: Contrast) -> tuple[FigureBlock, pl.DataFrame]:  # fmt: skip
    """Welch t-test per gene on log2 CPM, BH within each biotype; dots coloured by biotype."""
    mask = keep & genes["biotype"].is_in(tested_biotypes()).to_numpy()
    logs = logs[mask][:, con.columns]
    case, ref = logs[:, con.is_tested], logs[:, ~con.is_tested]
    frame = genes.filter(pl.Series(mask)).with_columns(
        lfc=pl.Series(case.mean(axis=1) - ref.mean(axis=1)),
        p=pl.Series(ttest_ind(case, ref, axis=1, equal_var=False).pvalue),
    )
    padj = np.empty(frame.height)
    for b in tested_biotypes():
        within = (frame["biotype"] == b).to_numpy()
        padj[within] = false_discovery_control(frame["p"].to_numpy()[within], method="bh")
    frame = (
        frame.with_columns(padj=pl.Series(padj), neg_log10_p=-pl.col("p").log10())
        .with_columns(call=classify(pl.col("padj"), pl.col("lfc")))
        .with_columns(level=pl.when(pl.col("call") == "ns").then(pl.lit("ns")).otherwise("biotype"))
    )
    n = {b: frame.filter(pl.col("biotype") == b).height for b in PAIR}
    called = {b: frame.filter(pl.col("level") == b).height for b in PAIR}
    top = float(frame["neg_log10_p"].max()) * 1.05
    name = text.word(PAGE, "figures.volcano.threshold")
    data = ScatterData(
        frame=frame, columns=_columns(frame), x="lfc", y="neg_log10_p",
        x_label=text.fill(PAGE, "figures.volcano.x", **con.words()),
        y_label=text.word(PAGE, "figures.volcano.y"),
        colour="level", colours={**{b: COLOURS[b] for b in tested_biotypes()}, "ns": None},
        names={"ns": text.word(PAGE, "figures.volcano.ns")},
        caption=_caption("volcano", **con.words(), n_tested=frame.height,
                         n_protein_coding=n[CODING], n_lncrna=n[LNCRNA]),
        hover=("gene_id", "padj"),
        lines=tuple(Line(name, (x, x), (0.0, top), dash=True) for x in (-LFC_MIN, LFC_MIN)),
        width="double",
    )  # fmt: skip
    words = _text("volcano", **con.words(), called_protein_coding=called[CODING],
                  called_lncrna=called[LNCRNA], alpha=ALPHA, lfc_min=LFC_MIN)  # fmt: skip
    return writer.figure(data, "volcano_by_biotype", words), frame


def neighbour(writer: Writer, sheet: pl.DataFrame, genes: pl.DataFrame, logs: Any,
              con: Contrast) -> FigureBlock:  # fmt: skip
    """The lncRNA against its neighbour, per sample; Spearman on values and on residuals."""
    rows = genes.with_row_index().filter("pair")
    row = dict(zip(rows["biotype"], rows["index"], strict=True))
    frame = (
        sheet.select("sample", "cell_type", "condition")
        .with_columns(lncrna=logs[row[LNCRNA]], neighbour=logs[row[CODING]])
        .filter(pl.col("condition").is_in([con.reference, con.tested]))
    )
    frame = frame.with_columns(
        (pl.col(c) - pl.col(c).mean().over("cell_type", "condition")).alias(f"{c}_residual")
        for c in ("lncrna", "neighbour")
    )
    rho, p = spearmanr(frame["lncrna"], frame["neighbour"])
    rho_within, p_within = spearmanr(frame["lncrna_residual"], frame["neighbour_residual"])
    names = {"lncrna": text.word(PAGE, "preview.lncrna"),
             "neighbour": text.word(PAGE, "preview.neighbour")}  # fmt: skip
    data = ScatterData(
        frame=frame, columns=_columns(frame), x="lncrna", y="neighbour",
        x_label=text.fill(PAGE, "figures.neighbour.x", **names),
        y_label=text.fill(PAGE, "figures.neighbour.y", **names), colour="condition",
        colours={con.reference: ROLES["reference"], con.tested: ROLES["comparison"]},
        symbol="cell_type", hover=("sample",),
        caption=_caption("neighbour", n_samples=frame.height, rho=rho, p=p,
                         rho_within=rho_within, p_within=p_within),
    )  # fmt: skip
    words = _text("neighbour", **names, rho=rho, rho_within=rho_within)
    return writer.figure(data, "lncrna_neighbour", words)


def target_shift(writer: Writer, results: pl.DataFrame, con: Contrast) -> FigureBlock:
    """Cumulative fold changes of tested protein-coding genes by target class; KS, BH."""
    labels = text.words(PAGE, "figures.targets.classes")
    classes = dict(zip((*TARGETS, "non"), labels, strict=True))
    frame = (results.filter(pl.col("biotype") == CODING)
             .with_columns(target_class=pl.col("target").replace_strict(classes))
             .sort("target", "lfc")
             .with_columns(share=((pl.int_range(pl.len()) + 1) / pl.len()).over("target"))
             .select("gene_id", "target_class", "lfc", "share"))  # fmt: skip
    lfc = {k: frame.filter(pl.col("target_class") == v)["lfc"].to_numpy()
           for k, v in classes.items()}  # fmt: skip
    tests = {k: ks_2samp(lfc[k], lfc["non"]) for k in TARGETS}
    adjusted = false_discovery_control([t.pvalue for t in tests.values()])
    padj = dict(zip(TARGETS, adjusted, strict=True))
    shift = {k: float(np.median(lfc[k]) - np.median(lfc["non"])) for k in TARGETS}
    data = CurveData(
        frame=frame, columns=_columns(frame), series="target_class", x="lfc", y="share",
        x_label=text.fill(PAGE, "figures.targets.x", **con.words()),
        y_label=text.word(PAGE, "figures.targets.y"),
        colour="target_class", hover=("gene_id",),
        colours={classes["validated"]: ROLES["comparison"], classes["predicted"]: DATA_COLOURS[0],
                 classes["non"]: ROLES["reference"]},
        caption=_caption("targets", n_validated=lfc["validated"].size,
                         n_predicted=lfc["predicted"].size, n_non=lfc["non"].size,
                         d_validated=tests["validated"].statistic, p_validated=padj["validated"],
                         d_predicted=tests["predicted"].statistic, p_predicted=padj["predicted"]),
    )  # fmt: skip
    words = _text("targets", mirna=text.word(PAGE, "preview.mirna"),
                  **{f"shift_{k}": v for k, v in shift.items()})  # fmt: skip
    return writer.figure(data, "target_shift", words)


def cards() -> tuple[ParameterCard, ...]:
    """The options ``ctx ncrna`` will have, their text from ``planned`` in the content file."""
    planned = text.raw(PAGE, "planned")
    return tuple(
        ParameterCard(flag, text.word(PAGE, f"planned.{flag}.meaning"),
                      text.fill(PAGE, f"planned.{flag}.default", min_count=MIN_COUNT),
                      text.word(PAGE, f"planned.{flag}.effect"))
        for flag in planned
    )  # fmt: skip


def key_numbers(genes: pl.DataFrame, results: pl.DataFrame) -> tuple[KeyNumber, ...]:
    """Genes tested, lncRNAs tested and genes called."""
    in_tested = genes.filter(pl.col("biotype").is_in(tested_biotypes())).height
    lncrna = results.filter(pl.col("biotype") == LNCRNA).height
    all_lncrna = genes.filter(pl.col("biotype") == LNCRNA).height
    called = results.filter(pl.col("call") != "ns").height
    return (
        KeyNumber(text.word(PAGE, "keys.tested"), f"{results.height:,}",
                  text.fill(PAGE, "keys.tested_detail", n_genes=in_tested)),
        KeyNumber(text.word(PAGE, "keys.lncrna"), f"{lncrna:,}",
                  text.fill(PAGE, "keys.lncrna_detail", pct=100 * lncrna / all_lncrna)),
        KeyNumber(text.word(PAGE, "keys.called"), f"{called:,}",
                  text.fill(PAGE, "keys.called_detail", alpha=ALPHA, lfc_min=LFC_MIN)),
    )  # fmt: skip


def page(root: Path, dataset: str, day: date) -> Path:
    """Write the M8 preview page and its figures under ``root/outputs/08_ncrna/``.

    Args:
        root: A ctx root holding the made-up dataset's gather outputs.
        dataset: The dataset slug, which names every output.
        day: The date every output name carries.

    Returns:
        The page.
    """
    started = perf_counter()
    gathered = root / "outputs" / "01_gather"
    sheet_path = newest(gathered, dataset, "samplesheet")
    sheet = pl.read_parquet(sheet_path).unique("sample", keep="first", maintain_order=True)
    con = Contrast.of(sheet)
    genes, counts = made_up(sheet, con.tested, np.random.default_rng(SEED))
    logs = _log_cpm(counts)
    out = root / "outputs" / f"{dict(MODULES)[PAGE]:02d}_{PAGE}"
    writer = Writer(out, dataset, day, {"input:samplesheet": sheet_path.name,
                                        "made-up data": f"numpy default_rng({SEED})"})  # fmt: skip
    shares = composition(writer, sheet, genes, counts)
    curves, keep = expression_filter(writer, genes, counts, con)
    dots, results = volcano(writer, genes, logs, keep, con)
    groups = (
        Group(text.word(PAGE, "groups.library"), (shares, curves)),
        Group(text.word(PAGE, "groups.de"), (dots,)),
        Group(text.word(PAGE, "groups.lncrna"), (neighbour(writer, sheet, genes, logs, con),)),
        Group(text.word(PAGE, "groups.targets"), (target_shift(writer, results, con),)),
    )
    usage, budget = runtime.measure(started), runtime.budget(PAGE)
    profile = runtime.resolve_profile(None, os.environ)
    report = ModulePage(
        title=text.fill(PAGE, "title", dataset=dataset), **text.skeleton(PAGE),
        parameters=cards(), results=(text.word("ui", "preview.note"), *groups),
        key_numbers=key_numbers(genes, results),
        reproduce=Reproduction(
            command=text.fill(PAGE, "preview.command", dataset=dataset),
            wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=text.fill(PAGE, "preview.profile", profile=profile),
            budget=text.fill(PAGE, "preview.budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=(sheet_path.name,),
        ),
    )  # fmt: skip
    return write(report, out / f"{stem(PAGE, f'{dataset}_overview', 1, day=day)}.html")
