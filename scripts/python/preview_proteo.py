#!/usr/bin/env python3
# /// codebase-meta
# exempt: preview builder; the M7 RNA versus protein page on made-up data, for the course preview
# ///
r"""Course preview of M7 RNA versus protein (``ctx proteo``), on made-up data.

``ctx proteo`` draws the real page; here ``page()`` draws what its page will show from a seeded
generator: protein groups quantified per sample by TMT plex, per-gene RNA-protein Spearman
rho by gene class, RNA against protein fold change in one contrast, and paired z-scores of
the top genes. Samples, cell types and conditions come from the gather sample sheet, so the
preview pages agree on who the samples are. Every sentence is in ``content/proteo_preview.toml``,
apart from ctx proteo's own ``content/proteo.toml``.

The contrast is the sheet's second condition against its first, within its first cell type,
both in sheet order: the preview holds no cohort value.

Usage:
    pixi run python scripts/python/build_preview.py --root <dir> --only proteo
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from time import perf_counter

import numpy as np
import polars as pl
from scipy.special import expit
from scipy.stats import false_discovery_control, kruskal, spearmanr, ttest_ind

from ctx import runtime
from ctx.de.run import Writer
from ctx.de.stats import ALPHA
from ctx.io.naming import MODULES, stem
from ctx.report import text
from ctx.report.builder import (
    FigureBlock,
    FigureText,
    ModulePage,
    ParameterCard,
    Reproduction,
    write,
)
from ctx.viz.bar import BarData
from ctx.viz.curve import CurveData
from ctx.viz.heatmap import HeatmapData, Track, diverging_levels
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import DATA_COLOURS, PROJECTOR_UNSAFE

# The made-up dataset's name on the page, one word for every preview (content/ui.toml).
SHOWN = ("ui", "preview.dataset")

MODULE, PAGE = "proteo", "proteo_preview"
SEED = 7
N_GENES = 3000
# A TMT 11-plex: ten sample channels beside one pooled reference channel.
PER_PLEX = 10
# One colour per plex, the projector-unsafe yellow left out.
PLEX_COLOURS = tuple(c for c in DATA_COLOURS if c not in PROJECTOR_UNSAFE)
# CONVENTIONS.md, Statistics rules: proteins quantified in at least 70 % of samples.
MIN_QUANTIFIED = 0.7
LFC_MIN = 0.5
CI_PCT = 95
BOOTSTRAPS = 1000
TOP = 6
# Share of genes per class, and how far each class's protein follows its RNA: subunits made
# in excess are degraded unless assembled, secreted proteins leave the cell before lysis.
CLASSES = {"complex": (0.12, 0.25), "enzyme": (0.20, 0.85), "secreted": (0.08, 0.30),
           "other": (0.60, 0.65)}  # fmt: skip
# log2 CPM: abundance, biological spread, RNA and protein measurement noise.
ABUNDANCE, SPREAD, RNA_NOISE, PROTEIN_NOISE = (5.0, 2.0), 0.3, 0.25, 0.35
TRACKING_SD = 0.2
CELL_SHARE, CELL_SD = 0.3, 1.0
CHANGED_SHARE, CHANGE = 0.12, (1.5, 0.4)
PROTEIN_ONLY_SHARE, PROTEIN_ONLY = 0.03, (1.2, 0.3)
# Co-isolated peptides compress TMT ratios toward zero.
COMPRESSION = 0.7
PLEX_SD = 0.05
# Identification and per-plex quantification rise with abundance (log2 CPM midpoint, scale).
IDENTIFIED, IN_PLEX, SPORADIC = (3.5, 0.8), (2.5, 0.6), 0.005
# The options ctx proteo will take (TODO.md, M7), with their defaults.
PLANNED = {
    "--method": "spearman",
    "--min-quantified": f"{MIN_QUANTIFIED:g}",
    "--lfc-min": f"{LFC_MIN:g}",
    "--gene-class": "all",
    "--evidence": "validated",
}


@dataclass(frozen=True)
class Design:
    """The samples and the contrast: who is compared with whom, and which plex each ran in."""

    samples: pl.DataFrame
    reference: str
    tested: str
    in_contrast: np.ndarray
    is_tested: np.ndarray
    plex: np.ndarray
    n_plex: int


@dataclass(frozen=True)
class Measured:
    """Made-up RNA and protein per gene and sample; NaN is a protein value not measured."""

    genes: pl.DataFrame
    rna: np.ndarray
    protein: np.ndarray


def _columns(*names: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of each named column, from the content file."""
    table = text.raw(PAGE, "columns")
    return {name: (table[name][0], table[name][1]) for name in names}


def _lines(key: str, **values: object) -> tuple[str, ...]:
    """A figure's reading guide, each line filled from the run's numbers."""
    return tuple(line.format_map(values) for line in text.words(PAGE, key))


def design_sheet(root: Path) -> Path:
    """The gather sheet the design is read from: the newest by name."""
    folder = root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    return sorted(folder.glob("*_samplesheet_v1.parquet"))[-1]


def design(root: Path) -> Design:
    """The newest gather sheet, one row per sample, with a plex and the contrast marked.

    Samples are dealt to plexes in turn after sorting by cell type and condition, so every
    plex holds every group: a group in one plex alone could not be told from the plex.
    """
    sheet = pl.read_parquet(design_sheet(root))
    samples = sheet.unique("sample", keep="first", maintain_order=True).select(
        "sample", "cell_type", "condition"
    )
    reference, tested = samples["condition"].unique(maintain_order=True).to_list()[:2]
    cell = samples["cell_type"][0]
    n_plex = min(math.ceil(samples.height / PER_PLEX), len(PLEX_COLOURS))
    order = np.lexsort((samples["condition"].to_numpy(), samples["cell_type"].to_numpy()))
    plex = np.empty(samples.height, dtype=int)
    plex[order] = np.arange(samples.height) % n_plex
    in_contrast = (samples["cell_type"] == cell).to_numpy() & samples["condition"].is_in(
        [reference, tested]
    ).to_numpy()
    is_tested = (samples["condition"] == tested).to_numpy()
    return Design(samples, reference, tested, in_contrast, is_tested, plex, n_plex)


def simulate(rng: np.random.Generator, plan: Design) -> Measured:
    """RNA and protein for every gene, protein following RNA as far as its class allows."""
    n = plan.samples.height
    classes = rng.choice(list(CLASSES), N_GENES, p=[share for share, _ in CLASSES.values()])
    abundance = rng.normal(*ABUNDANCE, N_GENES)
    other_cell = (plan.samples["cell_type"] != plan.samples["cell_type"][0]).to_numpy()
    cell = np.where(rng.random(N_GENES) < CELL_SHARE, rng.normal(0.0, CELL_SD, N_GENES), 0.0)
    signs = rng.choice([-1.0, 1.0], N_GENES)
    change = np.where(rng.random(N_GENES) < CHANGED_SHARE, signs * rng.normal(*CHANGE, N_GENES), 0)
    spread = rng.normal(0.0, SPREAD, (N_GENES, n))
    biological = (abundance[:, None] + cell[:, None] * other_cell
                  + change[:, None] * plan.is_tested + spread)  # fmt: skip
    rna = biological + rng.normal(0.0, RNA_NOISE, (N_GENES, n))
    tracking = np.clip(rng.normal([CLASSES[c][1] for c in classes], TRACKING_SD), 0.0, 1.2)
    own = np.where((rng.random(N_GENES) < PROTEIN_ONLY_SHARE) & (change == 0),
                   signs * rng.normal(*PROTEIN_ONLY, N_GENES), 0.0)  # fmt: skip
    deviation = biological - biological.mean(axis=1, keepdims=True)
    protein = COMPRESSION * (tracking[:, None] * deviation + own[:, None] * plan.is_tested
                             + rng.normal(0.0, PROTEIN_NOISE, (N_GENES, n)))  # fmt: skip
    protein += rng.normal(0.0, PLEX_SD, (N_GENES, plan.n_plex))[:, plan.plex]
    # Each value is a log2 ratio to the pooled reference channel, so a gene centres on zero.
    protein -= protein.mean(axis=1, keepdims=True)
    identified = rng.random(N_GENES) < expit((abundance - IDENTIFIED[0]) / IDENTIFIED[1])
    in_plex = rng.random((N_GENES, plan.n_plex)) < expit(
        (abundance[:, None] - IN_PLEX[0]) / IN_PLEX[1]
    )
    seen = identified[:, None] & in_plex[:, plan.plex] & (rng.random((N_GENES, n)) >= SPORADIC)
    names = text.raw(PAGE, "classes")
    genes = pl.DataFrame({
        "gene": [f"g{i + 1:04d}" for i in range(N_GENES)],
        "gene_class": [names[c] for c in classes],
        "identified": identified,
    })  # fmt: skip
    return Measured(genes, rna, np.where(seen, protein, np.nan))


def correlate(made: Measured) -> pl.DataFrame:
    """Per-gene Spearman rho of RNA and protein, for genes quantified in enough samples."""
    share = (~np.isnan(made.protein)).mean(axis=1)
    kept = np.flatnonzero(share >= MIN_QUANTIFIED)
    rho, p = [], []
    for g in kept:
        present = ~np.isnan(made.protein[g])
        result = spearmanr(made.rna[g, present], made.protein[g, present])
        rho.append(float(result.statistic))
        p.append(float(result.pvalue))
    padj = false_discovery_control(p, method="bh")
    return made.genes[kept].select("gene", "gene_class").with_columns(
        quantified_share=pl.Series(share[kept]), rho=pl.Series(rho), rho_p=pl.Series(p),
        rho_padj=pl.Series(padj),
    )  # fmt: skip


def fold_changes(made: Measured, plan: Design, tested: pl.DataFrame) -> pl.DataFrame:
    """Welch t-test per gene and layer in the contrast, BH per layer; a call per gene."""
    index = {g: i for i, g in enumerate(made.genes["gene"])}
    rows = [index[g] for g in tested["gene"]]
    case = plan.in_contrast & plan.is_tested
    control = plan.in_contrast & ~plan.is_tested
    out = tested.select("gene", "gene_class")
    for layer, values in (("rna", made.rna[rows]), ("protein", made.protein[rows])):
        test = ttest_ind(values[:, case], values[:, control], axis=1, equal_var=False,
                         nan_policy="omit")  # fmt: skip
        lfc = np.nanmean(values[:, case], axis=1) - np.nanmean(values[:, control], axis=1)
        padj = false_discovery_control(np.asarray(test.pvalue), method="bh")
        out = out.with_columns(**{f"{layer}_lfc": pl.Series(lfc), f"{layer}_padj": pl.Series(padj)})
    changed = {layer: (pl.col(f"{layer}_padj") < ALPHA) & (pl.col(f"{layer}_lfc").abs() >= LFC_MIN)
               for layer in ("rna", "protein")}  # fmt: skip
    word = {k: text.word(PAGE, f"figures.fold.{k}") for k in ("both", "rna", "protein", "neither")}
    call = (pl.when(changed["rna"] & changed["protein"]).then(pl.lit(word["both"]))
            .when(changed["rna"]).then(pl.lit(word["rna"]))
            .when(changed["protein"]).then(pl.lit(word["protein"]))
            .otherwise(pl.lit(word["neither"])))  # fmt: skip
    return out.with_columns(call=call).with_columns(keep=pl.col("call") != word["neither"])


def quantified_figure(writer: Writer, made: Measured, plan: Design) -> FigureBlock:
    """Protein groups with a value in each sample, bars grouped and coloured by plex."""
    quantified = ~np.isnan(made.protein)
    names = [text.fill(PAGE, "figures.quantified.plex", n=k + 1) for k in range(plan.n_plex)]
    frame = (
        plan.samples.with_columns(
            plex=pl.Series([names[k] for k in plan.plex]),
            n_quantified=pl.Series(quantified.sum(axis=0).astype(np.int64)),
            order=pl.Series(plan.plex),
        )
        .sort("order", maintain_order=True)
        .drop("order")
    )
    n_ids = int(made.genes["identified"].sum())
    n_kept = int((quantified.mean(axis=1) >= MIN_QUANTIFIED).sum())
    data = BarData(
        frame=frame, columns=_columns("sample", "cell_type", "condition", "plex", "n_quantified"),
        category="sample", value="n_quantified",
        category_label=text.word(PAGE, "figures.quantified.category"),
        value_label=text.word(PAGE, "figures.quantified.value"), series="plex",
        colours=dict(zip(names, PLEX_COLOURS, strict=False)),
        caption=text.fill(PAGE, "figures.quantified.caption", n_samples=plan.samples.height,
                          n_plex=plan.n_plex, n_ids=n_ids,
                          n_all=int(quantified.all(axis=1).sum())),
        hover=("cell_type", "condition"),
    )  # fmt: skip
    return writer.figure(data, "quantified", FigureText(
        title=text.word(PAGE, "figures.quantified.title"),
        lead=text.fill(PAGE, "figures.quantified.lead", n_kept=n_kept, n_ids=n_ids,
                       min_pct=round(100 * MIN_QUANTIFIED)),
        reading=_lines("figures.quantified.reading")))  # fmt: skip


def rho_figure(writer: Writer, tested: pl.DataFrame, n_samples: int) -> FigureBlock:
    """Per-gene rho as one cumulative curve per gene class; the dashed line marks the medians.

    Curves, not dots: thousands of genes in a few slots would cover one another.
    """
    names = text.raw(PAGE, "classes")
    order = [names[c] for c in CLASSES]
    frame = (
        tested.select("gene", "gene_class", "rho")
        .with_columns(slot=pl.col("gene_class").replace_strict(order, list(range(len(order)))))
        .sort("slot", "rho")
        .with_columns(share=(pl.col("rho").cum_count().over("slot") / pl.len().over("slot")))
        .drop("slot")
    )
    by_class = [frame.filter(pl.col("gene_class") == c)["rho"].to_numpy() for c in order]
    h, p = kruskal(*by_class)
    counts = ", ".join(text.fill(PAGE, "classes.count", name=c, n=v.size)
                       for c, v in zip(order, by_class, strict=True))  # fmt: skip
    values = {"n_samples": n_samples, "n_genes": frame.height, "alpha": ALPHA}
    ends = (float(frame["rho"].min()), float(frame["rho"].max()))  # type: ignore[arg-type]
    data = CurveData(
        frame=frame, columns=_columns(*frame.columns), series="gene_class", x="rho", y="share",
        x_label=text.word(PAGE, "figures.rho.x"), y_label=text.word(PAGE, "figures.rho.y"),
        colour="gene_class", colours={}, distinct=True, hover=("gene",),
        caption=text.fill(PAGE, "figures.rho.caption", **values, counts=counts, h=h,
                          p=f"{p:.2g}", eps=h / (frame.height - 1)),
        thresholds=(Line(text.word(PAGE, "figures.rho.median"), ends, (0.5, 0.5), dash=True),),
    )  # fmt: skip
    n_sig = tested.filter((pl.col("rho_padj") < ALPHA) & (pl.col("rho") > 0)).height
    return writer.figure(data, "rho_by_class", FigureText(
        title=text.word(PAGE, "figures.rho.title"),
        lead=text.fill(PAGE, "figures.rho.lead", **values, median=frame["rho"].median(),
                       n_sig=n_sig),
        reading=_lines("figures.rho.reading", **values)))  # fmt: skip


def concordance(rng: np.random.Generator, folds: pl.DataFrame) -> tuple[float, float, float]:
    """Spearman rho of RNA against protein fold change across genes, and its bootstrap CI."""
    x, y = folds["rna_lfc"].to_numpy(), folds["protein_lfc"].to_numpy()
    rho = float(spearmanr(x, y).statistic)
    draws = rng.integers(0, x.size, (BOOTSTRAPS, x.size))
    boot = [float(spearmanr(x[i], y[i]).statistic) for i in draws]
    tail = (100 - CI_PCT) / 2
    return rho, float(np.percentile(boot, tail)), float(np.percentile(boot, 100 - tail))


def fold_figure(writer: Writer, folds: pl.DataFrame, plan: Design,
                ci: tuple[float, float, float]) -> FigureBlock:  # fmt: skip
    """RNA against protein fold change per gene, coloured by the layers it changed at."""
    rho, low, high = ci
    x, y = folds["rna_lfc"], folds["protein_lfc"]
    signs = (("up", 1), ("down", -1))
    quadrant = {f"{a}_{b}": int(((sx * x > 0) & (sy * y > 0)).sum())
                for a, sx in signs for b, sy in signs}  # fmt: skip
    top = float(max(x.abs().max(), y.abs().max())) * 1.05  # type: ignore[arg-type]
    zero = text.word(PAGE, "figures.fold.zero")
    lines = (Line(text.word(PAGE, "figures.fold.identity"), (-top, top), (-top, top)),
             Line(zero, (-top, top), (0.0, 0.0), dash=True),
             Line(zero, (0.0, 0.0), (-top, top), dash=True))  # fmt: skip
    counts = {"n_ref": int((plan.in_contrast & ~plan.is_tested).sum()),
              "n_tested": int((plan.in_contrast & plan.is_tested).sum()),
              "reference": plan.reference, "tested": plan.tested}  # fmt: skip
    word = {k: text.word(PAGE, f"figures.fold.{k}") for k in ("both", "rna", "protein", "neither")}
    data = ScatterData(
        frame=folds, columns=_columns(*folds.columns), x="rna_lfc", y="protein_lfc",
        x_label=text.word(PAGE, "figures.fold.x"), y_label=text.word(PAGE, "figures.fold.y"),
        colour="call", colours={word["both"]: DATA_COLOURS[2], word["rna"]: DATA_COLOURS[0],
                                word["protein"]: DATA_COLOURS[4], word["neither"]: None},
        caption=text.fill(PAGE, "figures.fold.caption", **counts, **quadrant, rho=rho,
                          ci=CI_PCT, lo=low, hi=high),
        hover=("gene", "gene_class", "rna_padj", "protein_padj"), lines=lines, keep="keep",
    )  # fmt: skip
    return writer.figure(data, "fold_change", FigureText(
        title=text.fill(PAGE, "figures.fold.title", **counts),
        lead=text.fill(PAGE, "figures.fold.lead", rho=rho, n_genes=folds.height,
                       n_same=quadrant["up_up"] + quadrant["down_down"]),
        reading=_lines("figures.fold.reading", alpha=ALPHA, lfc=LFC_MIN)))  # fmt: skip


def _z(values: np.ndarray) -> np.ndarray:
    """Each row centred and scaled over its present values."""
    centred = values - np.nanmean(values, axis=1, keepdims=True)
    return centred / np.nanstd(values, axis=1, keepdims=True)


def paired_figure(writer: Writer, made: Measured, plan: Design,
                  tables: tuple[pl.DataFrame, pl.DataFrame]) -> FigureBlock:  # fmt: skip
    """z-scores of the top genes, each gene's RNA row above its protein row, samples by group."""
    tested, folds = tables
    ranked = folds.join(tested.select("gene", "rho"), on="gene").sort("rna_padj")
    picks = [ranked.filter(pl.col("call") == text.word(PAGE, f"figures.fold.{k}")).head(TOP)
             for k in ("both", "rna")]  # fmt: skip
    genes = pl.concat(picks)["gene"].to_list()
    rows = [made.genes["gene"].to_list().index(g) for g in genes]
    columns = np.flatnonzero(plan.in_contrast)
    columns = columns[np.argsort(plan.is_tested[columns], kind="stable")]
    names = plan.samples["sample"].to_list()
    shown = [names[j] for j in columns]
    labels, values = [], []
    for g, row in zip(genes, rows, strict=True):
        for key, layer in (("rna_row", made.rna), ("protein_row", made.protein)):
            labels.append(text.fill(PAGE, f"figures.paired.{key}", gene=g))
            values.append(layer[row, columns])
    z = _z(np.array(values))
    flat = [None if math.isnan(v) else float(v) for v in z.ravel()]
    levels, swatches = diverging_levels(flat)
    frame = pl.DataFrame({"row": [r for r in labels for _ in shown], "sample": shown * len(labels),
                          "z": flat, "level": levels}, strict=False)  # fmt: skip
    rho = dict(ranked.select("gene", "rho").iter_rows())
    groups = dict(zip(names, plan.samples["condition"], strict=True))
    counts = {"n_ref": int((~plan.is_tested[columns]).sum()),
              "n_tested": int(plan.is_tested[columns].sum()), "reference": plan.reference,
              "tested": plan.tested, "n_samples": plan.samples.height}  # fmt: skip
    data = HeatmapData(
        frame=frame, columns=_columns("row", "sample", "z", "level"), row="row",
        column="sample", value="z", level="level",
        row_label=text.word(PAGE, "figures.paired.row"),
        column_label=text.word(PAGE, "figures.paired.column"), swatches=swatches,
        caption=text.fill(PAGE, "figures.paired.caption", **counts, n_shown=len(shown)),
        row_order=tuple(labels), column_order=tuple(shown),
        tracks=(Track(text.word(PAGE, "figures.paired.track"),
                      {s: groups[s] for s in shown},
                      {plan.reference: DATA_COLOURS[2], plan.tested: DATA_COLOURS[6]}),),
        row_stats={label: rho[genes[i // 2]] for i, label in enumerate(labels)},
        row_stats_label=text.word(PAGE, "figures.paired.rho"), row_stats_format=".2f",
    )  # fmt: skip
    return writer.figure(data, "paired_zscores", FigureText(
        title=text.word(PAGE, "figures.paired.title"),
        lead=text.fill(PAGE, "figures.paired.lead", n_both=picks[0].height,
                       n_rna=picks[1].height),
        reading=_lines("figures.paired.reading", **counts)))  # fmt: skip


def page(root: Path, dataset: str, day: date) -> Path:
    """Write the M7 preview page and its figures under ``outputs/07_proteo``.

    Args:
        root: The ctx root holding the gather outputs of the made-up dataset.
        dataset: The dataset's short name, in every output name.
        day: The date in every output name.

    Returns:
        The page's path.
    """
    started = perf_counter()
    rng = np.random.default_rng(SEED)
    plan = design(root)
    made = simulate(rng, plan)
    tested = correlate(made)
    folds = fold_changes(made, plan, tested)
    ci = concordance(rng, folds)
    writer = Writer(root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset, day, {})
    merged = tested.join(folds.drop("gene_class", "keep"), on="gene")
    writer.table(merged, "genes", _columns(*merged.columns))
    figures = (quantified_figure(writer, made, plan),
               rho_figure(writer, tested, plan.samples.height),
               fold_figure(writer, folds, plan, ci),
               paired_figure(writer, made, plan, (tested, folds)))  # fmt: skip
    cards = tuple(ParameterCard(flag, text.word(PAGE, f"planned.{flag}.meaning"), default,
                                text.word(PAGE, f"planned.{flag}.effect"))
                  for flag, default in PLANNED.items())  # fmt: skip
    usage, budget = runtime.measure(started), runtime.budget(MODULE)
    track = "A" if MODULE in runtime.TRACK_A_MODULES else "B"
    built = ModulePage(
        title=text.fill(PAGE, "title", dataset=text.word(*SHOWN)), **text.skeleton(PAGE),
        parameters=cards, results=(text.word("ui", "preview.note"), *figures),
        reproduce=Reproduction(
            command=" ".join(["ctx", MODULE, "--dataset", dataset,
                              *(part for pair in PLANNED.items() for part in pair)]),
            wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children, profile=text.word(PAGE, "profile"),
            budget=text.fill(PAGE, "budget", track=track, minutes=budget.minutes,
                             ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=(design_sheet(root).name,),
        ),
    )  # fmt: skip
    return write(built, writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html")
