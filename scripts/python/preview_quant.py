#!/usr/bin/env python3
# /// codebase-meta
# exempt: preview builder; draws the M4 Quantification page on made-up data until ctx quant exists
# ///
"""The M4 Quantification preview: the page ``ctx quant`` will write, drawn from made-up data.

Every number is drawn from one seeded generator over the made-up dataset's own samples, read
from its gather sheet, so the preview pages name the same samples. The text is
``content/quant_preview.toml``; the sample correlation's bins are ``ctx quant``'s own.
Called by ``build_preview.py``.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu, spearmanr

from ctx import runtime
from ctx.de.run import Writer
from ctx.io.naming import MODULES
from ctx.quant.figures import rho_levels
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
from ctx.viz.heatmap import HeatmapData, Track
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.strip import StripData
from ctx.viz.style import DATA_COLOURS, category_colours

# The made-up dataset's name on the page, one word for every preview (content/ui.toml).
SHOWN = ("ui", "preview.dataset")

MODULE, PAGE = "quant", "quant_preview"
SEED = 4
# A gene is detected at this many reads; the figures only, the real module has no such option.
MIN_READS = 10
PERMUTATIONS = 9_999
# Made-up but plausible for bulk RNA-seq: Salmon maps 80 to 90 % of reads, and a library of
# about 20 million reads detects 15,000 to 20,000 genes at 10 reads.
RATE = ((83.5, 2.0), (86.5, 2.0))  # mapping rate %, mean and SD, per batch in sheet order
DEPTH = (24.0, 0.35)  # million reads sequenced: median and SD of its log
CEILING = ((17_800.0, 500.0), (20_300.0, 500.0))  # genes at full depth, per cell type
HALF = (3.5, 0.3)  # million reads at which half the ceiling is detected
STEPS = 15  # subsampled depths per detection curve
GENES = 3_000  # genes behind the correlation matrix
# (sheet column, share of genes it shifts, SD of the log2 shift) for samples not in its
# first level: cell type dominates, as it does in immune-cell data.
EFFECTS = (("cell_type", 0.25, 1.6), ("batch", 0.08, 0.7), ("condition", 0.04, 1.0))
UNIVERSE, DRAWN = 20_000, 2_500  # genes in one sample, and those the units figure draws
LENGTH = (1_800.0, 0.7, 150.0, 12_000.0)  # effective length nt: median, log SD, floor, ceiling
SHORT_NT, LONG_NT = 1_000, 4_000
# Genes on the tool comparison; a tenth are paralogs whose reads map to several loci, and a
# twentieth overlap another gene, which alignment counting calls ambiguous.
TOOL_GENES, PARALOGS, OVERLAPPING = 30_000, 0.10, 0.05
UNITS = {
    "sample": "none", "cell_type": "none", "condition": "none", "batch": "none",
    "mapping_rate": "%", "assigned": "million reads", "genes_at_depth": "genes",
    "depth": "million reads", "genes": "genes", "sample_row": "none", "sample_column": "none",
    "rho": "dimensionless", "level": "none", "gene": "none", "tpm": "TPM", "reads": "reads",
    "effective_length": "nt", "length_class": "none", "region": "none", "salmon": "none",
    "n_genes": "genes", "median_reads": "reads",
}  # fmt: skip


@dataclass(frozen=True)
class Design:
    """The made-up dataset's samples, and the two levels each comparison sets side by side.

    Attributes:
        sheet: One row per sample: sample, cell_type, condition, batch, in sheet order.
        cells: (reference, tested) cell types, in sheet order.
        batches: (reference, tested) read-length batches, in sheet order.
    """

    sheet: pl.DataFrame
    cells: tuple[str, str]
    batches: tuple[str, str]


def columns(*names: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of each column, the meaning from the content file."""
    return {name: (UNITS[name], text.word(PAGE, f"columns.{name}")) for name in names}


def reading(key: str, **values: Any) -> tuple[str, ...]:
    """The reading guide of a figure, each line filled as strictly as ``text.fill``."""
    lines = []
    for line in text.words(PAGE, f"figures.{key}.reading"):
        if missing := text.placeholders(line) - values.keys():
            raise text.TextError(f"content/{PAGE}.toml figures.{key}.reading needs {missing}")
        lines.append(line.format_map(values))
    return tuple(lines)


def figure_text(key: str, lead: dict[str, Any], guide: dict[str, Any] | None = None) -> FigureText:
    """Title, filled lead and reading guide of one figure."""
    return FigureText(title=text.word(PAGE, f"figures.{key}.title"),
                      lead=text.fill(PAGE, f"figures.{key}.lead", **lead),
                      reading=reading(key, **(guide or {})))  # fmt: skip


def mann_whitney(tested: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    """Two-sided Mann-Whitney U, and the rank-biserial r, above 0 when tested is higher."""
    u, p = mannwhitneyu(tested, reference, alternative="two-sided")
    return {"u": float(u), "p": float(p), "r": 2 * float(u) / (tested.size * reference.size) - 1}


def load_design(root: Path) -> tuple[Design, Path]:
    """The newest gather sheet of the made-up dataset, one row per sample."""
    path = sorted((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))[-1]
    sheet = (pl.read_parquet(path).unique("sample", keep="first", maintain_order=True)
             .select("sample", "cell_type", "condition", "batch"))  # fmt: skip
    cells, batches = (tuple(sheet[c].unique(maintain_order=True)) for c in ("cell_type", "batch"))
    if len(cells) != 2 or len(batches) != 2:  # noqa: PLR2004 - each test compares two levels
        raise ValueError(f"{path.name}: needs two cell types and two batches")
    return Design(sheet, (cells[0], cells[1]), (batches[0], batches[1])), path


def sample_metrics(rng: np.random.Generator, design: Design) -> pl.DataFrame:
    """Mapping rate, assigned reads and the detection model of every sample."""
    n = design.sheet.height
    late = (design.sheet["batch"] == design.batches[1]).to_numpy()
    tested = (design.sheet["cell_type"] == design.cells[1]).to_numpy()
    rate = np.where(late, rng.normal(*RATE[1], n), rng.normal(*RATE[0], n))
    sequenced = np.exp(rng.normal(np.log(DEPTH[0]), DEPTH[1], n))
    ceiling = np.where(tested, rng.normal(*CEILING[1], n), rng.normal(*CEILING[0], n))
    return design.sheet.with_columns(
        mapping_rate=pl.Series(rate.round(2)),
        assigned=pl.Series((sequenced * rate / 100).round(3)),
        ceiling=pl.Series(ceiling), half=pl.Series(rng.normal(*HALF, n)),
    )  # fmt: skip


def curves(rng: np.random.Generator, metrics: pl.DataFrame) -> tuple[pl.DataFrame, float]:
    """Genes detected at subsampled depths per sample, and the depth every sample reaches."""
    common = float(metrics["assigned"].min())
    rows = []
    for row in metrics.iter_rows(named=True):
        # Rounded before the union: the smallest library's own depth is the common depth.
        depths = np.union1d(np.round(np.geomspace(0.02, 1.0, STEPS) * row["assigned"], 3), [common])
        found = row["ceiling"] * depths / (depths + row["half"]) + rng.normal(0, 30, depths.size)
        found = np.maximum.accumulate(found).round().astype(int)
        labels = (row["sample"], row["cell_type"], row["condition"], row["batch"])
        rows += [(*labels, float(d), int(g)) for d, g in zip(depths, found, strict=True)]
    schema = ["sample", "cell_type", "condition", "batch", "depth", "genes"]
    frame = pl.DataFrame(rows, orient="row", schema=schema)
    return frame, common


def mapping_figure(writer: Writer, design: Design, metrics: pl.DataFrame) -> FigureBlock:
    """Mapping rate per sample by read-length batch, tested between the batches."""
    order = {b: i for i, b in enumerate(design.batches)}
    frame = (metrics.select("sample", "cell_type", "condition", "batch", "mapping_rate")
             .sort(pl.col("batch").replace_strict(order), maintain_order=True))  # fmt: skip
    rates = [frame.filter(pl.col("batch") == b)["mapping_rate"].to_numpy() for b in design.batches]
    test = mann_whitney(rates[1], rates[0])
    medians = "; ".join(
        text.fill(PAGE, "figures.mapping.median", batch=b, rate=float(np.median(r)), n=r.size)
        for b, r in zip(design.batches, rates, strict=True)
    )
    data = StripData(
        frame=frame, columns=columns(*frame.columns), category="batch", value="mapping_rate",
        category_label=text.word(PAGE, "figures.mapping.x"),
        value_label=text.word(PAGE, "figures.mapping.y"), colour="cell_type",
        colours=category_colours(list(design.cells)), hover=("sample", "condition"),
        caption=text.fill(PAGE, "figures.mapping.caption", n=frame.height, tested=design.batches[1],
                          reference=design.batches[0], **test),
    )  # fmt: skip
    return writer.figure(data, "mapping_rate", figure_text("mapping", {"medians": medians}))


def detection_figure(writer: Writer, design: Design, frame: pl.DataFrame, common: float,
                     at_depth: pl.DataFrame) -> FigureBlock:  # fmt: skip
    """Detection curves, one per sample, and a test of detection at the common depth."""
    found = [at_depth.filter(pl.col("cell_type") == c)["genes_at_depth"].to_numpy()
             for c in design.cells]  # fmt: skip
    test = mann_whitney(found[1], found[0])
    top = float(frame["genes"].max())
    names = {"reference": design.cells[0], "tested": design.cells[1], "depth": common}
    data = CurveData(
        frame=frame, columns=columns(*frame.columns), series="sample", x="depth", y="genes",
        x_label=text.word(PAGE, "figures.detection.x"),
        y_label=text.word(PAGE, "figures.detection.y"), colour="cell_type",
        colours=category_colours(list(design.cells)), hover=("condition", "batch"),
        thresholds=(Line(text.word(PAGE, "figures.detection.depth_line"), (common, common),
                         (0.0, top)),),
        caption=text.fill(PAGE, "figures.detection.caption", min_reads=MIN_READS,
                          n=at_depth.height, **names, **test),
    )  # fmt: skip
    lead = {**names, "ref_genes": int(np.median(found[0])), "test_genes": int(np.median(found[1]))}
    block = figure_text("detection", lead, {"min_reads": MIN_READS})
    return writer.figure(data, "gene_detection", block)


def correlation(rng: np.random.Generator, design: Design, assigned: np.ndarray) -> tuple:
    """Spearman rho between samples over the detected genes of a made-up count matrix."""
    n = design.sheet.height
    log_mean = rng.normal(7.0, 2.2, GENES)
    shifts = np.zeros((GENES, n))
    for column, share, sd in EFFECTS:
        effect = np.where(rng.random(GENES) < share, rng.normal(0, sd, GENES), 0.0)
        shifts += np.outer(effect, (design.sheet[column] != design.sheet[column][0]).to_numpy())
    noise = rng.normal(0, 0.3, (GENES, n))
    counts = rng.poisson(2 ** (log_mean[:, None] + shifts + noise) * assigned / np.median(assigned))
    detected = counts[counts.mean(axis=1) >= MIN_READS]
    return spearmanr(detected).statistic, detected.shape[0]


def permutation(rho: np.ndarray, labels: np.ndarray, rng: np.random.Generator) -> tuple:
    """Mean rho within minus between label groups, and its label-permutation p."""
    upper = np.triu_indices(labels.size, 1)
    values = rho[upper]

    def within_minus_between(shuffled: np.ndarray) -> float:
        """The statistic for one labelling."""
        same = shuffled[upper[0]] == shuffled[upper[1]]
        return float(values[same].mean() - values[~same].mean())

    observed = within_minus_between(labels)
    null = np.array([within_minus_between(rng.permutation(labels)) for _ in range(PERMUTATIONS)])
    same = labels[upper[0]] == labels[upper[1]]
    medians = {"within": float(np.median(values[same])), "between": float(np.median(values[~same]))}
    return observed, (1 + int((null >= observed).sum())) / (1 + PERMUTATIONS), medians


def correlation_figure(writer: Writer, design: Design, metrics: pl.DataFrame,
                       rng: np.random.Generator) -> FigureBlock:  # fmt: skip
    """The sample correlation heatmap, clustered by average linkage with optimal leaf order."""
    rho, genes = correlation(rng, design, metrics["assigned"].to_numpy())
    diff, p, medians = permutation(rho, design.sheet["cell_type"].to_numpy(), rng)
    samples = design.sheet["sample"].to_list()
    distance = squareform(1 - rho, checks=False)
    tree = linkage(distance, "average", optimal_ordering=True)
    order = tuple(samples[i] for i in leaves_list(tree))
    levels, swatches = rho_levels(rho)
    frame = pl.DataFrame({"sample_row": [s for s in samples for _ in samples],
                          "sample_column": samples * len(samples),
                          "rho": rho.ravel().round(4), "level": levels})  # fmt: skip
    batch_colours = dict(zip(design.batches, (DATA_COLOURS[2], DATA_COLOURS[6]), strict=True))
    tracks = tuple(
        Track(text.word(PAGE, f"figures.correlation.{column}"),
              dict(zip(samples, design.sheet[column], strict=True)), colours)
        for column, colours in (("cell_type", category_colours(list(design.cells))),
                                ("batch", batch_colours))
    )  # fmt: skip
    axis = text.word(PAGE, "figures.correlation.axis")
    data = HeatmapData(
        frame=frame, columns=columns(*frame.columns), row="sample_row", column="sample_column",
        value="rho", level="level", row_label=axis, column_label=axis, swatches=swatches,
        caption=text.fill(PAGE, "figures.correlation.caption", n=len(samples), genes=genes,
                          diff=diff, p=p, perms=PERMUTATIONS),
        row_order=order, column_order=order, tracks=tracks,
    )  # fmt: skip
    return writer.figure(data, "sample_correlation", figure_text("correlation", medians))


def units_frame(rng: np.random.Generator, depth_millions: float) -> tuple[pl.DataFrame, int]:
    """One sample's genes: effective length, estimated reads and TPM; a random subset drawn."""
    median, sd, floor, ceiling = LENGTH
    length = np.clip(np.exp(rng.normal(np.log(median), sd, UNIVERSE)), floor, ceiling).round()
    weight = np.exp(rng.normal(0, 2.0, UNIVERSE)) * length
    reads = rng.poisson(depth_millions * 1e6 * weight / weight.sum()).astype(float)
    tpm = 1e6 * (reads / length) / (reads / length).sum()
    seen = np.flatnonzero(reads > 0)
    drawn = np.sort(rng.choice(seen, size=min(DRAWN, seen.size), replace=False))
    return pl.DataFrame({
        "gene": [f"G{i:05d}" for i in drawn], "tpm": tpm[drawn].round(4), "reads": reads[drawn],
        "effective_length": length[drawn],
    }), int(seen.size)  # fmt: skip


def units_figure(writer: Writer, metrics: pl.DataFrame, rng: np.random.Generator) -> FigureBlock:
    """Estimated reads against TPM in the median-depth sample, coloured by effective length."""
    row = metrics.sort("assigned").row(metrics.height // 2, named=True)
    frame, total = units_frame(rng, row["assigned"])
    bounds = {"short": SHORT_NT, "long": LONG_NT}
    names = {k: text.fill(PAGE, f"figures.units.{k}", **bounds)
             for k in ("short", "medium", "long")}  # fmt: skip
    length = pl.col("effective_length")
    frame = frame.with_columns(
        length_class=pl.when(length < SHORT_NT).then(pl.lit(names["short"]))
        .when(length > LONG_NT).then(pl.lit(names["long"])).otherwise(pl.lit(names["medium"]))
    )  # fmt: skip
    per_tpm = frame.group_by("length_class").agg((pl.col("reads") / pl.col("tpm")).median())
    ratio = dict(per_tpm.iter_rows())
    colours = dict(zip(names.values(), (DATA_COLOURS[2], DATA_COLOURS[4], DATA_COLOURS[6]),
                       strict=True))  # fmt: skip
    data = ScatterData(
        frame=frame, columns=columns(*frame.columns), x="tpm", y="reads",
        x_label=text.word(PAGE, "figures.units.x"), y_label=text.word(PAGE, "figures.units.y"),
        colour="length_class", colours=colours, hover=("gene", "effective_length"),
        log_x=True, log_y=True,
        caption=text.fill(PAGE, "figures.units.caption", sample=row["sample"], n=frame.height,
                          total=total),
    )  # fmt: skip
    lead = {"sample": row["sample"], "ratio": ratio[names["long"]] / ratio[names["short"]]}
    return writer.figure(data, "reads_against_tpm", figure_text("units", lead))


def tool_means(rng: np.random.Generator, n: int) -> np.ndarray:
    """Mean reads per sample of every gene by Salmon, STAR and featureCounts, one row per tool.

    Salmon keeps multi-mapping reads; STAR's gene counts and featureCounts drop them and the
    reads on overlapping genes, each with its own small noise.
    """
    mean = np.exp(rng.normal(2.5, 2.5, TOOL_GENES))
    multi = np.where(rng.random(TOOL_GENES) < PARALOGS, rng.beta(3, 2, TOOL_GENES),
                     rng.beta(1, 40, TOOL_GENES))  # fmt: skip
    ambiguous = np.where(rng.random(TOOL_GENES) < OVERLAPPING,
                         rng.uniform(0.1, 0.6, TOOL_GENES), 0.0)  # fmt: skip
    shares = (
        np.ones(TOOL_GENES),
        (1 - multi) * (1 - ambiguous),
        (1 - multi) * (1 - ambiguous * rng.uniform(0.8, 1.2, TOOL_GENES)),
    )
    return np.array([rng.poisson(n * mean * s * np.exp(rng.normal(0, 0.1, TOOL_GENES))) / n
                     for s in shares])  # fmt: skip


def tools_figure(writer: Writer, n: int, rng: np.random.Generator) -> tuple[FigureBlock, int, int]:
    """Genes detected by one or two of the three tools; return the block, shared and total."""
    means = tool_means(rng, n)
    found = means >= MIN_READS
    tools = text.words(PAGE, "figures.tools.tools")
    salmon = {True: text.word(PAGE, "figures.tools.salmon_in"),
              False: text.word(PAGE, "figures.tools.salmon_out")}  # fmt: skip
    rows = []
    for combo in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)):
        genes = np.all(found == np.array(combo, dtype=bool)[:, None], axis=0)
        named = [t for t, keep in zip(tools, combo, strict=True) if keep]
        region = (
            text.fill(PAGE, "figures.tools.only", tool=named[0]) if len(named) == 1
            else text.fill(PAGE, "figures.tools.pair", first=named[0], second=named[1])
        )  # fmt: skip
        typical = float(np.median(means[0, genes])) if genes.any() else None
        rows.append((region, salmon[bool(combo[0])], int(genes.sum()), typical))
    schema = ["region", "salmon", "n_genes", "median_reads"]
    frame = pl.DataFrame(rows, schema=schema, orient="row").sort("n_genes", descending=True,
                                                                maintain_order=True)  # fmt: skip
    shared, total = int(found.all(axis=0).sum()), int(found.any(axis=0).sum())
    data = BarData(
        frame=frame, columns=columns(*frame.columns), category="region", value="n_genes",
        category_label=text.word(PAGE, "figures.tools.x"),
        value_label=text.word(PAGE, "figures.tools.y"), series="salmon",
        colours=dict(zip(salmon.values(), (DATA_COLOURS[4], DATA_COLOURS[5]), strict=True)),
        caption=text.fill(PAGE, "figures.tools.caption", min_reads=MIN_READS, n=n,
                          shared=shared, total=total),
        horizontal=True, stacked=True, hover=("median_reads",),
    )  # fmt: skip
    lead = {"disagree": int(frame["n_genes"].sum()), "salmon_only": int(rows[0][2])}
    return writer.figure(data, "tool_detection", figure_text("tools", lead)), shared, total


def cards() -> tuple[ParameterCard, ...]:
    """One card per option ``ctx quant`` will have."""
    planned = text.raw(PAGE, "planned")
    return tuple(ParameterCard(flag, card["meaning"], card["default"], card["effect"])
                 for flag, card in planned.items())  # fmt: skip


def page(root: Path, dataset: str, day: date) -> Path:
    """Write the M4 preview page, its figures and its sample table; return the page.

    Args:
        root: The ctx root holding the made-up dataset's gather outputs.
        dataset: The made-up dataset's slug.
        day: The date every output name carries.

    Returns:
        The page, ``outputs/04_quant/<day>_quant_<dataset>_overview_v1.html``.
    """
    started = time.perf_counter()
    rng = np.random.default_rng(SEED)
    design, sheet_path = load_design(root)
    releases = {"input:samplesheet": sheet_path.name,
                "made-up data": f"numpy default_rng({SEED}), {Path(__file__).name}"}  # fmt: skip
    writer = Writer(root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset, day,
                    releases)  # fmt: skip
    metrics = sample_metrics(rng, design)
    frame, common = curves(rng, metrics)
    at_depth = frame.filter(pl.col("depth") == common).select(
        "sample", "cell_type", genes_at_depth="genes")  # fmt: skip
    if at_depth["sample"].to_list() != design.sheet["sample"].to_list():
        raise ValueError("the common depth is not one point of every curve")
    table = metrics.join(at_depth.drop("cell_type"), on="sample", maintain_order="left").select(
        "sample", "cell_type", "condition", "batch", "mapping_rate", "assigned", "genes_at_depth"
    )
    writer.table(table, "sample_metrics", columns(*table.columns))
    samples = Group(text.word(PAGE, "groups.samples"), (
        mapping_figure(writer, design, metrics),
        detection_figure(writer, design, frame, common, at_depth),
        correlation_figure(writer, design, metrics, rng)))  # fmt: skip
    tools, shared, total = tools_figure(writer, design.sheet.height, rng)
    genes = Group(text.word(PAGE, "groups.genes"), (units_figure(writer, metrics, rng), tools))
    keys = (
        KeyNumber(text.word(PAGE, "keys.samples"), f"{design.sheet.height}"),
        KeyNumber(text.word(PAGE, "keys.mapping"), f"{metrics['mapping_rate'].median():.1f} %",
                  text.fill(PAGE, "keys.mapping_detail", n=design.sheet.height)),
        KeyNumber(text.word(PAGE, "keys.shared"), f"{shared:,}",
                  text.fill(PAGE, "keys.shared_detail", total=total)),
    )  # fmt: skip
    usage, budget = runtime.measure(started), runtime.budget(MODULE)
    report = ModulePage(
        title=text.fill(PAGE, "title", dataset=text.word(*SHOWN)), **text.skeleton(PAGE),
        parameters=cards(), results=(text.word("ui", "preview.note"), samples, genes),
        key_numbers=keys,
        reproduce=Reproduction(
            command=f"ctx quant --dataset {dataset}", wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self, peak_rss_children=usage.peak_rss_children,
            profile=text.word(PAGE, "profile"),
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=(sheet_path.name,),
        ),
    )  # fmt: skip
    return write(report, writer.out / f"{writer.stem('overview')}.html")
