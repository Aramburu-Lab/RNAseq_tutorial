"""M7's figures, each from one contract: missingness first, then correlation or fold change.

Every sentence comes from ``content/proteo.toml``; numbers are filled from the run.
"""

import os
from collections.abc import Mapping
from dataclasses import replace

import numpy as np
import polars as pl

from ctx.de.run import Writer
from ctx.proteo import folds, measure, tables
from ctx.report import text
from ctx.report.builder import FigureBlock, FigureText
from ctx.viz.bar import BarData
from ctx.viz.contour import DENSITY_MIN_POINTS, HDR_MASSES
from ctx.viz.curve import CurveData
from ctx.viz.density import DensityData
from ctx.viz.heatmap import HeatmapData, Track, diverging_levels
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import DATA_COLOURS, PROJECTOR_UNSAFE, ROLES, SPENT_GREYS

PAGE = "proteo"
# Category colours, the projector-unsafe yellow left out; at most seven levels (CONVENTIONS.md).
COLOURS = tuple(c for c in DATA_COLOURS if c not in PROJECTOR_UNSAFE)
# One colour per gene class on every figure, by role: Other is the reference population.
# Complex and metabolic take the palette's orange and green, which no role holds.
CLASS_COLOURS = {"granule": ROLES["comparison"], "complex": DATA_COLOURS[0],
                 "metabolic": DATA_COLOURS[2], measure.OTHER: ROLES["reference"]}  # fmt: skip
# F4's groups by where a gene was called; neither is the reference population.
KIND_COLOURS = {"neither": ROLES["reference"], "rna": DATA_COLOURS[0],
                "protein": DATA_COLOURS[6], "both": DATA_COLOURS[2]}  # fmt: skip
# F5's track colours: never the heat's two role colours, or a PAM50 Normal-like cell read as a
# low z (C3.4); a fifth level takes the lighter grey, two lightnesses being the colours left.
TRACK_COLOURS = (*(c for c in COLOURS if c not in ROLES.values()), SPENT_GREYS[1])
# The two highest-density regions, as the captions name them.
HDR = {"inner": f"{HDR_MASSES[0] * 100:.0f} %", "outer": f"{HDR_MASSES[1] * 100:.0f} %"}


def columns(*names: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of each named column, from the content file."""
    table = text.raw(PAGE, "columns")
    return {name: (table[name][0], table[name][1]) for name in names}


def lines(key: str, **values: object) -> tuple[str, ...]:
    """A figure's reading guide, each line filled from the run's numbers."""
    return tuple(line.format_map(values) for line in text.words(PAGE, key))


def class_name(name: str) -> str:
    """A class's shown name."""
    return text.word(PAGE, f"classes.{name}")


def class_colours(frame: pl.DataFrame) -> dict[str, str]:
    """Each shown class name of the frame to its colour."""
    names = {class_name(k): v for k, v in CLASS_COLOURS.items()}
    return {c: names[c] for c in frame["gene_class"].unique().sort()}


def short(names: tuple[str, ...]) -> list[str]:
    """Sample names without the words every one of them shares at either end.

    ``CopyNumber_Neutrophil_01_steady-state`` becomes ``Neutrophil_01``: whole underscore-
    separated words only, so ``XA01`` and ``XA02`` keep their shared ``X``.
    """
    parts = [n.split("_") for n in names]
    head = len(os.path.commonprefix(parts))
    tail = len(os.path.commonprefix([p[::-1] for p in parts]))
    if len(names) < 2 or head + tail >= min(len(p) for p in parts):  # noqa: PLR2004 - a pair
        return list(names)
    return ["_".join(p[head : len(p) - tail]) for p in parts]


def quantified(writer: Writer, layer: tables.Layer, colour: list[str], hover: list[str],
               minimum: float) -> FigureBlock:  # fmt: skip
    """F1: genes with a protein value in each sample, bars coloured by group."""
    present = ~np.isnan(layer.values)
    frame = pl.DataFrame({"sample": short(layer.samples), "group": colour, "track": hover,
                          "n_quantified": present.sum(axis=0).astype(np.int64)})  # fmt: skip
    frame = frame.sort("track", "group", maintain_order=True)
    levels = list(dict.fromkeys(sorted(set(colour))))
    seen = present.any(axis=1)
    n_kept = int(measure.enough(present.sum(axis=1), present.shape[1], minimum).sum())
    data = BarData(
        frame=frame, columns=columns("sample", "group", "track", "n_quantified"),
        category="sample", value="n_quantified",
        category_label=text.word(PAGE, "figures.quantified.category"),
        value_label=text.word(PAGE, "figures.quantified.value"), series="group",
        colours=dict(zip(levels, COLOURS, strict=False)) if len(levels) <= len(COLOURS)
        else dict.fromkeys(levels),
        caption=text.fill(PAGE, "figures.quantified.caption", n_samples=frame.height,
                          n_groups=len(levels), n_ids=int(seen.sum()),
                          n_all=int(present.all(axis=1).sum())),
        hover=("track",),
    )  # fmt: skip
    return writer.figure(data, "quantified", FigureText(
        title=text.word(PAGE, "figures.quantified.title"),
        lead=text.fill(PAGE, "figures.quantified.lead", n_kept=n_kept, n_ids=int(seen.sum()),
                       min_pct=round(100 * minimum)),
        reading=lines("figures.quantified.reading")))  # fmt: skip


def _cumulative(frame: pl.DataFrame, group: str, value: str) -> pl.DataFrame:
    """Each group's values sorted, with the share at or below each one."""
    return frame.sort(group, value).with_columns(
        share=pl.col(value).cum_count().over(group) / pl.len().over(group)
    )


def coverage(writer: Writer, genes: pl.DataFrame, abundance: np.ndarray, minimum: float,
             abundance_label: str) -> FigureBlock:  # fmt: skip
    """F2: the quantified share per protein, one cumulative curve per RNA abundance tertile."""
    tertile = measure.tertiles(abundance)
    names = [text.fill(PAGE, "figures.coverage.tertile", n=t) for t in tertile]
    frame = _cumulative(genes.select("gene", "quantified_share").with_columns(
        tertile=pl.Series(names)), "tertile", "quantified_share")  # fmt: skip
    tested = genes["tested"].to_numpy()
    passing = {t: float(tested[tertile == t].mean()) for t in (1, 3)}
    data = CurveData(
        frame=frame, columns=columns("gene", "quantified_share", "tertile", "share"),
        series="tertile", x="quantified_share", y="share",
        x_label=text.word(PAGE, "figures.coverage.x"),
        y_label=text.word(PAGE, "figures.coverage.y"),
        colour="tertile", colours={}, distinct=True, hover=("gene",),
        caption=text.fill(PAGE, "figures.coverage.caption", n_genes=genes.height,
                          abundance=abundance_label),
        thresholds=(Line(text.word(PAGE, "figures.coverage.line"), (minimum, minimum), (0.0, 1.0),
                         dash=True),),
    )  # fmt: skip
    return writer.figure(data, "coverage", FigureText(
        title=text.word(PAGE, "figures.coverage.title"),
        lead=text.fill(PAGE, "figures.coverage.lead", low_pct=100 * passing[1],
                       high_pct=100 * passing[3]),
        reading=lines("figures.coverage.reading")))  # fmt: skip


def _counts(frame: pl.DataFrame) -> str:
    """Genes per class, as the caption lists them; the frame holds the shown class names."""
    found = frame.group_by("gene_class", maintain_order=True).len().rows()
    return ", ".join(text.fill(PAGE, "classes.count", name=c, n=n) for c, n in found)


def shown(frame: pl.DataFrame, gene_class: str) -> pl.DataFrame:
    """The genes a class figure draws: one class against Other, or every class."""
    kept = frame.filter(pl.col("gene_class") != "")
    if gene_class != "all":
        kept = kept.filter(pl.col("gene_class").is_in([gene_class, measure.OTHER]))
    return kept.with_columns(gene_class=pl.col("gene_class").replace_strict(
        {c: class_name(c) for c in kept["gene_class"].unique()}))  # fmt: skip


def rho_by_class(writer: Writer, tested: pl.DataFrame, values: Mapping[str, object]
                 ) -> FigureBlock:  # fmt: skip
    """F3: per-gene rho, one cumulative curve per gene class."""
    frame = _cumulative(
        tested.select("gene", "gene_class", "rho").drop_nulls(), "gene_class", "rho"
    )
    ends = (float(frame["rho"].to_numpy().min()), float(frame["rho"].to_numpy().max()))
    data = CurveData(
        frame=frame, columns=columns(*frame.columns), series="gene_class", x="rho", y="share",
        x_label=text.word(PAGE, "figures.rho.x"), y_label=text.word(PAGE, "figures.rho.y"),
        colour="gene_class", colours=class_colours(frame), hover=("gene",),
        caption=text.fill(PAGE, "figures.rho.caption", counts=_counts(frame), **values),
        thresholds=(Line(text.word(PAGE, "figures.rho.median"), ends, (0.5, 0.5), dash=True),),
    )  # fmt: skip
    return writer.figure(data, "rho_by_class", FigureText(
        title=text.word(PAGE, "figures.rho.title"),
        lead=text.fill(PAGE, "figures.rho.lead", **values),
        reading=lines("figures.rho.reading")))  # fmt: skip


def worked(writer: Writer, gene: str, frame: pl.DataFrame, values: Mapping[str, object]
           ) -> FigureBlock:  # fmt: skip
    """F6: one worked gene's RNA against protein per sample, coloured by the annotation."""
    levels = sorted(set(frame["track"]))
    data = ScatterData(
        frame=frame, columns=columns("sample", "rna", "protein", "track"), x="rna", y="protein",
        x_label=text.word(PAGE, "figures.worked.x"), y_label=text.word(PAGE, "figures.worked.y"),
        colour="track", colours=dict(zip(levels, COLOURS, strict=False)), hover=("sample",),
        caption=text.fill(PAGE, "figures.worked.caption", gene=gene, **values),
    )  # fmt: skip
    return writer.figure(data, f"worked_{gene.lower()}", FigureText(
        title=text.fill(PAGE, "figures.worked.title", gene=gene),
        lead=text.fill(PAGE, "figures.worked.lead", **values),
        reading=lines("figures.worked.reading", **values)))  # fmt: skip


def fold(writer: Writer, tested: pl.DataFrame, values: Mapping[str, object]) -> FigureBlock:
    """F7: RNA against protein fold change per gene, a density per class, the class marked."""
    frame = tested.select(
        "gene", "gene_class", "rna_lfc", "protein_lfc", "rna_padj", "protein_padj"
    ).with_columns(keep=pl.col("gene_class") != measure.OTHER)
    frame = frame.with_columns(gene_class=pl.col("gene_class").replace_strict(
        {c: class_name(c) for c in frame["gene_class"].unique()}))  # fmt: skip
    data = DensityData(
        frame=frame, columns=columns(*frame.columns), x="rna_lfc", y="protein_lfc",
        group="gene_class", x_label=text.word(PAGE, "figures.fold.x"),
        y_label=text.word(PAGE, "figures.fold.y"), colours=class_colours(frame), caption="",
        highlight="keep", highlight_label=", ".join(
            c for c in class_colours(frame) if c != class_name(measure.OTHER)),
        hover=("gene", "rna_padj", "protein_padj"),
    )  # fmt: skip
    data = replace(data, caption=text.fill(PAGE, "figures.fold.caption", **values, **HDR,
                                           bandwidth=data.bandwidth()))  # fmt: skip
    return writer.figure(data, "fold_change", FigureText(
        title=text.fill(PAGE, "figures.fold.title", **values),
        lead=text.fill(PAGE, "figures.fold.lead", **values),
        reading=lines("figures.fold.reading", minimum=DENSITY_MIN_POINTS, **HDR)))  # fmt: skip


def contrast_density(writer: Writer, made: folds.Folds, values: Mapping[str, object]
                     ) -> FigureBlock:  # fmt: skip
    """F4: one contrast's RNA against protein fold change on the same tumours.

    A density per group of where a gene was called; the genes called at both layers are drawn
    as points over the contours.
    """
    names = {k: text.word(PAGE, f"figures.contrast.kinds.{k}") for k in folds.KINDS}
    frame = made.genes.select(
        "gene", "rna_lfc", "protein_lfc", "rna_padj", "protein_padj", "rna_call", "protein_call",
        called=pl.col("kind") == folds.KINDS[0], kind=pl.col("kind").replace_strict(names),
    )  # fmt: skip
    present = set(frame["kind"])
    colours = {names[k]: c for k, c in KIND_COLOURS.items() if names[k] in present}
    data = DensityData(
        frame=frame, columns=columns(*frame.columns), x="rna_lfc", y="protein_lfc", group="kind",
        x_label=text.word(PAGE, "figures.contrast.x"),
        y_label=text.word(PAGE, "figures.contrast.y"), colours=colours, caption="",
        highlight="called", highlight_label=names[folds.KINDS[0]],
        hover=("gene", "rna_padj", "protein_padj"),
    )  # fmt: skip
    data = replace(data, caption=text.fill(PAGE, "figures.contrast.caption", **values, **HDR,
                                           bandwidth=data.bandwidth()))  # fmt: skip
    return writer.figure(data, f"contrast_{made.contrast.name}", FigureText(
        title=text.fill(PAGE, "figures.contrast.title", **values),
        lead=text.fill(PAGE, "figures.contrast.lead", **values),
        reading=lines("figures.contrast.reading", minimum=DENSITY_MIN_POINTS, **HDR,
                      **values)))  # fmt: skip


def zscores(values: np.ndarray) -> np.ndarray:
    """Each row centred and scaled over its present values; a constant row is 0 / 0, NaN."""
    with np.errstate(invalid="ignore", divide="ignore"):
        centred = values - np.nanmean(values, axis=1, keepdims=True)
        return centred / np.nanstd(values, axis=1, keepdims=True)


def heatmap_genes(genes: pl.DataFrame, named: list[str], top: int) -> list[str]:
    """F5's genes: the named ones that were tested, then the ``top`` highest and lowest rho."""
    tested = genes.filter(pl.col("tested") & pl.col("rho").is_not_null())
    first = [g for g in dict.fromkeys(named) if g in set(tested["gene"])]
    rest = tested.filter(~pl.col("gene").is_in(first)).sort("rho", "gene", descending=[True, False])
    ends = [*rest["gene"].head(top), *rest["gene"].tail(top)]
    return [*first, *dict.fromkeys(ends)]


def track_colours(levels: list[str]) -> dict[str, str]:
    """A heatmap track's colours, the heat's own two role colours last.

    A track with more levels than the palette holds (CPTAC's 22 plexes) takes two greys in turn
    (D-063): lightnesses, not hues, so no plex reads as a category of the data colours. Its
    levels must be one run each in the column order (``sort_keys``), or a boundary between two
    plexes of one grey is lost; the hover names the plex.
    """
    if len(levels) <= len(TRACK_COLOURS):
        return dict(zip(levels, TRACK_COLOURS, strict=False))
    return {level: SPENT_GREYS[i % 2] for i, level in enumerate(levels)}


def track(name: str, values: Mapping[str, str], levels: list[str]) -> Track:
    """F5's strip for one annotation column: a greyed one names its alternation in one entry."""
    colours = track_colours(levels)
    if len(levels) <= len(TRACK_COLOURS):
        return Track(name, values, colours)
    words = text.fill(PAGE, "figures.paired.alternating", track=name, n=len(levels))
    return Track(name, values, colours, key=((words, SPENT_GREYS[0]),))


def by_label(values: pl.Series) -> list[str]:
    """A track's levels in display order: plex 2 before plex 10, by length, then text."""
    return sorted(set(values.cast(pl.String)), key=lambda v: (len(v), v))


def sort_order(annotation: pl.DataFrame, tracks: tuple[tuple[str, str], ...]
               ) -> list[tuple[str, str]]:  # fmt: skip
    """The tracks the columns sort by, greyed ones first, so each of their levels is one run.

    Sorted subtype first, CPTAC's plexes broke into 63 runs over 122 tumours, and no two greys
    given per plex could mark every boundary (the plexes' neighbour graph has an odd cycle).
    """
    greyed = [t for t in tracks if len(by_label(annotation[t[1]])) > len(TRACK_COLOURS)]
    return [*greyed, *(t for t in tracks if t not in greyed)]


def paired(writer: Writer, measured: tuple[pl.DataFrame, np.ndarray, np.ndarray],
           annotation: pl.DataFrame, chosen: list[str],
           tracks: tuple[tuple[str, str], ...]) -> FigureBlock:  # fmt: skip
    """F5: each gene's RNA row above its protein row, z-scored, tumour by tumour.

    Args:
        writer: The module's writer.
        measured: The per-gene table, and RNA and protein, genes by paired samples, in its order.
        annotation: The paired samples' annotation, in the matrices' column order.
        chosen: The genes drawn, top to bottom.
        tracks: (shown name, annotation column) per track; columns sort by them in turn.
    """
    genes, rna, protein = measured
    row_of = {g: i for i, g in enumerate(genes["gene"])}
    samples = short(tuple(annotation["sample"]))
    sort = sort_order(annotation, tracks)
    by = [k for _, c in sort for k in (pl.col(c).cast(pl.String).str.len_chars(),
                                       pl.col(c).cast(pl.String))]  # fmt: skip
    order = (annotation.with_columns(shown=pl.Series(samples), at=pl.int_range(pl.len()))
             .sort([*by, pl.col("shown")]))  # fmt: skip
    at, shown = order["at"].to_list(), order["shown"].to_list()
    labels, values = [], []
    for gene in chosen:
        for key, layer in (("rna_row", rna), ("protein_row", protein)):
            labels.append(text.fill(PAGE, f"figures.paired.{key}", gene=gene))
            values.append(layer[row_of[gene], at])
    z = zscores(np.array(values)).round(4)
    flat = [None if np.isnan(v) else float(v) for v in z.ravel()]
    levels, swatches = diverging_levels(flat)
    frame = pl.DataFrame({"row": [r for r in labels for _ in shown], "sample": shown * len(labels),
                          "z": flat, "level": levels}, strict=False)  # fmt: skip
    rho = dict(genes.select("gene", "rho").iter_rows())
    strips = []
    for name, column in tracks:
        values_of = dict(zip(shown, order[column].cast(pl.String), strict=True))
        strips.append(track(name, values_of, by_label(order[column])))
    counts = {"n_genes": len(chosen), "n_samples": len(shown),
              "missing": sum(v is None for v in flat),
              "tracks": ", ".join(n for n, _ in tracks),
              "sorted": ", ".join(n for n, _ in sort)}  # fmt: skip
    data = HeatmapData(
        frame=frame, columns=columns("row", "sample", "z", "level"), row="row", column="sample",
        value="z", level="level", row_label=text.word(PAGE, "figures.paired.row"),
        column_label=text.word(PAGE, "figures.paired.column"), swatches=swatches,
        caption=text.fill(PAGE, "figures.paired.caption", **counts), row_order=tuple(labels),
        column_order=tuple(shown), tracks=tuple(strips),
        row_stats={label: rho[chosen[i // 2]] for i, label in enumerate(labels)},
        row_stats_label=text.word(PAGE, "figures.paired.rho"), row_stats_format=".2f",
        column_title=text.fill(PAGE, "figures.paired.column_title", **counts),
    )  # fmt: skip
    return writer.figure(data, "paired_zscores", FigureText(
        title=text.word(PAGE, "figures.paired.title"),
        lead=text.fill(PAGE, "figures.paired.lead", **counts),
        reading=lines("figures.paired.reading", **counts)))  # fmt: skip


def shift(writer: Writer, tested: pl.DataFrame, values: Mapping[str, object]) -> FigureBlock:
    """F8: protein minus RNA fold change, one cumulative curve per class."""
    frame = _cumulative(tested.select("gene", "gene_class", "difference"), "gene_class",
                        "difference")  # fmt: skip
    data = CurveData(
        frame=frame, columns=columns(*frame.columns), series="gene_class", x="difference",
        y="share", x_label=text.word(PAGE, "figures.shift.x"),
        y_label=text.word(PAGE, "figures.shift.y"), colour="gene_class",
        colours=class_colours(frame), hover=("gene",),
        caption=text.fill(PAGE, "figures.shift.caption", counts=_counts(frame), **values),
    )  # fmt: skip
    return writer.figure(data, "shift_by_class", FigureText(
        title=text.word(PAGE, "figures.shift.title"),
        lead=text.fill(PAGE, "figures.shift.lead", **values),
        reading=lines("figures.shift.reading")))  # fmt: skip
