"""M7's figures, each from one contract: missingness first, then correlation or fold change.

Every sentence comes from ``content/proteo.toml``; numbers are filled from the run.
"""

import os
from collections.abc import Mapping

import numpy as np
import polars as pl

from ctx.de.run import Writer
from ctx.proteo import measure, tables
from ctx.report import text
from ctx.report.builder import FigureBlock, FigureText
from ctx.viz.bar import BarData
from ctx.viz.curve import CurveData
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import DATA_COLOURS, PROJECTOR_UNSAFE, ROLES

PAGE = "proteo"
# Category colours, the projector-unsafe yellow left out; at most seven levels (CONVENTIONS.md).
COLOURS = tuple(c for c in DATA_COLOURS if c not in PROJECTOR_UNSAFE)
# One colour per gene class on every figure, by role: Other is the reference population.
# Complex and metabolic take the palette's orange and green, which no role holds.
CLASS_COLOURS = {"granule": ROLES["comparison"], "complex": DATA_COLOURS[0],
                 "metabolic": DATA_COLOURS[2], measure.OTHER: ROLES["reference"]}  # fmt: skip


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
    """F7: RNA against protein fold change per gene, the declared class marked."""
    picked = tested.select("gene", "gene_class", "rna_lfc", "protein_lfc", "rna_padj",
                           "protein_padj")  # fmt: skip
    frame = picked.with_columns(keep=pl.col("gene_class") != measure.OTHER)
    frame = frame.with_columns(gene_class=pl.col("gene_class").replace_strict(
        {c: class_name(c) for c in frame["gene_class"].unique()}))  # fmt: skip
    spread = np.abs(frame.select("rna_lfc", "protein_lfc").to_numpy())
    top = float(np.nanmax(spread)) * 1.05
    zero = text.word(PAGE, "figures.fold.zero")
    marked = [c for c in frame["gene_class"].unique().sort() if c != class_name(measure.OTHER)]
    data = ScatterData(
        frame=frame, columns=columns(*frame.columns), x="rna_lfc", y="protein_lfc",
        x_label=text.word(PAGE, "figures.fold.x"), y_label=text.word(PAGE, "figures.fold.y"),
        colour="gene_class",
        colours={**{c: v for c, v in class_colours(frame).items() if c in marked},
                 class_name(measure.OTHER): None},
        caption=text.fill(PAGE, "figures.fold.caption", **values),
        hover=("gene", "rna_padj", "protein_padj"), keep="keep",
        lines=(Line(text.word(PAGE, "figures.fold.identity"), (-top, top), (-top, top)),
               Line(zero, (-top, top), (0.0, 0.0), dash=True),
               Line(zero, (0.0, 0.0), (-top, top), dash=True)),
    )  # fmt: skip
    return writer.figure(data, "fold_change", FigureText(
        title=text.fill(PAGE, "figures.fold.title", **values),
        lead=text.fill(PAGE, "figures.fold.lead", **values),
        reading=lines("figures.fold.reading")))  # fmt: skip


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
