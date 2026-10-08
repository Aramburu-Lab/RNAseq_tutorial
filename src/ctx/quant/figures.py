"""M4's figure contracts, built from ``ctx.quant.tables``' frames; every word from content.

F1 Salmon's fragment fates per sample, F2 genes detected against depth, F3 Spearman between
samples, F4 reads against TPM, F6 the featured contrast's fold change from estimated and from
length-scaled counts, F7 isoforms expressed against annotated. F1 and F7 filter by sheet
column and gene type (D-051); F2, F3, F4 and F6 have no control yet, as a curve, heatmap or
scatter takes none (M4.14 in the M4 plan). F5 and F8 need runs that do not exist yet.
"""

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import polars as pl
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from ctx.quant import tables as tab
from ctx.report import text
from ctx.viz.annotate import Caption
from ctx.viz.bar import MAX_TICK_LABELS, BarData
from ctx.viz.curve import CurveData
from ctx.viz.heatmap import HeatmapData, Swatch, Track
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import BAND_COLOURS, CIVIDIS, DATA_COLOURS, category_colours

PAGE = "quant"
# Each fate by its role: a mapped fragment is good, one with no hit bad; a decoy's and a
# below-score fragment take the colours the QC bands do not use for good or bad.
FATE_COLOURS = {"mapped": BAND_COLOURS["good"], "decoy": DATA_COLOURS[1],
                "below_score": BAND_COLOURS["usable"], "no_hit": BAND_COLOURS["bad"]}  # fmt: skip
BINS = 6  # rho bins of equal width over the off-diagonal range
SHORT_NT, LONG_NT = 1_000, 4_000
UNITS = {
    "sample": "none", "cell_type": "none", "condition": "none", "batch": "none",
    "fate": "none", "reads": "fragments", "percent": "%", "processed": "fragments",
    "depth": "reads", "genes": "genes", "common": "none", "sample_row": "none",
    "sample_column": "none", "rho": "dimensionless", "level": "none", "gene_id": "none",
    "gene_name": "none", "tpm": "TPM", "counts": "reads", "length": "nt", "length_class": "none",
    "estimated": "log2", "length_scaled": "log2", "difference": "log2", "gene_type": "none",
    "annotated_bin": "none", "share": "%", "n_genes": "genes", "mapped": "fragments",
    "decoy": "fragments", "below_score": "fragments", "no_hit": "fragments",
    "strand_bias": "fraction", "fragment_mean": "nt", "fragment_sd": "nt", "targets": "sequences",
    "decoys": "sequences", "length_scaled_sum": "reads", "relative_difference": "fraction",
    "n_reference": "samples", "n_tested": "samples", "u": "dimensionless", "p": "probability",
    "p_bh": "probability", "r": "dimensionless", "annotated": "transcripts",
    "expressed": "transcripts",
}  # fmt: skip


def columns(frame: pl.DataFrame) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of a figure's frame, the meaning from content."""
    return {c: (UNITS[c], text.word(PAGE, f"columns.{c}")) for c in frame.columns}


def fate_bars(frame: pl.DataFrame, groups: Sequence[str], caption: str) -> BarData:
    """F1: every sample's fragments by fate, as a share of those processed, filtered by group."""
    names = {f: text.word(PAGE, f"fates.{f}") for f in tab.FATES}
    labels = {c: text.word(PAGE, f"labels.{c}") for c in ("fate", "reads", "processed", *groups)}
    return BarData(
        frame=frame, columns=columns(frame), category="sample", value="percent",
        category_label=text.word(PAGE, "figures.fates.category"),
        value_label=text.word(PAGE, "figures.fates.value"), series="fate",
        colours={f: FATE_COLOURS[f] for f in tab.FATES}, caption=caption,
        horizontal=frame["sample"].n_unique() <= MAX_TICK_LABELS,
        stacked=True, hover=("reads", "processed"), names=names, filters=tuple(groups),
        labels=labels, summary_median=False,
    )  # fmt: skip


def detection_curves(frame: pl.DataFrame, cells: Sequence[str], common: float,
                     caption: str) -> CurveData:  # fmt: skip
    """F2: genes detected as each sample is thinned; the dashed line is the common depth."""
    top = tab.stat(frame["genes"], "max")
    return CurveData(
        frame=frame, columns=columns(frame), series="sample", x="depth", y="genes",
        x_label=text.word(PAGE, "figures.detection.x"),
        y_label=text.word(PAGE, "figures.detection.y"), colour="cell_type",
        colours=category_colours(list(cells)), hover=("condition", "batch"),
        thresholds=(Line(text.word(PAGE, "figures.detection.depth_line"), (common, common),
                         (0.0, top)),),
        caption=caption,
    )  # fmt: skip


def rho_levels(rho: np.ndarray) -> tuple[list[str | None], tuple[Swatch, ...]]:
    """Each cell's rho bin, in cividis from low to high; no bin on the diagonal."""
    off = rho[np.triu_indices(rho.shape[0], 1)]
    edges = np.unique(np.round(np.linspace(off.min(), off.max(), BINS + 1)[1:-1], 3))
    key = "figures.correlation.bins"
    pairs = zip(edges[:-1], edges[1:], strict=True)
    names = [text.fill(PAGE, f"{key}.below", edge=edges[0]),
             *(text.fill(PAGE, f"{key}.between", low=a, high=b) for a, b in pairs),
             text.fill(PAGE, f"{key}.above", edge=edges[-1])]  # fmt: skip
    stops = np.linspace(0, len(CIVIDIS) - 1, len(names)).round().astype(int)
    swatches = tuple(Swatch(name, CIVIDIS[i]) for name, i in zip(names, stops, strict=True))
    n = rho.shape[0]
    levels = [None if i == j else names[int(np.searchsorted(edges, rho[i, j], side="right"))]
              for i in range(n) for j in range(n)]  # fmt: skip
    return levels, swatches


def leaf_order(rho: np.ndarray, samples: Sequence[str]) -> tuple[str, ...]:
    """Samples in the leaf order of average linkage on 1 - rho, leaves optimally ordered."""
    tree = linkage(squareform(1 - rho, checks=False), "average", optimal_ordering=True)
    return tuple(samples[i] for i in leaves_list(tree))


def correlation_heatmap(rho: np.ndarray, sheet: pl.DataFrame,
                        caption: str) -> HeatmapData:  # fmt: skip
    """F3: the samples' Spearman rho, clustered, with cell type and batch strips on top.

    Only the lower triangle in drawn order is kept: rho is symmetric, and the full 210 by 210
    matrix drew an 18,579,509-byte figure page against the 15 MB cap (Step-1, 2026-10-08).
    """
    samples = sheet["sample"].to_list()
    order = leaf_order(rho, samples)
    levels, swatches = rho_levels(rho)
    position = {s: i for i, s in enumerate(order)}
    frame = pl.DataFrame({"sample_row": [s for s in samples for _ in samples],
                          "sample_column": samples * len(samples),
                          "rho": rho.ravel().round(4), "level": levels}).filter(
        pl.col("sample_row").replace_strict(position, return_dtype=pl.Int64)
        > pl.col("sample_column").replace_strict(position, return_dtype=pl.Int64))  # fmt: skip
    tracks = []
    for column in ("cell_type", "batch"):
        values = dict(zip(samples, sheet[column].cast(pl.String), strict=True))
        levels_seen = list(dict.fromkeys(values.values()))
        tracks.append(Track(text.word(PAGE, f"labels.{column}"), values,
                            category_colours(levels_seen)))  # fmt: skip
    axis = text.word(PAGE, "figures.correlation.axis")
    return HeatmapData(
        frame=frame, columns=columns(frame), row="sample_row", column="sample_column",
        value="rho", level="level", row_label=axis, column_label=axis, swatches=swatches,
        caption=caption, row_order=order, column_order=order, tracks=tuple(tracks),
    )  # fmt: skip


def length_classes(frame: pl.DataFrame) -> tuple[pl.DataFrame, dict[str, str]]:
    """The frame with each gene's effective-length class, and the class names in order."""
    bounds = {"short": SHORT_NT, "long": LONG_NT}
    names = {k: text.fill(PAGE, f"figures.units.{k}", **bounds) for k in ("short", "medium",
                                                                          "long")}  # fmt: skip
    length = pl.col("length")
    return frame.with_columns(
        length_class=pl.when(length < SHORT_NT).then(pl.lit(names["short"]))
        .when(length > LONG_NT).then(pl.lit(names["long"])).otherwise(pl.lit(names["medium"]))
    ), names  # fmt: skip


def units_scatter(frame: pl.DataFrame, names: Mapping[str, str], caption: str) -> ScatterData:
    """F4: one sample's genes, estimated reads against TPM, coloured by effective length."""
    colours = dict(zip(names.values(), (DATA_COLOURS[2], DATA_COLOURS[4], DATA_COLOURS[6]),
                       strict=True))  # fmt: skip
    return ScatterData(
        frame=frame, columns=columns(frame), x="tpm", y="counts",
        x_label=text.word(PAGE, "figures.units.x"), y_label=text.word(PAGE, "figures.units.y"),
        colour="length_class", colours=colours, hover=("gene_name", "length"), log_x=True,
        log_y=True, caption=caption,
    )  # fmt: skip


def type_names(frame: pl.DataFrame) -> pl.DataFrame:
    """Gene types as the page names them: the two it compares, every other as one word."""
    other = text.word(PAGE, "gene_types.other")
    named = {t: text.word(PAGE, f"gene_types.{t}") for t in tab.GENE_TYPES}
    return frame.with_columns(gene_type=pl.col("gene_type").replace_strict(
        named, default=pl.lit(other), return_dtype=pl.String))  # fmt: skip


def type_colours() -> dict[str, str]:
    """A colour per gene type the page names, in legend order."""
    levels = [*(text.word(PAGE, f"gene_types.{t}") for t in tab.GENE_TYPES),
              text.word(PAGE, "gene_types.other")]  # fmt: skip
    return dict(zip(levels, (DATA_COLOURS[4], DATA_COLOURS[0], DATA_COLOURS[6]), strict=True))


def estimates_scatter(frame: pl.DataFrame, caption: str) -> ScatterData:
    """F6: each gene's fold change from estimated counts against from length-scaled counts."""
    both = pl.concat([frame["estimated"], frame["length_scaled"]])
    low, high = tab.stat(both, "min"), tab.stat(both, "max")
    return ScatterData(
        frame=frame, columns=columns(frame), x="length_scaled", y="estimated",
        x_label=text.word(PAGE, "figures.estimates.x"),
        y_label=text.word(PAGE, "figures.estimates.y"), colour="gene_type",
        colours=type_colours(), hover=("gene_name", "difference"),
        lines=(Line(text.word(PAGE, "figures.estimates.identity"), (low, high), (low, high),
                    dash=True),),
        caption=caption,
    )  # fmt: skip


def isoform_bars(frame: pl.DataFrame, caption: str) -> BarData:
    """F7: median share of a gene's isoforms expressed, by isoforms annotated and gene type."""
    colours = {k: v for k, v in type_colours().items() if k in set(frame["gene_type"])}
    labels = {c: text.word(PAGE, f"labels.{c}") for c in ("gene_type", "n_genes")}
    return BarData(
        frame=frame, columns=columns(frame), category="annotated_bin", value="share",
        category_label=text.word(PAGE, "figures.isoforms.x"),
        value_label=text.word(PAGE, "figures.isoforms.y"), series="gene_type", colours=colours,
        caption=caption, hover=("n_genes",), filters=("gene_type",), labels=labels,
        summary_median=False,
    )  # fmt: skip


def test_caption(found: Caption) -> dict[str, Any]:
    """The placeholders a caption fills from ``ctx.viz.annotate``'s tests, one item per test.

    Each item names its feature, n per group, BH-adjusted p and rank-biserial r; the groups'
    names are said once, as ``reference`` and ``tested``, since a batch name runs long.
    """
    items = "; ".join(text.fill(PAGE, "tests_item", feature=t.feature, n_ref=t.n_reference,
                                n_test=t.n_comparison, p=t.p_adj, r=t.r)
                      for t in found.tests)  # fmt: skip
    return {"tests": items, "n_tests": len(found.tests), "reference": found.groups[0],
            "tested": found.groups[1]}  # fmt: skip
