"""M7's two density figures with the views a reader looks up (P7.10, D-051, D-073).

F4 (one contrast's RNA against protein fold change, a density per group of where a gene was
called) is drawn at each fold change the call rule may take, each gene's
call at each layer made by ``classify`` again: the figure controls plan expects protein calls,
whose TMT ratios are compressed, to fall first as the threshold rises. F7 (Step-1's fold change
against copy numbers) can group its genes by where they were called instead of by class.
Nothing is recomputed in the page; every view is one DensityData drawn in Python.
"""

import base64
from collections.abc import Mapping
from dataclasses import replace
from typing import Any

import numpy as np
import polars as pl

from ctx.de.run import Writer
from ctx.de.stats import ALPHA, LFC_MIN, classify
from ctx.proteo import folds, measure
from ctx.proteo.figures import HDR, KIND_COLOURS, PAGE, class_colours, class_name, columns, lines
from ctx.report import text
from ctx.report.builder import FigureBlock, FigureText
from ctx.report.lookups import lookup
from ctx.viz.contour import DENSITY_MIN_POINTS
from ctx.viz.density import DensityData

# F4's fold changes, the ones the controls plan names for M5's lessons (D-051 amendment), at
# D-012's cutoff. A cutoff axis as well took CPTAC breast's page to 23,093,159 B.
FOLDS = (0.5, 1.0, 1.5, 2.0)
# The steps of a contour grid sent as 16-bit integers, its maximum the last.
QUANTA = 65535
CALLS_TABLE = {
    "lfc_min": ("log2 FC", "The smallest absolute log2 fold change a call needs, at both layers."),
    "alpha": ("dimensionless", "The BH-adjusted p a call needs to be below, at both layers."),
    **{f"{layer}_{d}": ("genes", f"Genes called {d} at the {layer} layer.")
       for layer in ("rna", "protein") for d in ("up", "down")},
    **{k: ("genes", f"Genes of kind {k}: where they were called.") for k in folds.KINDS},
    "opposite": ("genes", "Genes called at both layers in opposite directions."),
}  # fmt: skip


def kinds(frame: pl.DataFrame, alpha: float, lfc_min: float) -> pl.DataFrame:
    """Each gene's calls at both layers at one setting, by the one rule, and its kind."""
    return frame.with_columns(
        rna_call=classify(pl.col("rna_padj"), pl.col("rna_lfc"), alpha, lfc_min),
        protein_call=classify(pl.col("protein_padj"), pl.col("protein_lfc"), alpha, lfc_min),
    ).with_columns(kind_code=folds.kind(pl.col("rna_call"), pl.col("protein_call")))


def counted(frame: pl.DataFrame) -> dict[str, int]:
    """Calls per layer and direction, per kind, and both-layer calls in opposite ways."""
    out = {f"{layer}_{d}": int((frame[f"{layer}_call"] == d).sum())
           for layer in ("rna", "protein") for d in ("up", "down")}  # fmt: skip
    out |= {k: int((frame["kind_code"] == k).sum()) for k in folds.KINDS}
    out["opposite"] = int(((frame["kind_code"] == folds.KINDS[0])
                           & (frame["rna_call"] != frame["protein_call"])).sum())  # fmt: skip
    return out


def _contrast_data(genes: pl.DataFrame, values: Mapping[str, object], alpha: float,
                   lfc_min: float) -> tuple[DensityData, FigureText]:  # fmt: skip
    """F4 at one setting of the call rule: its contract and its words."""
    names = {k: text.word(PAGE, f"figures.contrast.kinds.{k}") for k in folds.KINDS}
    called = kinds(genes, alpha, lfc_min)
    frame = called.select(
        "gene", "rna_lfc", "protein_lfc", "rna_padj", "protein_padj", "rna_call", "protein_call",
        called=pl.col("kind_code") == folds.KINDS[0],
        kind=pl.col("kind_code").replace_strict(names),
    )  # fmt: skip
    present = set(frame["kind"])
    colours = {names[k]: c for k, c in KIND_COLOURS.items() if names[k] in present}
    shown = {**values, **counted(called), "alpha": alpha, "lfc": lfc_min}
    data = DensityData(
        frame=frame, columns=columns(*frame.columns), x="rna_lfc", y="protein_lfc", group="kind",
        x_label=text.word(PAGE, "figures.contrast.x"),
        y_label=text.word(PAGE, "figures.contrast.y"), colours=colours, caption="",
        highlight="called", highlight_label=names[folds.KINDS[0]],
        hover=("gene", "rna_padj", "protein_padj"),
    )  # fmt: skip
    data = replace(data, caption=text.fill(PAGE, "figures.contrast.caption", **shown, **HDR,
                                           bandwidth=data.bandwidth()))  # fmt: skip
    words = FigureText(
        title=text.fill(PAGE, "figures.contrast.title", **shown),
        lead=text.fill(PAGE, "figures.contrast.lead", **shown),
        reading=lines("figures.contrast.reading", minimum=DENSITY_MIN_POINTS, **HDR,
                      **{**values, "alpha": ALPHA, "lfc": LFC_MIN}),
    )  # fmt: skip
    return data, words


def contrast_density(writer: Writer, made: folds.Folds, values: Mapping[str, object]
                     ) -> FigureBlock:  # fmt: skip
    """F4: one contrast's RNA against protein fold change on the same tumours, every setting.

    A density per group of where a gene was called; the genes called at both layers are drawn
    as points over the contours. The build's setting, D-012's, is the first view.
    """
    options, rows = [], []
    for lfc_min in sorted(FOLDS, key=lambda f: f != LFC_MIN):
        data, words = _contrast_data(made.genes, values, ALPHA, lfc_min)
        options.append((text.fill(PAGE, "figures.contrast.option", alpha=ALPHA, lfc=lfc_min),
                        data, words))  # fmt: skip
        rows.append({"lfc_min": lfc_min, "alpha": ALPHA,
                     **counted(kinds(made.genes, ALPHA, lfc_min))})  # fmt: skip
    writer.table(pl.DataFrame(rows), f"contrast_{made.contrast.name}_calls", CALLS_TABLE)
    block = writer.figure(options[0][1], f"contrast_{made.contrast.name}", options[0][2])
    return packed(lookup(block, options, text.word(PAGE, "figures.contrast.view")))


def _grid(z: Any) -> np.ndarray:
    """A contour's grid as plotly.py hands it over: a typed array in base64, or nested lists."""
    if isinstance(z, dict):
        shape = tuple(int(n) for n in str(z["shape"]).split(","))
        flat = np.frombuffer(base64.b64decode(z["bdata"]), dtype=np.dtype(z["dtype"]))
        return flat.reshape(shape).astype(float)
    return np.asarray(z, dtype=float)


def _typed(trace: dict[str, Any]) -> dict[str, Any]:
    """A contour trace with its grid in 16 bits of its maximum, and its level on that scale.

    The contour is where the grid crosses the level, so scaling both by one number moves
    nothing; 16 bits put each value within 1/131,070 of the grid's maximum of where float64
    has it.
    """
    if trace.get("type") != "contour":
        return trace
    z = _grid(trace["z"])
    top = float(z.max()) or 1.0
    q = np.rint(z / top * QUANTA).astype("<u2")
    contours = {**trace["contours"], "start": trace["contours"]["start"] / top * QUANTA,
                "end": trace["contours"]["end"] / top * QUANTA}  # fmt: skip
    grid = {"dtype": "u2", "shape": ",".join(map(str, q.shape)),
            "bdata": base64.b64encode(q.tobytes()).decode("ascii")}  # fmt: skip
    return {**trace, "z": grid, "contours": contours}


def packed(block: FigureBlock) -> FigureBlock:
    """A density view's options with each contour grid in 16 bits, in base64.

    A contour carries its whole grid: a view of a contrast's four groups took about 470 kB
    with its grids as float32, and CPTAC breast's page 15,981,038 B with four fold changes
    for each of its two contrasts. The figure the page draws at load keeps its float grids;
    the build's view, drawn again after another, takes the 16-bit one.
    """
    assert block.controls is not None  # noqa: S101 - lookup sets them
    controls = dict(block.controls)
    controls["options"] = [
        {**o, "figure": {**o["figure"], "data": [_typed(t) for t in o["figure"]["data"]]}}
        for o in controls["options"]
    ]  # fmt: skip
    return replace(block, controls=controls)


def _fold_data(tested: pl.DataFrame, by_kind: bool) -> DensityData:
    """F7's contract: a density per gene class, the declared classes marked, or per kind."""
    frame = tested.select(
        "gene", "gene_class", "rna_lfc", "protein_lfc", "rna_padj", "protein_padj"
    ).with_columns(keep=pl.col("gene_class") != measure.OTHER)
    frame = frame.with_columns(gene_class=pl.col("gene_class").replace_strict(
        {c: class_name(c) for c in frame["gene_class"].unique()}))  # fmt: skip
    if not by_kind:
        return DensityData(
            frame=frame, columns=columns(*frame.columns), x="rna_lfc", y="protein_lfc",
            group="gene_class", x_label=text.word(PAGE, "figures.fold.x"),
            y_label=text.word(PAGE, "figures.fold.y"), colours=class_colours(frame), caption="",
            highlight="keep", highlight_label=", ".join(
                c for c in class_colours(frame) if c != class_name(measure.OTHER)),
            hover=("gene", "rna_padj", "protein_padj"),
        )  # fmt: skip
    names = {k: text.word(PAGE, f"figures.contrast.kinds.{k}") for k in folds.KINDS}
    called = kinds(frame, ALPHA, LFC_MIN).with_columns(
        kind=pl.col("kind_code").replace_strict(names)).drop("kind_code")  # fmt: skip
    present = set(called["kind"])
    return DensityData(
        frame=called, columns=columns(*called.columns), x="rna_lfc", y="protein_lfc",
        group="kind", x_label=text.word(PAGE, "figures.fold.x"),
        y_label=text.word(PAGE, "figures.fold.y"),
        colours={names[k]: c for k, c in KIND_COLOURS.items() if names[k] in present},
        caption="", highlight="keep", highlight_label=", ".join(
            c for c in class_colours(frame) if c != class_name(measure.OTHER)),
        hover=("gene", "rna_padj", "protein_padj"),
    )  # fmt: skip


def fold(writer: Writer, tested: pl.DataFrame, values: Mapping[str, object]) -> FigureBlock:
    """F7: RNA against protein fold change per gene, a density per class, the class marked.

    A view draws a density per kind of call instead, at D-012's rule, the class still marked.
    """
    words = FigureText(
        title=text.fill(PAGE, "figures.fold.title", **values),
        lead=text.fill(PAGE, "figures.fold.lead", **values),
        reading=lines("figures.fold.reading", minimum=DENSITY_MIN_POINTS, **HDR))  # fmt: skip
    options = []
    for by_kind in (False, True):
        data = _fold_data(tested, by_kind)
        data = replace(data, caption=text.fill(PAGE, "figures.fold.caption", **values, **HDR,
                                               bandwidth=data.bandwidth()))  # fmt: skip
        options.append((text.word(PAGE, f"figures.fold.view_{int(by_kind)}"), data, words))
    block = writer.figure(options[0][1], "fold_change", words)
    return packed(lookup(block, options, text.word(PAGE, "figures.fold.view")))
