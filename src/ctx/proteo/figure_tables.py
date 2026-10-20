"""The tables M7 writes so ``ctx figure`` can redraw its figures from tables alone (F.1).

Each holds what a figure drew and no other table holds: the run's options; which protein cells
hold a value, and each protein sample's groups; the RNA abundance the coverage tertiles are cut
on; the RNA and protein of the genes the worked and paired figures draw, beside the paired
samples' annotation; each contrast's genes tested at both layers and its correlation. They are
written before any figure, so a figure's spec finds them, and read back here as the
``Measured`` the page's block functions take.

The protein layer is written as presence alone, 1 or 0 per cell: F1 reads nothing else of it,
and its values took 3,526,528 B as parquet and 8,676,633 B as CSV on CPTAC breast. The RNA and
protein matrices are written for the drawn genes only: the full pair took 8,221,928 B as
parquet there, and no figure reads another row.
"""

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import polars as pl

from ctx.dataset import Dataset, ProteoSpec
from ctx.de.run import Writer
from ctx.figures.kinds_m7 import TABLES, folds_what, summary_what
from ctx.proteo import figures
from ctx.proteo.folds import Folds
from ctx.proteo.run import Measured, Options
from ctx.proteo.tables import Layer
from ctx.report import text

PAGE = "proteo"
# The drawn rows of the correlation matrices, and of the paired design's both-tissue ones.
PAIRED, WORKED = "paired", "worked"
ROW_SCHEMA = {"gene": pl.String, "view": pl.String, "column": pl.Int64, "rna": pl.Float64,
              "protein": pl.Float64}  # fmt: skip
SUMMARY_SCHEMA = {"n_reference": pl.Int64, "n_tested": pl.Int64, "fold_rho": pl.Float64,
                  "fold_rho_low": pl.Float64, "fold_rho_high": pl.Float64,
                  "paired": pl.Boolean}  # fmt: skip


def drawn_genes(measured: Measured, p: ProteoSpec, top: int) -> list[str]:
    """The genes F5 and F6 draw, as ``page._across`` chooses them.

    The worked genes that were tested, then F5's: the named genes and the ``top`` highest and
    lowest rho (``figures.heatmap_genes``).
    """
    tested = set(measured.genes.filter("tested")["gene"])
    named = [*(k for k in p.published.get("mertins", {}) if k.isupper()), *p.genes]
    chosen = figures.heatmap_genes(measured.genes, named, top)
    return list(dict.fromkeys([*(g for g in p.genes if g in tested), *chosen]))


def _rows(genes: list[str], order: dict[str, int], layers: tuple[np.ndarray, np.ndarray],
          view: str) -> pl.DataFrame:  # fmt: skip
    """One row per drawn gene and matrix column: its RNA and its protein, NaN kept."""
    rna, protein = layers
    at, n = [order[g] for g in genes], rna.shape[1]
    return pl.DataFrame({"gene": np.repeat(np.array(genes, dtype=object), n).tolist(),
                         "view": [view] * (len(genes) * n),
                         "column": np.tile(np.arange(n, dtype=np.int64), len(genes)),
                         "rna": rna[at].ravel(), "protein": protein[at].ravel()},
                        schema=ROW_SCHEMA)  # fmt: skip


def _correlation(measured: Measured, p: ProteoSpec, top: int) -> list[tuple[str, pl.DataFrame]]:
    """The paired samples' annotation and the drawn genes' rows, across samples only."""
    extra = measured.extra
    order = {g: i for i, g in enumerate(measured.genes["gene"])}
    genes = drawn_genes(measured, p, top)
    rows = [_rows(genes, order, (extra["rna"], extra["values"]), PAIRED)]
    found = [(TABLES["paired_annotation"], extra["annotation"])]
    if WORKED in extra:
        rna, protein, shown = extra[WORKED]
        rows.append(_rows(genes, order, (rna, protein), WORKED))
        found.append((TABLES["worked_samples"], shown))
    return [*found, (TABLES["drawn_rows"], pl.concat(rows))]


def _folds(made: Folds) -> list[tuple[str, pl.DataFrame]]:
    """One contrast's genes tested at both layers, and its sample counts and correlation."""
    summary = pl.DataFrame({"n_reference": [made.n_protein[0]], "n_tested": [made.n_protein[1]],
                            "fold_rho": [made.rho[0]], "fold_rho_low": [made.rho[1]],
                            "fold_rho_high": [made.rho[2]], "paired": [made.paired]},
                           schema=SUMMARY_SCHEMA)  # fmt: skip
    name = made.contrast.name
    return [(folds_what(name), made.genes), (summary_what(name), summary)]


def tables(measured: Measured, spec: Dataset, options: Options, top: int
           ) -> list[tuple[str, pl.DataFrame]]:  # fmt: skip
    """Every table a figure of this run reads, by its output description."""
    p, layer = spec.proteo, measured.protein
    assert p is not None  # noqa: S101 - run() checks
    present = ~np.isnan(layer.values)
    found = [
        (TABLES["proteo_options"], pl.DataFrame({
            "method": [options.method], "min_quantified": [float(options.min_quantified)],
            "lfc_min": [float(options.lfc_min)], "shown_class": [options.gene_class]})),
        (TABLES["protein_samples"], pl.DataFrame({
            "sample": list(layer.samples), "colour": measured.extra.get("colour", measured.groups),
            "group": measured.groups}, schema=dict.fromkeys(("sample", "colour", "group"),
                                                           pl.String))),
        (TABLES["protein_present"], pl.DataFrame({"gene": list(layer.genes),
                                          **{s: present[:, j].astype(np.uint8)
                                             for j, s in enumerate(layer.samples)}})),
        (TABLES["abundance"], pl.DataFrame({"gene": measured.genes["gene"],
                                    "abundance": np.asarray(measured.abundance, dtype=float)})),
    ]  # fmt: skip
    if "rna" in measured.extra:
        found += _correlation(measured, p, top)
    for made in measured.extra.get("folds", []):
        found += _folds(made)
    return found


def _columns(what: str, frame: pl.DataFrame) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of each column: the content file's, a sample's or an annotation's."""
    if what == TABLES["protein_present"]:
        return {"gene": figures.columns("gene")["gene"],
                **{s: ("none", text.fill(PAGE, "figure_tables.present", sample=s))
                   for s in frame.columns[1:]}}  # fmt: skip
    if what == TABLES["paired_annotation"]:
        return {c: ("none", text.fill(PAGE, "figure_tables.annotation", column=c))
                for c in frame.columns}  # fmt: skip
    return figures.columns(*frame.columns)


def write(writer: Writer, measured: Measured, spec: Dataset, options: Options, top: int) -> None:
    """Write every table a figure reads, each with its README, before any figure."""
    for what, frame in tables(measured, spec, options, top):
        writer.table(frame, what, _columns(what, frame))


def options_from(path: Path, dataset_file: Path) -> Options:
    """The run's options, as ``ctx proteo`` was given them."""
    row = pl.read_parquet(path).row(0, named=True)
    return Options(dataset_file, row["method"], row["min_quantified"], row["lfc_min"],
                   row["shown_class"])  # fmt: skip


def _matrices(rows: pl.DataFrame, view: str, order: dict[str, int], n: int
              ) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """RNA and protein, genes by columns, the drawn rows filled and every other row NaN."""
    shown = rows.filter(pl.col("view") == view)
    at = np.array([order[g] for g in shown["gene"]], dtype=np.int64)
    column = shown["column"].to_numpy()
    out = []
    for layer in ("rna", "protein"):
        matrix = np.full((len(order), n), np.nan)
        matrix[at, column] = shown[layer].to_numpy()
        out.append(matrix)
    return out[0], out[1]


def measured_from(paths: Mapping[str, Path]) -> Measured:
    """What the page's block functions read of a run, from the tables a figure's inputs name.

    The genes always; the protein layer as F1 reads it, its samples' groups and the abundance
    F2 cuts on where ``protein_present`` is given; the drawn genes' matrices where
    ``drawn_rows`` is. A present cell reads back as 0.0 and a missing one as NaN: F1 counts
    present cells only. A part not given is left empty, as no figure of the kind reads it.
    """
    genes = pl.read_parquet(paths["proteo_genes"])
    measured = Measured(Layer((), (), np.empty((0, 0))), [], genes, np.empty(0))
    if "protein_samples" in paths:
        samples = pl.read_parquet(paths["protein_samples"])
        measured.groups = samples["group"].to_list()
        measured.extra["colour"] = samples["colour"].to_list()
        if "protein_present" in paths:
            present = pl.read_parquet(paths["protein_present"])
            names = tuple(samples["sample"])
            values = np.where(present.select(names).to_numpy() == 1, 0.0, np.nan)
            measured.protein = Layer(tuple(present["gene"]), names, values)
            measured.abundance = pl.read_parquet(paths["abundance"])["abundance"].to_numpy()
    if "drawn_rows" in paths:
        _correlation_from(paths, measured)
    return measured


def _correlation_from(paths: Mapping[str, Path], measured: Measured) -> None:
    """The paired samples' annotation and the drawn genes' matrices, into ``measured``."""
    annotation = pl.read_parquet(paths["paired_annotation"])
    rows = pl.read_parquet(paths["drawn_rows"])
    order = {g: i for i, g in enumerate(measured.genes["gene"])}
    rna, values = _matrices(rows, PAIRED, order, annotation.height)
    # The annotation holds one row per paired sample, as ``n_samples`` counts them.
    measured.extra |= {"rna": rna, "values": values, "annotation": annotation,
                       "n_samples": annotation.height}  # fmt: skip
    if "worked_samples" in paths:
        shown = pl.read_parquet(paths["worked_samples"])
        measured.extra[WORKED] = (*_matrices(rows, WORKED, order, shown.height), shown)
        # The tumours then the normals, each named by tissue: the normals' class figure's n.
        measured.extra["n_normal"] = int((shown["tissue"] == "normal").sum())


def folds_from(folds: Path, summary: Path, spec: Dataset, contrast: str) -> Folds:
    """One contrast's F4 as ``ctx proteo`` made it, its genes and numbers read back."""
    row = pl.read_parquet(summary).row(0, named=True)
    declared = next(c for c in spec.contrasts if c.name == contrast)
    return Folds(declared, pl.read_parquet(folds), (row["n_reference"], row["n_tested"]),
                 (row["fold_rho"], row["fold_rho_low"], row["fold_rho_high"]),
                 paired=row["paired"])  # fmt: skip
