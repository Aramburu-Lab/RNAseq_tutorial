"""M10's figures redrawn from the tables ``ctx clinical`` wrote, by the page's own builders.

The run is read back as ``ctx clinical`` held it when it drew: the patient table, the options
as resolved, and the stored values of the carried genes, the build list's among them. Each
drawer calls the builder that drew the figure on the page, with a writer that writes nothing
(``Recorder``), so the figure and its words are M10's own. The score, its scan with its seeded
permutations, and the per-gene screen the live panel lists are computed again by M10's code
from those values, as the run computed them.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

import numpy as np
import polars as pl

from ctx import runtime
from ctx.clinical import figure_tables, harmonise, page, panel
from ctx.clinical.run import Options
from ctx.dataset import Dataset
from ctx.figures.kinds_m10 import OneBox, OneEndpoint
from ctx.figures.recorder import Recorder
from ctx.report.builder import FigureText

FOLDER = "12_clinical"


def context(inputs: Mapping[str, Path], spec: Dataset) -> page.Context:
    """The page's context as ``ctx clinical`` built it, from its tables.

    Without ``panel_values`` the context holds no genes, for the figures that score none.
    The patient flow and Table 3's flags, the run's settings and its command feed only the
    page's tables and its Reproduce section, which no figure reads.
    """
    patients = pl.read_parquet(inputs["patients"])
    row = figure_tables.read_options(inputs["run_options"])
    options = Options(dataset_file=inputs["dataset_file"], endpoint=row["endpoint"],
                      split=row["split"], horizon=row["horizon"], covariates=row["covariates"],
                      strata=row["strata"])  # fmt: skip
    if "panel_values" in inputs:
        genes, values, chosen = figure_tables.read_values(inputs["panel_values"],
                                                          patients["sample"].to_list())  # fmt: skip
    else:
        genes = pl.DataFrame(schema={"gene_id": pl.String, "gene_name": pl.String})
        values, chosen = np.zeros((0, patients.height), np.int16), []
    carried = panel.Panel(np.arange(genes.height), panel.encode(values), genes.height,
                          patients.height)  # fmt: skip
    return page.Context(
        spec=spec, options=options, endpoint=row["endpoint"],
        clinical=harmonise.Clinical(patients, (), {}), values=values, gene_table=genes,
        carried=carried, chosen=chosen, source=row["source"], filtered=(0, 0),
        # No figure builder reads the settings; a None fails loudly if one ever does.
        settings=cast(runtime.Settings, None), command="", cards=(),
    )  # fmt: skip


def _recorder(spec: Dataset) -> Recorder:
    """A writer for M10's folder that keeps what the builders draw."""
    return Recorder(FOLDER, spec.slug)


def er(inputs: Mapping[str, Path], p: OneEndpoint, spec: Dataset) -> tuple[Any, FigureText]:
    """ER-negative against ER-positive under one endpoint."""
    recorder = _recorder(spec)
    page._er_check(recorder, context(inputs, spec))
    return recorder.one(f"er_{p.endpoint.lower()}")


def expression(inputs: Mapping[str, Path], p: OneBox, spec: Dataset) -> tuple[Any, FigureText]:
    """The build genes' log2 CPM by the declared two-level variable."""
    recorder = _recorder(spec)
    page._expression(recorder, context(inputs, spec))
    return recorder.one(f"expression_{p.box}")


def _score(inputs: Mapping[str, Path], spec: Dataset, what: str) -> tuple[Any, FigureText]:
    """The score's curves or its scan, drawn by ``page._scored`` with its live panel."""
    ctx = context(inputs, spec)
    result = page.scored(ctx)
    recorder = _recorder(spec)
    page._scored(recorder, ctx, result, page.screen_table(ctx, result))
    return recorder.one(what)


def km_score(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The build list's score split into groups, with the live survival panel."""
    return _score(inputs, spec, "km_score")


def scan_score(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """The cut-point scan of the build list's score, its permutation line."""
    return _score(inputs, spec, "scan_score")


def competing(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Death from the cancer and from other causes as competing risks."""
    recorder = _recorder(spec)
    page._competing(recorder, context(inputs, spec))
    return recorder.one("competing")


def power(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Power over events and hazard ratio, the endpoint's events marked."""
    ctx = context(inputs, spec)
    recorder = _recorder(spec)
    page._power(recorder, ctx, page.scored(ctx))
    return recorder.one("power")


def confounding(inputs: Mapping[str, Path], _: Any, spec: Dataset) -> tuple[Any, FigureText]:
    """Cramér's V between the patient variables."""
    recorder = _recorder(spec)
    page._context_groups(recorder, context(inputs, spec))
    return recorder.one("confounding")
