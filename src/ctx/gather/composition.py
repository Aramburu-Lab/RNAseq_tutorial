"""M1's Figure 1, samples by cell type and condition, and the controls that recount it.

``recount`` is the one definition of the bars: runs, or distinct individuals, per cell type and
condition and per level of a colour column. The figure as built is ``recount`` over every run
coloured by batch; the page recounts the same rows for a contrast, a colour and a count through
its twin ``ctxRecount.counts`` (templates ``_qc.js.j2``), which ``tests/test_gather_recount.py``
holds to this one under node (work-order item 5, D-051, D-012).

A contrast's runs are the samples M5 compares (``ctx.de.deseq.select_samples``), every run of
each, so the individuals under a contrast are the n M5 tests with.
"""

from collections.abc import Sequence
from typing import Any

import polars as pl

from ctx.dataset import Dataset, DatasetError
from ctx.de.deseq import select_samples
from ctx.report import text
from ctx.viz.axes import axis_title
from ctx.viz.style import (
    MAX_CATEGORIES,
    PLOTLY_BAR_BASE_PX,
    PLOTLY_BAR_ROW_PX,
    PLOTLY_HEIGHT_PX,
    category_colours,
)

PAGE = "gather"
# The columns that make a bar; a colour column is any other curated column.
GROUPS = ("cell_type", "condition")
COUNTS = ("runs", "individuals")
# Curated columns that are never a colour: the bars' own, and the individual each run counts.
NOT_COLOURS = (*GROUPS, "donor", "replicate")


def runs(frame: pl.DataFrame, colours: Sequence[str]) -> pl.DataFrame:
    """One row per run in bar order: ``group``, its two parts, ``donor``, each colour as text.

    Bar order is cell type, then condition, as the figure has always drawn them; an empty
    batch reads "all", as the build's legend names a dataset without batches.
    """
    return (
        frame.with_columns(batch=pl.col("batch").replace("", "all"))
        .sort(*GROUPS, "run_accession", maintain_order=True)
        .select(*GROUPS, *(pl.col(c).cast(pl.String).fill_null("") for c in colours),
                group=pl.format("{} {}", *GROUPS), donor=pl.col("donor").cast(pl.String))
    )  # fmt: skip


def built(frame: pl.DataFrame) -> pl.DataFrame:
    """The figure as built: runs and individuals per cell type, condition and batch.

    ``recount`` over every run coloured by batch, in the table's order of old: by cell type,
    condition and batch.
    """
    rows = runs(frame, ["batch"])
    bars = recount(rows, "batch", "runs")
    parts = rows.unique("group", maintain_order=True).select("group", *GROUPS)
    return (
        bars.join(parts, on="group", maintain_order="left")
        .select(*GROUPS, batch="level", n_samples="n", n_individuals="other", group="group")
        .sort(*GROUPS, "batch", maintain_order=True)
    )


def recount(rows: pl.DataFrame, colour: str, count: str) -> pl.DataFrame:
    """The bars: ``n`` runs or distinct individuals per group and colour level.

    Groups in first-appearance order of ``rows``, levels within a group in theirs; a group's
    ``total`` is its runs, or its distinct individuals, which do not add up across levels
    when one individual gave runs in two of them. ``other`` is the count the bar does not
    show, which its hover names.
    """
    measure = {"runs": pl.len(), "individuals": pl.col("donor").n_unique()}
    if count not in measure:
        raise ValueError(f"count is one of {COUNTS}, not {count!r}")
    shown, other = measure[count], measure[COUNTS[1 - COUNTS.index(count)]]
    bars = rows.group_by("group", colour, maintain_order=True).agg(
        n=shown.cast(pl.Int64), other=other.cast(pl.Int64))  # fmt: skip
    totals = rows.group_by("group", maintain_order=True).agg(total=shown.cast(pl.Int64))
    return bars.rename({colour: "level"}).join(totals, on="group", maintain_order="left")


def smallest(bars: pl.DataFrame) -> tuple[int, str, int]:
    """Groups, and the smallest group's name and total: fewest first, a tie by name."""
    totals = bars.unique("group", maintain_order=True).sort("total", "group")
    return totals.height, str(totals["group"][0]), int(totals["total"][0])


def colours(frame: pl.DataFrame, spec: Dataset) -> list[str]:
    """The colour columns offered: batch, then the ones ``[dataset] colours`` names.

    A dataset file names them, because a column nested in the bars misleads as a colour:
    Step-1's tissue and time point nest in cell type and condition (figure-controls critique,
    line 232; Rytis 2026-10-18). A file naming none offers batch alone.

    Raises:
        DatasetError: If a named colour is not a curated column of the file, or has fewer than
            two or more than a page's categories of levels.
    """
    named = list(dict.fromkeys(spec.info.get("colours", ())))
    known = {c.name for c in spec.columns} - {*NOT_COLOURS, "batch"}
    problems = [f"[dataset] colours names {c}, which is not a curated column" for c in named
                if c not in known or c not in frame.columns]  # fmt: skip
    problems += [f"[dataset] colours names {c}, which holds {frame[c].n_unique()} levels; "
                 f"a colour needs 2 to {MAX_CATEGORIES}" for c in named
                 if c in known and c in frame.columns
                 and not 1 < frame[c].n_unique() <= MAX_CATEGORIES]  # fmt: skip
    if problems:
        raise DatasetError(f"{spec.path.name}: " + "; ".join(problems))
    return ["batch", *named]


def levels(rows: pl.DataFrame, colour: str, spec: Dataset) -> list[str]:
    """A colour's levels in legend order: batch by the dataset's bins, others by name."""
    present = set(rows[colour])
    if colour == "batch":
        return [b for b in (*spec.batch_bins, "all") if b in present]
    return sorted(present)


def members(frame: pl.DataFrame, spec: Dataset) -> list[tuple[Any, list[int]]]:
    """Each contrast and the bar-order indices of its runs: every run of a sample M5 compares."""
    ordered = frame.sort(*GROUPS, "run_accession", maintain_order=True)
    out = []
    for contrast in spec.contrasts:
        used, _ = select_samples(frame, contrast, spec.order)
        kept = set(used["sample"])
        out.append((contrast, [i for i, s in enumerate(ordered["sample"]) if s in kept]))
    return out


def _contrast_line(contrast: Any, n: int) -> str:
    """What a contrast's bars hold, and what no filter of them shows."""
    line = text.fill(PAGE, "figures.composition.live_contrast", runs=n,
                     contrast=contrast.label or contrast.name, variable=contrast.variable,
                     tested=contrast.tested, reference=contrast.reference)  # fmt: skip
    notes = [key for key, held in (("paired", contrast.paired),
                                   ("split", contrast.split is not None),
                                   ("donors", "donor" in contrast.filters)) if held]  # fmt: skip
    return " ".join([line, *(text.word(PAGE, f"figures.composition.notes.{k}") for k in notes)])


def controls(frame: pl.DataFrame, spec: Dataset) -> dict[str, Any]:
    """The recount panel's settings: the runs as columns, each contrast's runs, the words.

    ``colours`` holds each colour column's levels and their hex in legend order;
    ``contrasts`` each contrast's label, its runs' indices and its live line; ``words`` the
    selects' labels, the lead's templates and the hover's names.
    """
    offered = colours(frame, spec)
    rows = runs(frame, offered)

    def word(key: str) -> str:
        """The figure's words at ``key``."""
        return text.word(PAGE, f"figures.composition.{key}")

    contrasts = [{"label": c.label or c.name, "runs": kept, "line": _contrast_line(c, len(kept))}
                 for c, kept in members(frame, spec)]  # fmt: skip
    return {
        "kind": "recount",
        "rows": {c: rows[c].to_list() for c in ("group", "donor", *offered)},
        "colours": {c: [[level, hexa] for level, hexa in
                        category_colours(levels(rows, c, spec)).items()] for c in offered},
        "contrasts": contrasts,
        "all": text.fill(PAGE, "figures.composition.live_all", runs=rows.height),
        "words": {
            "contrast": word("contrast"), "every": word("every"), "colour": word("colour"),
            "count": word("count"), "counts": {c: word(f"counts.{c}") for c in COUNTS},
            "titles": {c: axis_title(word(f"counts.{c}").capitalize(), c) for c in COUNTS},
            "shown": text.raw(PAGE, "figures.composition.live_shown"),
            "lead": {f"{size}_{c}": text.raw(PAGE, f"figures.composition.lead_{size}_{c}")
                     for size in ("one", "many") for c in COUNTS},
            "category": word("category"), "empty": word("empty"), "noscript": word("noscript"),
        },
        "height": {"base": PLOTLY_BAR_BASE_PX, "row": PLOTLY_BAR_ROW_PX,
                   "least": PLOTLY_HEIGHT_PX["bar"]},
    }  # fmt: skip
