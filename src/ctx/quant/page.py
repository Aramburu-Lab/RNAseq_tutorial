"""The M4 page and its files: tables, figures, key numbers and the nine sections (D-045).

Every word is ``content/quant.toml``'s. Figures that need copy Q1 are drawn only when it is
there; the page says which are missing and why, and names the figures no run has made yet.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import polars as pl

from ctx import runtime, tools
from ctx.de.run import Writer
from ctx.quant import NUM_READS_DECIMALS, figure_tables
from ctx.quant import blocks as blk
from ctx.quant import tables as tab
from ctx.quant.figures import UNITS
from ctx.report import text
from ctx.report.builder import (
    FigureBlock,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
)
from ctx.viz.style import MAX_CATEGORIES

if TYPE_CHECKING:
    from ctx.quant.run import Inputs, Measured, Options

PAGE = MODULE = "quant"
# Sheet columns the bars may filter by, where the samples hold two levels or more.
GROUPS = ("cell_type", "condition", "batch")


@dataclass(frozen=True)
class Context:
    """Everything one run hands the page; ``made_up`` for the course preview's made-up data."""

    writer: Writer
    inputs: "Inputs"
    options: "Options"
    measured: "Measured"
    cards: tuple[ParameterCard, ...]
    command: str
    settings: runtime.Settings
    made_up: bool = False


def shown(sheet: pl.DataFrame) -> blk.Shown:
    """The sheet as the figures name it; a missing value is a word."""
    missing = text.word(PAGE, "no_value")
    kept = sheet.select("sample", *(pl.col(c).cast(pl.String).fill_null(missing)
                                    for c in (*GROUPS, "donor")))  # fmt: skip
    usable = [g for g in GROUPS if 1 < kept[g].n_unique() <= MAX_CATEGORIES]
    batches = tuple(kept["batch"].unique(maintain_order=True).to_list())
    return blk.Shown(kept, kept.select("sample", *usable), batches if len(batches) == 2 else ())  # noqa: PLR2004


def described(
    frame: pl.DataFrame, units: dict[str, str] | None = None
) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of a table, meanings from ``columns``."""
    units = units or {}
    return {c: (units.get(c, UNITS.get(c, "none")), text.word(PAGE, f"columns.{c}"))
            for c in frame.columns}  # fmt: skip


def library_table(library: pl.DataFrame) -> ResultTable:
    """One row per library type Salmon detected: samples, strand bias range, fragment lengths."""
    rows = []
    for (kind, source), part in library.group_by("library_type", "fragment_source",
                                                  maintain_order=True):  # fmt: skip
        bias = part["strand_bias"]
        span = f"{tab.stat(bias, 'min'):.3f} to {tab.stat(bias, 'max'):.3f}"
        rows.append((str(kind), f"{part.height:,}", span,
                     str(source), f"{tab.stat(part['fragment_mean'], 'median'):.1f}",
                     ", ".join(sorted(set(part["optimiser"].cast(pl.String))))))  # fmt: skip
    return ResultTable(text.word(PAGE, "tables.library.title"),
                       tuple(text.words(PAGE, "tables.library.columns")), tuple(rows),
                       numeric=(1, 4))  # fmt: skip


def notes(ctx: Context) -> list[Result]:
    """What the page says before any figure: the sums check, Q1's state, what no run made."""
    sums = ctx.measured.sums
    worst = tab.stat(sums["relative_difference"].abs(), "max")
    said: list[Result] = [text.fill(PAGE, "notes.sums", n=sums.height, worst=worst)]
    if ctx.made_up:
        return [text.word(PAGE, "preview.note"), *said]
    if (ref := ctx.measured.ref) is not None and ctx.measured.port is not None:
        said.append(
            text.fill(
                PAGE,
                "laptop.note",
                contigs=len(ref.genome.contigs),
                genes=len(ctx.measured.port.genes),
                atg=100 * (ref.start_codons.share or 0.0),
            )
        )
    elif ctx.inputs.salmon_dir is None:
        said.append(text.fill(PAGE, "notes.no_q1", dataset=ctx.inputs.spec.slug))
    else:
        said.append(text.fill(PAGE, "notes.q1", n=len(ctx.measured.port.samples)
                              if ctx.measured.port else 0))  # fmt: skip
    said.append(text.word(PAGE, "notes.not_yet"))
    return said


def salmon_group(ctx: Context, found: blk.Shown) -> Group:
    """F1, its tests and the library table; F3."""
    writer, measured = ctx.writer, ctx.measured
    fates_block, tests = blk.fates_figure(writer, measured.fates, found)
    items: list[FigureBlock | ResultTable | str] = [fates_block]
    if tests is not None:
        items.append(blk.test_table(tests, "fate_tests"))
        writer.table(blk.tests_frame(tests), "fate_tests", described(blk.tests_frame(tests)))
    items.append(library_table(tab.library(measured.salmon)))
    matrix = ctx.measured.counts.matrix.select(found.sheet["sample"].to_list()).to_numpy()
    items.append(blk.correlation_figure(writer, matrix.astype(float), found))
    return Group(text.word(PAGE, "groups.salmon"), tuple(items))


def q1_groups(ctx: Context, found: blk.Shown) -> list[Group]:
    """F2, F4, F6 and F7, from copy Q1 through the port; none without it."""
    port, folder = ctx.measured.port, ctx.measured.salmon_dir
    if port is None or folder is None:
        return []
    writer, genes = ctx.writer, ctx.inputs.genes
    names = dict(genes.select("gene_id", "gene_name").iter_rows())
    detection, curves = blk.detection_figure(writer, port, found)
    units, unit_frame = blk.units_figure(writer, port, ctx.measured.fates, names)
    writer.table(curves, "detection", described(curves))
    writer.table(unit_frame, "units", described(unit_frame))
    gene_items: list[FigureBlock | ResultTable | str] = [detection, units]
    if ctx.inputs.contrast is not None:
        used = blk.contrast_samples(ctx.inputs)
        block, _ = blk.estimates_figure(
            writer,
            port,
            used,
            genes,
            ctx.inputs.contrast.name,
            lambda f: writer.table(f, "fold_changes", described(f)),
        )
        gene_items.append(block)
    out = [Group(text.word(PAGE, "groups.genes"), tuple(gene_items))]
    isoforms = blk.isoforms_figure(writer, folder, port.samples, genes)
    if isoforms is not None:
        writer.table(isoforms[1], "isoforms", described(isoforms[1]))
        out.append(Group(text.word(PAGE, "groups.transcripts"), (isoforms[0],)))
    counts = port.frame(port.counts_from(ctx.options.counts_from), names)
    columns = {c: ("reads", text.fill(PAGE, "columns.sample_counts", sample=c,
                                      mode=str(ctx.options.counts_from)))
               for c in port.samples}  # fmt: skip
    writer.table(
        counts,
        "gene_counts",
        described(counts.select("gene_id", "gene_name")) | columns,
        csv_decimals=NUM_READS_DECIMALS,
    )
    return out


def key_numbers(ctx: Context) -> tuple[KeyNumber, ...]:
    """Samples, median mapped share, the commonest library type, median decoy share."""
    fates, library = ctx.measured.fates, tab.library(ctx.measured.salmon)
    share = fates.select(mapped=100 * pl.col("mapped") / pl.col("processed"),
                         decoy=100 * pl.col("decoy") / pl.col("processed"))  # fmt: skip
    kind, n_kind = library["library_type"].value_counts(sort=True).row(0)
    bias = library.filter(pl.col("library_type") == kind)["strand_bias"]
    return (
        KeyNumber(text.word(PAGE, "keys.samples"), f"{fates.height:,}"),
        KeyNumber(text.word(PAGE, "keys.mapping"), f"{tab.stat(share['mapped'], 'median'):.1f} %",
                  text.fill(PAGE, "keys.mapping_detail", low=tab.stat(share["mapped"], "min"),
                            high=tab.stat(share["mapped"], "max"))),
        KeyNumber(text.word(PAGE, "keys.library"), str(kind),
                  text.fill(PAGE, "keys.library_detail", n_kind=n_kind, n=fates.height,
                            low=tab.stat(bias, "min"), high=tab.stat(bias, "max"))),
        KeyNumber(text.word(PAGE, "keys.decoy"), f"{tab.stat(share['decoy'], 'median'):.1f} %",
                  text.fill(PAGE, "keys.decoy_detail", low=tab.stat(share["decoy"], "min"),
                            high=tab.stat(share["decoy"], "max"))),
    )  # fmt: skip


def write_tables(ctx: Context) -> None:
    """The per-sample tables: fates, library and the sums check."""
    measured = ctx.measured
    fates = measured.fates.with_columns(
        **{f"{f}_pct": 100 * pl.col(f) / pl.col("processed") for f in tab.FATES}
    )
    for what, frame in (("fates", fates), ("library", tab.library(measured.salmon)),
                        ("sums", measured.sums)):  # fmt: skip
        units = {f"{f}_pct": "%" for f in tab.FATES}
        ctx.writer.table(frame, what, described(frame, units))


def sections(ctx: Context) -> dict[str, Any]:
    """The fixed sections, Input filled from this run."""
    provenance, n = ctx.measured.counts.provenance, ctx.inputs.samples.height
    if ctx.made_up:
        filled = {"input": text.fill(PAGE, "preview.input", n=n)}
    elif ctx.measured.ref is not None:
        genome = ctx.options.genome.name if ctx.options.genome else ""
        filled = {
            "input": text.fill(
                PAGE, "laptop.input", n=n, genome=genome, salmon=provenance.get("salmon", "")
            ),
            "method": text.fill(
                PAGE,
                "laptop.method",
                decoys=text.word(
                    PAGE, f"laptop.decoys_{'on' if ctx.options.salmon.decoys else 'off'}"
                ),
            ),
        }
    else:
        filled = {"input": text.fill(PAGE, "input", n=n,
                                     release=provenance.get("nf-core/rnaseq", ""),
                                     salmon=provenance.get("salmon", ""))}  # fmt: skip
    return text.skeleton(PAGE) | filled


def build(ctx: Context) -> ModulePage:
    """Write every table and figure, and return the page."""
    found = shown(ctx.inputs.samples)
    write_tables(ctx)
    # What the figures draw that no table above holds, so ctx figure redraws them (F.1).
    figure_tables.write(ctx, found)
    salmon = salmon_group(ctx, found)
    body: Sequence[Result] = (*notes(ctx), salmon, *q1_groups(ctx, found))
    usage, budget = runtime.measure(ctx.settings.started), runtime.budget(MODULE)
    inputs = [ctx.inputs.sheet_path.name, ctx.inputs.gtf.name]
    inputs += [p.name for p in (ctx.inputs.rnaseq_dir, ctx.inputs.trim, ctx.options.genome,
                                ctx.measured.salmon_dir) if p is not None]  # fmt: skip
    return ModulePage(
        title=text.fill(PAGE, "title", dataset=ctx.inputs.spec.name), **sections(ctx),
        parameters=ctx.cards, results=tuple(body), key_numbers=key_numbers(ctx),
        reproduce=Reproduction(
            command=ctx.command, wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self, peak_rss_children=usage.peak_rss_children,
            profile=text.word(PAGE, "preview.profile") if ctx.made_up
            else f"{ctx.settings.profile}, {ctx.settings.threads} threads",
            budget=text.fill(PAGE, "preview.budget" if ctx.made_up else "budget",
                             minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions={**runtime.versions(), **ctx.measured.counts.provenance,
                      "multiqc": tools.TOOLS["multiqc"].version}, inputs=tuple(inputs),
        ),
    )  # fmt: skip
