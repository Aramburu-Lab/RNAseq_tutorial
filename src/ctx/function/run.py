"""M6: ``ctx function``. Enrichment of every contrast's results, then tables, figures, page.

Reads the newest M5 results of each declared contrast and the gene sets ``ctx setup``
fetched. GSEA ranks by the Wald statistic, one gene per symbol; ORA counts the up and the
down calls separately against the genes DESeq2 tested.
"""

import hashlib
import re
from collections.abc import Mapping
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl

from ctx import runtime
from ctx.dataset import DatasetError, load
from ctx.de.run import Writer, newest
from ctx.de.stats import ALPHA
from ctx.function import lite
from ctx.io.naming import MODULES, stem
from ctx.report.builder import ModulePage, ParameterCard, Reproduction, SelfCheck, write
from ctx.viz.bar import BarData
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import ROLES, category_colours

MODULE = "function"
TOP_SETS = 15
GSEA_COLUMNS = {
    "gene_set": ("none", "Gene set name."),
    "size": ("genes", "Set members among the ranked genes."),
    "es": ("dimensionless", "Enrichment score: the running sum's extreme."),
    "nes": ("dimensionless", "ES over the mean null ES of the same sign."),
    "p": ("dimensionless", "Nominal p from gene-set permutations; never 0."),
    "fdr": ("dimensionless", "FDR from the pooled normalised nulls (GSEA)."),
}
ORA_COLUMNS = {
    "direction": ("none", "up or down: which calls were counted."),
    "gene_set": ("none", "Gene set name."),
    "size": ("genes", "Set members among the tested genes."),
    "hits": ("genes", "Calls in the set."),
    "expected": ("genes", "Calls expected in the set by chance."),
    "p": ("dimensionless", "Hypergeometric P(X >= hits)."),
    "fdr": ("dimensionless", "BH over the sets tested."),
    "genes": ("none", "The calls in the set; ; separated."),
}
QUESTION = "Which pathways move together, beyond single genes?"
CONCEPT = (
    "Over-representation asks whether a pathway holds more of the called genes than chance, "
    "counted against the genes that could have been called. GSEA needs no cutoff: it walks "
    "down all genes ranked by the Wald statistic and asks whether a pathway's genes gather "
    "at one end."
)
PRACTICE = (
    "Use the tested genes as the background, rank by a signed statistic, and read the NES "
    "with its FDR. A pathway can move with no single gene significant, and one strong gene "
    "does not make a pathway."
)
THINK = (
    SelfCheck(
        "Why does the gene-set union background make the same hits look more significant?",
        "It adds genes that could never have been called, such as ones not expressed in the "
        "cells, so the expected overlap shrinks and ordinary overlaps look rare.",
    ),
    SelfCheck(
        "Why rank by the Wald statistic rather than by the fold change?",
        "The statistic weighs the change by its certainty; a large fold change on few "
        "counts would otherwise sit at the top.",
    ),
)
MISTAKES = (
    "Using the whole genome, or every gene in the gene sets, as the ORA background.",
    "Reading a GSEA p-value of 1/1001 as exact: it is limited by the permutations.",
    "Mixing gene IDs and symbols: duplicated symbols must be collapsed first.",
)


def _ranking(results: pl.DataFrame) -> tuple[list[str], np.ndarray]:
    """Symbols and Wald statistics, one gene per symbol (highest baseMean)."""
    rows = results.filter(pl.col("stat").is_not_null())
    return lite.collapse(
        rows["gene_id"].to_list(),
        rows["gene_name"].fill_null("").to_list(),
        rows["stat"].to_list(),
        rows["baseMean"].to_list(),
    )


def _symbols(results: pl.DataFrame, call: str | None = None) -> list[str]:
    """Symbols of the tested genes, or of the genes with one call."""
    tested = set(lite.background(results["gene_id"].to_list(), results["padj"].to_list()))
    rows = results.filter(pl.col("gene_id").is_in(list(tested)))
    if call:
        rows = rows.filter(pl.col("call") == call)
    return sorted({s for s in rows["gene_name"].drop_nulls() if s})


def _gsea(results: pl.DataFrame, sets: Mapping[str, tuple[str, ...]]) -> list[lite.Enriched]:
    """Preranked GSEA of one contrast."""
    names, scores = _ranking(results)
    return lite.prerank(names, scores, sets)


def _ora(
    results: pl.DataFrame, sets: Mapping[str, tuple[str, ...]], background: list[str]
) -> pl.DataFrame:
    """ORA of the up and the down calls, against the given background."""
    rows = []
    for direction in ("up", "down"):
        for r in lite.ora(_symbols(results, direction), background, sets):
            rows.append(
                (direction, r.name, r.size, r.hits, r.expected, r.p, r.fdr, ";".join(r.genes))
            )
    return pl.DataFrame(rows, schema=list(ORA_COLUMNS), orient="row")


def _gsea_bars(name: str, frame: pl.DataFrame) -> BarData:
    """NES of the sets with the smallest FDR."""
    top = (
        frame.sort("fdr", "gene_set")
        .head(TOP_SETS)
        .with_columns(
            label=pl.col("gene_set")
            .str.replace("^HALLMARK_", "")
            .str.replace_all("_", " ")
            .str.to_lowercase(),
            significance=pl.when(pl.col("fdr") < ALPHA)
            .then(pl.lit(f"FDR < {ALPHA}"))
            .otherwise(pl.lit("not significant")),
        )
    )
    columns = {
        **GSEA_COLUMNS,
        "label": ("none", "Set name, shortened."),
        "significance": ("none", f"FDR below {ALPHA} or not."),
    }
    return BarData(
        frame=top,
        columns=columns,
        category="label",
        value="nes",
        category_label="Gene set",
        value_label="NES",
        series="significance",
        colours={f"FDR < {ALPHA}": ROLES["comparison"], "not significant": None},
        caption=(
            f"{name}: preranked GSEA on the Wald statistic, {lite.PERMUTATIONS} gene-set "
            f"permutations, seed {lite.SEED}; the {TOP_SETS} sets with the smallest FDR."
        ),
        horizontal=True,
        hover=("fdr", "size"),
        width="double",
    )


def _ora_bars(name: str, frame: pl.DataFrame, background: int) -> BarData:
    """-log10 FDR of the top ORA sets per direction."""
    top = (
        frame.sort("fdr", "gene_set")
        .group_by("direction", maintain_order=True)
        .head(TOP_SETS // 2)
        .with_columns(
            neg_log10_fdr=-pl.col("fdr").log10(),
            label=pl.col("gene_set")
            .str.replace("^HALLMARK_", "")
            .str.replace_all("_", " ")
            .str.to_lowercase(),
        )
        .unique(["label", "direction"], maintain_order=True)
    )
    columns = {
        **ORA_COLUMNS,
        "neg_log10_fdr": ("-log10 FDR", "Minus log10 of the FDR."),
        "label": ("none", "Set name, shortened."),
    }
    return BarData(
        frame=top,
        columns=columns,
        category="label",
        value="neg_log10_fdr",
        category_label="Gene set",
        value_label="Significance",
        series="direction",
        colours={"up": ROLES["comparison"], "down": ROLES["reference"]},
        caption=(
            f"{name}: hypergeometric ORA of the up and the down calls (padj < {ALPHA}), "
            f"BH per direction, against {background:,} tested genes; the dashed line is "
            f"FDR = {ALPHA}."
        ),
        horizontal=True,
        hline=(-float(np.log10(ALPHA)), f"FDR {ALPHA}"),
        hover=("hits", "size"),
        width="double",
    )


def _featured_bars(featured: str, per_contrast: Mapping[str, pl.DataFrame]) -> BarData:
    """One set's NES in every contrast."""
    rows = []
    for name, frame in per_contrast.items():
        found = frame.filter(pl.col("gene_set") == featured)
        if found.height:
            nes, fdr = found.row(0, named=True)["nes"], found.row(0, named=True)["fdr"]
            rows.append((name, nes, fdr, f"FDR < {ALPHA}" if fdr < ALPHA else "not significant"))
    frame = pl.DataFrame(rows, schema=["contrast", "nes", "fdr", "significance"], orient="row")
    columns = {
        "contrast": ("none", "Contrast."),
        "nes": GSEA_COLUMNS["nes"],
        "fdr": GSEA_COLUMNS["fdr"],
        "significance": ("none", f"FDR below {ALPHA} or not."),
    }
    return BarData(
        frame=frame,
        columns=columns,
        category="contrast",
        value="nes",
        category_label="Contrast",
        value_label="NES",
        series="significance",
        colours={f"FDR < {ALPHA}": ROLES["comparison"], "not significant": None},
        caption=f"{featured}: NES in every contrast, each ranked and permuted separately.",
        horizontal=True,
        hover=("fdr",),
        width="double",
    )


def _running(name: str, result: lite.Enriched) -> ScatterData:
    """The running enrichment score, with the set's genes as ticks along the ranking."""
    ranks = tuple(float(i) for i in range(len(result.running)))
    frame = pl.DataFrame(
        {
            "rank": [float(p) for p in result.positions],
            "running": [result.running[p] for p in result.positions],
            "mark": ["set member"] * result.size,
        }
    )
    columns = {
        "rank": ("rank", "Position in the ranking, from the most positive statistic."),
        "running": ("dimensionless", "Running enrichment score at this member."),
        "mark": ("none", "Set member."),
    }
    caption = (
        f"{result.name} in {name}: ES = {result.es:.2f}, NES = {result.nes:.2f}, nominal "
        f"p = {result.p:.2g}, FDR = {result.fdr:.2g}; {result.size} members of "
        f"{len(ranks):,} ranked genes."
    )
    return ScatterData(
        frame=frame,
        columns=columns,
        x="rank",
        y="running",
        x_label="Rank",
        y_label="Enrichment score",
        colour="mark",
        colours={"set member": ROLES["comparison"]},
        caption=caption,
        lines=(Line("running score", ranks, result.running),),
        width="double",
    )


def _trap(name: str, results: pl.DataFrame, sets: Mapping[str, tuple[str, ...]]) -> BarData:
    """The same up calls tested against the tested genes and against the gene-set union."""
    tested = _symbols(results)
    union = sorted({g for genes in sets.values() for g in genes} | set(tested))
    hits = _symbols(results, "up")
    right = {r.name: r for r in lite.ora(hits, tested, sets)}
    wrong = {r.name: r for r in lite.ora(hits, union, sets)}
    shown = sorted(right, key=lambda s: right[s].p)[: TOP_SETS // 2]
    rows = [
        (s.replace("HALLMARK_", "").replace("_", " ").lower(), label, float(-np.log10(r[s].p)))
        for s in shown
        for label, r in (("tested genes", right), ("gene-set union", wrong))
    ]
    frame = pl.DataFrame(rows, schema=["gene_set", "background", "neg_log10_p"], orient="row")
    columns = {
        "gene_set": ("none", "Gene set, shortened."),
        "background": ("none", "ORA background."),
        "neg_log10_p": ("-log10 p", "Hypergeometric p of the up calls."),
    }
    return BarData(
        frame=frame,
        columns=columns,
        category="gene_set",
        value="neg_log10_p",
        category_label="Gene set",
        value_label="Significance",
        series="background",
        colours=category_colours(["tested genes", "gene-set union"]),
        caption=(
            f"{name}: the same {len(hits)} up calls against {len(tested):,} tested genes "
            f"and against {len(union):,} genes, every gene in the sets plus the tested "
            "ones (gseapy's default when no background is given). The larger background "
            "makes every overlap look rarer."
        ),
        horizontal=True,
        width="double",
    )


def run(  # noqa: PLR0913 - one call from the CLI, every argument is a named option or context
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset: str,
    *,
    dataset_file: Path,
) -> list[Path]:
    """Enrichment of every contrast; the M6 tables, figures and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: Slug naming the dataset.
        dataset_file: The dataset file, with ``[function]`` and ``[de]``.

    Returns:
        The page, then every table and figure file.

    Raises:
        DatasetError: If the dataset file lacks ``[function]``, or an input is missing.
    """
    spec = load(dataset_file)
    if spec.slug != dataset or spec.function is None or spec.de is None:
        raise DatasetError(f"{dataset_file}: needs slug {dataset!r}, [de] and [function]")
    fn = spec.function
    gmt = settings.root / spec.resource(fn.gene_sets).dest
    if not gmt.is_file():
        raise DatasetError(f"{gmt} is missing; run ctx setup --dataset-file {dataset_file}")
    sets = lite.read_gmt(gmt)
    de_dir = settings.root / "outputs" / f"{dict(MODULES)['de']:02d}_de"
    paths = {c.name: newest(de_dir, dataset, f"{c.name}_results") for c in spec.contrasts}
    results = {name: pl.read_parquet(path) for name, path in paths.items()}
    day = date.today()
    releases = {
        "gene sets": f"{fn.gene_sets} ({gmt.name}); {fn.licence}",
        **{
            f"sha256:{p.name}": hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (gmt, *paths.values())
        },
    }
    writer = Writer(
        settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset, day, releases
    )
    writer.blocks.append(f"Gene sets: {gmt.name}, {len(sets)} sets. {fn.licence}")
    gsea_tables, metrics, enriched = {}, [], {}
    for name, frame in results.items():
        found = _gsea(frame, sets)
        enriched[name] = {r.name: r for r in found}
        table = pl.DataFrame(
            [(r.name, r.size, r.es, r.nes, r.p, r.fdr) for r in found],
            schema=list(GSEA_COLUMNS),
            orient="row",
        )
        gsea_tables[name] = table
        writer.table(table, f"{name}_gsea", GSEA_COLUMNS)
        background = _symbols(frame)
        writer.table(_ora(frame, sets, background), f"{name}_ora", ORA_COLUMNS)
        hit = enriched[name].get(fn.featured_set)
        where = f"{fn.featured_set} in {name}"
        metrics += (
            [
                (f"gsea.{name}.{fn.featured_set}.nes", hit.nes, where),
                (f"gsea.{name}.{fn.featured_set}.fdr", hit.fdr, where),
            ]
            if hit
            else [(f"gsea.{name}.{fn.featured_set}.nes", float("nan"), f"{where}: not tested")]
        )
    featured = spec.de.featured
    writer.figure(
        _ora_bars(
            featured,
            _ora(results[featured], sets, _symbols(results[featured])),
            len(_symbols(results[featured])),
        ),
        f"ora_{featured}",
        "M6 ORA",
    )
    writer.figure(_gsea_bars(featured, gsea_tables[featured]), f"gsea_{featured}", "M6 GSEA")
    tested = [n for n, table in gsea_tables.items() if fn.featured_set in set(table["gene_set"])]
    if tested:
        featured_bars = _featured_bars(fn.featured_set, gsea_tables)
        writer.figure(featured_bars, f"nes_{_slug(fn.featured_set)}", "M6 one set, every contrast")
    else:
        writer.blocks.append(
            f"{fn.featured_set} was tested in no contrast: after the cut to the ranked genes its "
            f"size is outside {lite.MIN_SIZE} to {lite.MAX_SIZE}, or it is not in {gmt.name}."
        )
    contrast, gene_set = fn.running
    if gene_set in enriched[contrast]:
        running = _running(contrast, enriched[contrast][gene_set])
        writer.figure(running, f"running_{contrast}", "M6 running score")
    else:
        writer.blocks.append(f"No running score: {gene_set} was not tested in {contrast}.")
    writer.figure(
        _trap(featured, results[featured], sets), f"background_{featured}", "M6 background trap"
    )
    writer.table(
        pl.DataFrame(metrics, schema=["metric", "value", "detail"], orient="row"),
        "metrics",
        {
            "metric": ("none", "Name the sanity rules read."),
            "value": ("varies", "The value."),
            "detail": ("none", "What it was computed on."),
        },
    )
    usage = runtime.measure(settings.started)
    budget = runtime.budget(MODULE)
    page = ModulePage(
        title=f"M6 Function: {dataset}",
        question=QUESTION,
        concept=CONCEPT,
        parameters=cards,
        results=tuple(writer.blocks),
        practice=PRACTICE,
        think=THINK,
        mistakes=MISTAKES,
        reproduce=Reproduction(
            command=command,
            wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=f"Track B: under {budget.minutes} min and {budget.ram_gb} GB",
            versions=runtime.versions(),
        ),
    )
    return [
        write(page, writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"),
        *writer.files,
    ]


def _slug(name: str) -> str:
    """A gene set name as an output-name token."""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
