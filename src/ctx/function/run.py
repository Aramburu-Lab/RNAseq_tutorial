"""M6: ``ctx function``. Enrichment of every contrast's results, then tables, figures, page.

Reads the newest M5 results of each declared contrast and the gene sets ``ctx setup``
fetched. GSEA ranks by the Wald statistic, one gene per symbol; ORA counts the up and the
down calls separately against the genes DESeq2 tested.
"""

import hashlib
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, replace
from datetime import date
from pathlib import Path
from typing import Any

import polars as pl

from ctx import runtime
from ctx.dataset import Contrast, Dataset, DatasetError, load
from ctx.de.page import Words, every
from ctx.de.run import Writer, newest
from ctx.de.stats import ALPHA
from ctx.function import background, lite, live, running
from ctx.function.page import (
    CALLS,
    PAGE,
    fdr_words,
    figure,
    gsea_text,
    nes_text,
    ora_text,
    questions,
    set_words,
    trap_text,
    union_direction,
)
from ctx.gather import gdc
from ctx.io.naming import MODULES, stem
from ctx.qc.run import GTF_RESOURCE, gene_types
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
    SelfCheck,
    write,
)
from ctx.report.entry import write_entry
from ctx.report.lookups import lookup
from ctx.report.text import fill, word
from ctx.report.views import viewed
from ctx.viz.bar import BarData
from ctx.viz.style import DATA_COLOURS, ROLES, THEMES

MODULE = "function"
TOP_SETS = 15
# A metric a [sanity] rule names for one gene set: gsea.<contrast>.<SET>.nes or .fdr.
GSEA_METRIC = re.compile(r"gsea\.(?P<contrast>[^.]+)\.(?P<gene_set>[^.]+)\.(?P<what>nes|fdr)")
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


def _tested(results: pl.DataFrame, call: str | None = None) -> pl.DataFrame:
    """The tested genes, or those of them with one call."""
    tested = set(lite.background(results["gene_id"].to_list(), results["padj"].to_list()))
    rows = results.filter(pl.col("gene_id").is_in(list(tested)))
    return rows.filter(pl.col("call") == call) if call else rows


def _symbols(results: pl.DataFrame, call: str | None = None) -> list[str]:
    """Symbols of the tested genes, or of the genes with one call; two genes may share one."""
    return sorted({s for s in _tested(results, call)["gene_name"].drop_nulls() if s})


def _gsea(results: pl.DataFrame, sets: Mapping[str, tuple[str, ...]]) -> list[lite.Enriched]:
    """Preranked GSEA of one contrast."""
    names, scores = running.ranking(results)
    return lite.prerank(names, scores, sets)


def rule_sets(
    rules: Mapping[str, Any], enriched: Mapping[str, Mapping[str, lite.Enriched]]
) -> list[tuple[str, float, str]]:
    """The NES and FDR of each gene set a ``[sanity]`` rule names, as metrics.

    A rule names one as ``gsea.<contrast>.<SET>.nes`` or ``.fdr``. A value the run cannot give
    is NaN with its reason, which the hub reports as NOT MEASURED.
    """
    named = {
        m
        for rule in rules.get("rules", [])
        for c in rule.get("conditions", [])
        for m in (c.get("metric"), c.get("other"))
        if m and GSEA_METRIC.fullmatch(m)
    }
    metrics = []
    for metric in sorted(named):
        found = GSEA_METRIC.fullmatch(metric)
        assert found is not None  # noqa: S101 - kept by the filter above
        contrast, gene_set, what = found["contrast"], found["gene_set"], found["what"]
        hit = enriched.get(contrast, {}).get(gene_set)
        where = f"{gene_set} in {contrast}"
        reason = (
            f"no contrast {contrast}" if contrast not in enriched
            else f"{where}: not tested" if hit is None
            else where
        )  # fmt: skip
        metrics.append((metric, math.nan if hit is None else float(getattr(hit, what)), reason))
    return metrics


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


def _passing() -> tuple[str, str]:
    """The words of a set that passes FDR < ALPHA, and of one that does not."""
    return text.fill(PAGE, "significance.passing", alpha=ALPHA), word(PAGE, "significance.failing")


def _labelled(frame: pl.DataFrame) -> pl.DataFrame:
    """The frame with each gene set in words, as the bars name it."""
    return frame.with_columns(label=pl.col("gene_set").map_elements(set_words, pl.String))


def _gsea_bars(said: Words, frame: pl.DataFrame) -> BarData:
    """NES of the sets with the smallest FDR."""
    passing, failing = _passing()
    top = _labelled(frame.sort("fdr", "gene_set").head(TOP_SETS)).with_columns(
        significance=pl.when(pl.col("fdr") < ALPHA).then(pl.lit(passing))
        .otherwise(pl.lit(failing)),
    )  # fmt: skip
    columns = {
        **GSEA_COLUMNS,
        "label": ("none", "Set name in words."),
        "significance": ("none", f"FDR below {ALPHA} or not."),
    }
    # The sets drawn, which fall short of TOP_SETS where fewer were tested.
    key = "figures.gsea.caption_one" if top.height == 1 else "figures.gsea.caption"
    return BarData(
        frame=top,
        columns=columns,
        category="label",
        value="nes",
        category_label=figure("gsea", "x"),
        value_label=figure("gsea", "y"),
        series="significance",
        colours={passing: ROLES["comparison"], failing: None},
        caption=fill(
            PAGE,
            key,
            label=said.label,
            permutations=lite.PERMUTATIONS,
            seed=lite.SEED,
            shown=top.height,
        ),  # fmt: skip
        horizontal=True,
        hover=("gene_set", "fdr", "size"),
        width="double",
    )


def _ora_bars(said: Words, frame: pl.DataFrame, background: int) -> BarData | None:
    """-log10 FDR of the sets ``live.bars`` draws, a block per direction; None when none is.

    A direction with fewer than ``live.TOP`` passing sets adds its failing sets of smallest
    p in grey (PK.3), past which the dashed line marks FDR = ``ALPHA``. The length is FDR's
    (Rytis 2026-10-21, D-102), so a grey bar is short: Step-1's 5 grey down sets 1.6e-13.
    """
    top = _labelled(live.bars(frame)).unique(["label", "series"], maintain_order=True)
    if top.is_empty():
        return None
    columns = {
        **ORA_COLUMNS,
        "passes": ("none", f"FDR below {ALPHA} or not."),
        "series": ("none", "The calls counted; _failing where the set does not pass, grey."),
        live.LENGTH: ("-log10 FDR", "Minus log10 of the FDR, the bar's length."),
        "label": ("none", "Set name in words."),
    }
    families = text.raw(PAGE, "live.families")
    # The page's own legend words (ctxOra.blocks names its traces from the same keys).
    names = {}
    for d in ("up", "down"):
        names[d] = fill(PAGE, "live.legend.passing", family=families[d], alpha=ALPHA)
        names[f"{d}_failing"] = fill(PAGE, "live.legend.failing", family=families[d], alpha=ALPHA)
    return BarData(
        frame=top,
        columns=columns,
        category="label",
        value=live.LENGTH,
        category_label=figure("ora", "x"),
        value_label=figure("ora", "y"),
        series="series",
        colours={
            "up": ROLES["comparison"],
            "down": ROLES["reference"],
            "up_failing": None,
            "down_failing": None,
        },  # fmt: skip
        caption=figure(
            "ora",
            "caption",
            label=said.label,
            calls=CALLS,
            background=background,
            alpha=ALPHA,
            top=live.TOP,
        ),  # fmt: skip
        horizontal=True,
        hline=(-math.log10(ALPHA), figure("ora", "line", alpha=ALPHA)),
        hover=("gene_set", "p", "fdr", "hits", "size"),
        width="double",
        names=names,
        labels={"series": word(PAGE, "live.direction")},
    )


def _featured_bars(
    featured: str, per_contrast: Mapping[str, pl.DataFrame], said: Mapping[str, Words]
) -> BarData:
    """One set's NES in every contrast, each contrast by its label."""
    passing, failing = _passing()
    rows = []
    for name, frame in per_contrast.items():
        found = frame.filter(pl.col("gene_set") == featured)
        if found.height:
            nes, fdr = found.row(0, named=True)["nes"], found.row(0, named=True)["fdr"]
            rows.append((name, said[name].label, nes, fdr, passing if fdr < ALPHA else failing))
    frame = pl.DataFrame(rows, schema=["contrast", "label", "nes", "fdr", "significance"],
                         orient="row")  # fmt: skip
    columns = {
        "contrast": ("none", "Contrast."),
        "label": ("none", "The contrast in words, as the page names it."),
        "nes": GSEA_COLUMNS["nes"],
        "fdr": GSEA_COLUMNS["fdr"],
        "significance": ("none", f"FDR below {ALPHA} or not."),
    }
    return BarData(
        frame=frame,
        columns=columns,
        category="label",
        value="nes",
        category_label=figure("nes", "x"),
        value_label=figure("nes", "y"),
        series="significance",
        colours={passing: ROLES["comparison"], failing: None},
        caption=figure("nes", "caption", gene_set=set_words(featured)),
        horizontal=True,
        hover=("contrast", "fdr"),
        width="double",
    )


def _background_views(
    writer: Writer, name: str, said: Words, results: pl.DataFrame, e: "_Enrichment"
) -> tuple[FigureBlock, tuple[int, int, int, int, int]]:
    """The background lesson, with every background, call and set count as a view (WO7).

    The build draws the up calls against the tested genes and the gene-set union; the page
    looks every other view up.

    Returns:
        The figure, and the build's tested and union sizes, sets moved more and less, and
        sets drawn, for the page's words.
    """
    tested = _symbols(results)
    wide = background.backgrounds(tested, e.sets, e.types)
    found = {(b, d): background.tested_against(_symbols(results, d), tested, wide[b], e.sets, d, b)
             for b in background.BACKGROUNDS if b in wide
             for d in background.DIRECTIONS}  # fmt: skip
    options, at = [], []
    kept = [b for b in background.BACKGROUNDS if b in wide]
    for i, b in enumerate(kept):
        for j, d in enumerate(background.DIRECTIONS):
            for k, how in enumerate(background.SHOWN):
                t = found[(b, d)]
                calls = _tested(results, d).height
                bars = background.bars(said, calls, len(tested), t, how)
                words_ = trap_text(len(tested), t.size, t.moved(how)[:2], t.moved(how)[2],
                                   t.passing())  # fmt: skip
                options.append((background.option(b, d, how), bars, words_))
                at.append((i, j, k))
    first = options[0]
    writer.table(background.table(list(found.values())), f"background_{name}_views",
                 background.TABLE)  # fmt: skip
    block = writer.figure(first[1], f"background_{name}", first[2])
    built = found[(kept[0], background.DIRECTIONS[0])]
    more, less, shown = built.moved(background.SHOWN[0])
    sizes = (len(tested), built.size, more, less, shown)
    axes = background.axes(kept)
    return lookup(block, options, figure("background", "view"), axes, at), sizes


@dataclass(frozen=True)
class _Enrichment:
    """What one run computed, for the page: results, gene sets, GSEA tables and hits."""

    results: Mapping[str, pl.DataFrame]
    sets: Mapping[str, tuple[str, ...]]
    gsea: Mapping[str, pl.DataFrame]
    enriched: Mapping[str, Mapping[str, lite.Enriched]]
    gmt: Path
    types: pl.DataFrame | None = None


def _gsea_views(writer: Writer, featured: str, labelled: Mapping[str, Words],
                e: _Enrichment) -> FigureBlock:  # fmt: skip
    """The featured contrast's GSEA bars, with every other contrast's as a view.

    Work-order item 7 (D-051): passing sets run from 36 of 49 in mono_gcsf to 2 of 48 in the
    null split (figure controls plan); nothing is recomputed in the page, and each contrast's
    numbers are its own GSEA table beside the page.
    """
    options = []
    for name in [featured, *(n for n in e.gsea if n != featured)]:
        gsea = e.gsea[name]
        ranked = _gsea_bars(labelled[name], gsea)
        drawn = (int((ranked.frame["fdr"] < ALPHA).sum()), ranked.frame.height)
        words_ = gsea_text(labelled[name], int((gsea["fdr"] < ALPHA).sum()), gsea.height, drawn)
        options.append((labelled[name].label, ranked, words_))
    block = writer.figure(options[0][1], f"gsea_{featured}", options[0][2])
    return viewed(block, options, figure("gsea", "view")) if len(options) > 1 else block


def _featured(writer: Writer, contrast: Contrast, labelled: Mapping[str, Words],
              e: _Enrichment) -> Group:  # fmt: skip
    """The featured contrast's ORA, live on the reader's list (GS.5), and GSEA bars.

    A sentence stands for the ORA bars where no set passes on the build's calls; the live
    panel rides on that figure, so it is absent then too.
    """
    name = contrast.name
    background = _symbols(e.results[name])
    ora = _ora(e.results[name], e.sets, background)
    passing = dict(ora.filter(pl.col("fdr") < ALPHA)["direction"].value_counts().rows())
    said = labelled[name]
    bars = _ora_bars(said, ora, len(background))
    if bars is not None:
        labels = {n: w.label for n, w in labelled.items() if n in e.results}
        labels["union"] = word(PAGE, "live.union")
        panel = live.panel(e.results, name, e.sets, labels,
                           {s: set_words(s) for s in e.sets})  # fmt: skip
        # "all": the list as one family, in the palette's orange, which neither call takes;
        # "failing": the grey of a set that does not pass, whatever its family (PK.3).
        grey = THEMES["white"].reference
        panel["colours"] = {"up": ROLES["comparison"], "down": ROLES["reference"],
                            "all": DATA_COLOURS[0], "failing": grey}  # fmt: skip
        # The redraw may name any set: a phone's layout is rebuilt for its labels (PK.4).
        bars = replace(bars, live=panel, drawable=tuple(panel["labels"].values()))
    items = (
        writer.figure(bars, f"ora_{name}", ora_text(said, len(background), passing.get("up", 0),
                                                    passing.get("down", 0),
                                                    int((~bars.frame["passes"]).sum())))
        if bars else fill(PAGE, "figures.ora.none"),
        _gsea_views(writer, name, labelled, e),
    )  # fmt: skip
    # Tabs hold figures only: with no ORA figure its sentence stands above the GSEA bars.
    return Group(fill(PAGE, "groups.featured", label=said.label), items, tabs=bars is not None)


def _page(
    writer: Writer, spec: Dataset, e: _Enrichment
) -> tuple[list[Result], tuple[KeyNumber, ...], tuple[SelfCheck, ...]]:
    """Every M6 figure, grouped as the page shows them, and the key numbers."""
    assert spec.de is not None and spec.function is not None  # noqa: S101 - run() checked
    fn = spec.function
    featured = spec.de.featured
    contrasts = {c.name: c for c in spec.contrasts}
    said = every(spec.contrasts, spec.columns)
    gsea = e.gsea[featured]
    first = _featured(writer, contrasts[featured], said, e)
    one_set: list[FigureBlock | str] = []
    tested = [n for n, table in e.gsea.items() if fn.featured_set in set(table["gene_set"])]
    if tested:
        featured_bars = _featured_bars(fn.featured_set, e.gsea, said)
        one_set.append(writer.figure(featured_bars, f"nes_{_slug(fn.featured_set)}",
                                     nes_text(fn.featured_set, len(tested))))  # fmt: skip
    else:
        one_set.append(fill(PAGE, "groups.untested", gene_set=set_words(fn.featured_set),
                            low=lite.MIN_SIZE, high=lite.MAX_SIZE, file=e.gmt.name))  # fmt: skip
    contrast, gene_set = fn.running
    if gene_set in e.enriched[contrast]:
        runs = {c: (said[c], e.results[c], e.enriched[c]) for c in e.results}
        rest = [r for c, r in runs.items() if c != contrast]
        one_set.append(running.block(writer, contrast, runs[contrast], gene_set, rest))
    else:
        one_set.append(fill(PAGE, "groups.no_running", gene_set=set_words(gene_set),
                            label=said[contrast].label))  # fmt: skip
    trap, (n_tested, n_union, more, less, _) = _background_views(
        writer, featured, said[featured], e.results[featured], e
    )
    page_results: list[Result] = [
        *_acknowledged(spec),
        ResultTable(word(PAGE, "sets.title"), text.words(PAGE, "sets.columns"),
                    ((e.gmt.name, f"{len(e.sets):,}", fn.licence),), numeric=(1,)),
        first,
        Group(fill(PAGE, "groups.one_set", gene_set=set_words(fn.featured_set)), tuple(one_set)),
        Group(word(PAGE, "groups.background"), (trap,)),
    ]  # fmt: skip
    label = said[featured].label
    keys: tuple[KeyNumber, ...] = (
        KeyNumber(word(PAGE, "keys.sets"), f"{len(e.sets):,}", e.gmt.name),
        KeyNumber(word(PAGE, "keys.contrasts"), f"{len(e.results)}"),
        KeyNumber(fill(PAGE, "keys.passing", alpha=ALPHA), f"{int((gsea['fdr'] < ALPHA).sum())}",
                  fill(PAGE, "keys.passing_detail", label=label)),
    )  # fmt: skip
    shown = e.enriched[featured].get(fn.featured_set)
    if shown:
        detail = fill(PAGE, "keys.nes_detail", fdr=fdr_words(shown.fdr), label=label)
        keys += (KeyNumber(fill(PAGE, "keys.nes", gene_set=set_words(fn.featured_set)),
                           f"{shown.nes:.2f}", detail),)  # fmt: skip
    direction = None if n_union == n_tested else union_direction(more, less)
    return page_results, keys, questions(direction)


def run(  # noqa: PLR0913 - one call from the CLI, every argument is a named option or context
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    command: str,
    dataset: str,
    *,
    dataset_file: Path,
    genes_from: Path | None = None,
) -> list[Path]:
    """Enrichment of every contrast; the M6 tables, figures and page.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        command: The command line, as the reader should retype it.
        dataset: Slug naming the dataset.
        dataset_file: The dataset file, with ``[function]`` and ``[de]``.
        genes_from: A gene list (``ctx.report.genes.read_list``): its ORA table is written,
            as the page's live panel tests the same list (GS.5).

    Returns:
        The page, then every table and figure file.

    Raises:
        DatasetError: If the dataset file lacks ``[function]``, or an input is missing.
        GeneListError: If ``genes_from`` is no list file.
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
    # The background lesson's annotation backgrounds (WO7); without the GTF, the union alone.
    types, _ = gene_types(settings.root, spec)
    annotation = [settings.root / r.dest for r in spec.resources if r.name == GTF_RESOURCE]
    day = date.today()
    releases = {
        "gene sets": f"{fn.gene_sets} ({gmt.name}); {fn.licence}",
        **{
            f"sha256:{p.name}": hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (gmt, *paths.values(), *(annotation if types is not None else []))
        },
    }
    writer = Writer(
        settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset, day, releases
    )
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
    if genes_from is not None:
        writer.table(live.listed(genes_from, results[spec.de.featured], sets),
                     f"{spec.de.featured}_ora_list", ORA_COLUMNS)  # fmt: skip
    # The featured set's rows stay first and as they were; a rule adds only what they lack.
    recorded = {m for m, _, _ in metrics}
    metrics += [m for m in rule_sets(spec.sanity, enriched) if m[0] not in recorded]
    page_results, keys, checks = _page(
        writer, spec, _Enrichment(results, sets, gsea_tables, enriched, gmt, types)
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
        title=fill(PAGE, "title", dataset=spec.name),
        **(text.skeleton(PAGE, alpha=ALPHA) | {"questions": checks}),
        parameters=cards,
        results=tuple(page_results),
        key_numbers=keys,
        reproduce=Reproduction(
            command=command,
            wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(),
        ),
    )
    return [
        write(page, writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html"),
        write_entry(writer.out, MODULE),
        *writer.files,
    ]


def _acknowledged(spec: Dataset) -> tuple[str, ...]:
    """The GDC acknowledgement for a page drawn from GDC data, as M5's page prints it."""
    return gdc.acknowledgement(spec.gdc) if spec.gdc is not None else ()


def _slug(name: str) -> str:
    """A gene set name as an output-name token."""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
