"""The M7 page: notes, the missingness figures first, then the tests, then the sanity rules.

``page.missingness_first`` is measured on the page as built, not assumed: it is 1 only when the
first two figures are the quantified counts and the coverage curves (rule K3 and P3).
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from ctx import concordance, runtime
from ctx.dataset import Dataset
from ctx.de.run import Writer
from ctx.de.sanity import VERDICT_COLUMNS
from ctx.de.stats import ALPHA
from ctx.io.naming import stem
from ctx.proteo import figures, measure
from ctx.proteo.run import Measured, Options, verdicts
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
    write,
)

PAGE = MODULE = "proteo"
MISSINGNESS = ("quantified", "coverage")
# The page's figures before any test: the quantified counts, then the coverage curves.
FIRST = len(MISSINGNESS)
METRIC_COLUMNS = {
    "metric": ("none", "Name the [proteo.sanity] rules read."),
    "value": ("varies", "The measured value; the unit follows from the metric's name."),
    "detail": ("none", "What the value was computed on."),
}


@dataclass(frozen=True)
class Context:
    """Everything the page reads."""

    spec: Dataset
    options: Options
    measured: Measured
    settings: runtime.Settings
    command: str
    cards: tuple[ParameterCard, ...]


def _p(value: float) -> str:
    """A p-value as the captions print it."""
    return f"{value:.2g}"


def _notes(ctx: Context) -> list[Result]:
    """The source lines: cross-study label, the collapse, zeros read as missing, symbols."""
    m, p = ctx.measured, ctx.spec.proteo
    assert p is not None and m.protein.collapse is not None  # noqa: S101 - run() checks
    c = m.protein.collapse
    out: list[Result] = []
    if p.design == "cross_dataset":
        out.append(text.fill(PAGE, "cross_study", rna=ctx.spec.name,
                             protein=p.label or p.protein))  # fmt: skip
    else:
        out.append(text.word(PAGE, "unregistered"))
    out.append(text.fill(PAGE, "collapse", rows=c.rows, no_symbol=c.no_symbol,
                         several=c.several_genes, twice=c.genes_twice, kept=c.kept))  # fmt: skip
    if m.protein.zeros:
        out.append(text.fill(PAGE, "zeros", zeros=m.protein.zeros))
    sources = "; ".join(
        f"{figures.class_name(n)}: {', '.join(g.sets)}" for n, g in p.classes.items()
    )
    for key, values in m.notes:
        out.append(text.fill(PAGE, key, **values, sources=sources))
    return out


def _missingness(writer: Writer, ctx: Context) -> list[FigureBlock]:
    """F1 and F2, drawn before any test."""
    m, p = ctx.measured, ctx.spec.proteo
    assert p is not None  # noqa: S101 - run() checks
    colour = m.extra.get("colour", m.groups)
    assert isinstance(colour, list)  # noqa: S101 - across() sets a list of track values
    label = text.word(PAGE, f"abundance.{p.design}")
    return [figures.quantified(writer, m.protein, colour, m.groups, ctx.options.min_quantified),
            figures.coverage(writer, m.genes, m.abundance, ctx.options.min_quantified,
                             label)]  # fmt: skip


def _across(writer: Writer, ctx: Context) -> list[Result]:
    """F3 and the worked genes."""
    m, p = ctx.measured, ctx.spec.proteo
    assert p is not None  # noqa: S101 - run() checks
    tested = figures.shown(m.genes.filter(pl.col("tested")), ctx.options.gene_class)
    test = m.extra["class_test"]
    below = m.extra.get("complex_below", (np.nan, np.nan))
    values = {"n_samples": m.extra["n_samples"], "n_genes": tested.height, "alpha": ALPHA,
              "median": measure.median(tested["rho"]), "h": test.h, "p": _p(test.p),
              "eps": test.epsilon2, "p_class": _p(below[0]),
              "n_sig": tested.filter((pl.col("rho_padj") < ALPHA)
                                     & (pl.col("rho") > 0)).height}  # fmt: skip
    out: list[Result] = [Group(text.word(PAGE, "groups.correlation"),
                               (figures.rho_by_class(writer, tested, values),))]  # fmt: skip
    rna, protein = m.extra["rna"], m.extra["values"]
    annotation = m.extra["annotation"]
    rows = {g: i for i, g in enumerate(m.genes["gene"])}
    blocks = []
    tested_rows = set(m.genes.filter("tested")["gene"])
    for gene in [g for g in p.genes if g in tested_rows]:
        i = rows[gene]
        keep = ~np.isnan(protein[i])
        frame = pl.DataFrame(
            {
                "sample": annotation["sample"].filter(pl.Series(keep)),
                "rna": rna[i][keep],
                "protein": protein[i][keep],
                "track": annotation[p.amplified or p.plex].filter(pl.Series(keep)),
            }
        )
        row = m.genes.row(i, named=True)
        blocks.append(figures.worked(writer, gene, frame, {
            "rho": row["rho"], "r": row["r"], "n": frame.height, "track": p.amplified,
            "padj": _p(row["rho_padj"])}))  # fmt: skip
    if blocks:
        out.append(Group(text.word(PAGE, "groups.worked"), tuple(blocks)))
    return out


def _cross(writer: Writer, ctx: Context) -> list[Result]:
    """F7 and F8."""
    m = ctx.measured
    tested = m.extra["tested"]
    assert isinstance(tested, pl.DataFrame)  # noqa: S101 - fold_metrics sets it
    rho, low, high = m.extra["fold"]
    contrast = m.extra["contrast"]
    same = int(((tested["rna_lfc"] * tested["protein_lfc"]) > 0).sum())
    values = {"rho": rho, "low": low, "high": high, "n_genes": tested.height, "n_same": same,
              "tested": contrast.tested, "reference": contrast.reference,
              "n_rna": " and ".join(map(str, contrast.expected_n)),
              "n_protein": f"{m.extra['n_tested']} and {m.extra['n_reference']}",
              "ci": round(100 * concordance.LEVEL)}  # fmt: skip
    blocks: list[FigureBlock] = [figures.fold(writer, tested, values)]
    assert ctx.spec.proteo is not None  # noqa: S101 - run() checks
    named = [c for c in ctx.spec.proteo.classes if f"shift.{c}" in m.extra]
    if named:
        p, r, _, _ = m.extra[f"shift.{named[0]}"]
        inside = measure.median(tested.filter(pl.col("gene_class") == named[0])["difference"])
        other = measure.median(tested.filter(pl.col("gene_class") == measure.OTHER)["difference"])
        word = "above" if inside > other else "below"
        shown = figures.shown(tested, ctx.options.gene_class)
        blocks.append(figures.shift(writer, shown, {
            "p": _p(p), "r": r, "granule": inside, "other": other,
            "direction": text.word(PAGE, f"figures.shift.{word}")}))  # fmt: skip
    return [Group(text.word(PAGE, "groups.fold"), tuple(blocks))]


def _sanity(writer: Writer, ctx: Context) -> ResultTable:
    """The registered rules, decided on this run's metrics, written and shown."""
    assert ctx.spec.proteo is not None  # noqa: S101 - run() checks
    rows = verdicts(ctx.spec.proteo, ctx.measured.metrics)
    writer.table(pl.DataFrame(rows, schema=dict.fromkeys(VERDICT_COLUMNS, pl.String),
                              strict=False), "sanity", VERDICT_COLUMNS)  # fmt: skip
    shown = []
    for row in rows:
        left = "" if row["left"] is None else f"{row['left']:.3g}"
        shown.append((str(row["rule"]), str(row["text"]), f"{row['condition']}: {left}",
                       str(row["verdict"])))  # fmt: skip
    return ResultTable(text.word(PAGE, "sanity_title"), text.words(PAGE, "sanity_columns"),
                       tuple(shown))  # fmt: skip


def _keys(ctx: Context) -> tuple[KeyNumber, ...]:
    """Genes in both layers, tested genes, and the headline correlation."""
    m = ctx.measured
    tested = int(m.genes["tested"].sum())
    keys = [KeyNumber(text.word(PAGE, "keys.genes"), f"{m.genes.height:,}"),
            KeyNumber(text.word(PAGE, "keys.tested"), f"{tested:,}",
                      text.fill(PAGE, "keys.tested_detail",
                                pct=round(100 * ctx.options.min_quantified)))]  # fmt: skip
    if "fold" in m.extra:
        rho, low, high = m.extra["fold"]
        keys.append(KeyNumber(text.word(PAGE, "keys.fold"), f"{rho:.2f}",
                              text.fill(PAGE, "keys.fold_detail", low=low, high=high)))  # fmt: skip
    else:
        median = measure.median(m.genes.filter(pl.col("tested"))["rho"])
        keys.append(KeyNumber(text.word(PAGE, "keys.rho"), f"{median:.2f}"))
    return tuple(keys)


def build(writer: Writer, ctx: Context) -> list[Path]:
    """Write every table and figure, decide the rules, then the page; the page first."""
    m = ctx.measured
    missing = _missingness(writer, ctx)
    across = ctx.spec.proteo is not None and ctx.spec.proteo.design == "across_samples"
    tests = _across(writer, ctx) if across else _cross(writer, ctx)  # fmt: skip
    results: list[Result] = [*_notes(ctx), Group(text.word(PAGE, "groups.missing"), tuple(missing)),
                             *tests]  # fmt: skip
    blocks = [b for r in results for b in (r.items if isinstance(r, Group) else (r,))
              if isinstance(b, FigureBlock)]  # fmt: skip
    first = all(b.stem.endswith(f"_{w}_v1") for b, w in zip(blocks, MISSINGNESS, strict=False))
    m.metrics.append(("page.missingness_first", float(first and len(blocks) >= FIRST),
                      "the first two figures are the quantified counts and coverage"))  # fmt: skip
    writer.table(m.genes, "genes", figures.columns(*m.genes.columns))
    writer.table(pl.DataFrame(m.metrics, schema=list(METRIC_COLUMNS), orient="row"), "metrics",
                 METRIC_COLUMNS)  # fmt: skip
    results.append(_sanity(writer, ctx))
    usage, budget = runtime.measure(ctx.settings.started), runtime.budget(MODULE)
    track = "A" if MODULE in runtime.TRACK_A_MODULES else "B"
    built = ModulePage(
        title=text.fill(PAGE, "title", dataset=ctx.spec.name), **text.skeleton(PAGE),
        parameters=ctx.cards, results=tuple(results), key_numbers=_keys(ctx),
        reproduce=Reproduction(
            command=ctx.command, wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self, peak_rss_children=usage.peak_rss_children,
            profile=f"{ctx.settings.profile}, {ctx.settings.threads} threads",
            budget=text.fill(PAGE, "budget", track=track, minutes=budget.minutes,
                             ram_gb=budget.ram_gb),
            versions=runtime.versions(),
        ),
    )  # fmt: skip
    page = writer.out / f"{stem(PAGE, f'{writer.dataset}_overview', 1, day=writer.day)}.html"
    return [write(built, page), *writer.files]
