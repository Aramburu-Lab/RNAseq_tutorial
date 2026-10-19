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
from ctx.de.stats import ALPHA, LFC_MIN
from ctx.io.naming import stem
from ctx.proteo import densities, figures, folds, measure
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
# F5 draws the named genes and this many of the highest and of the lowest per-gene rho.
PAIRED_TOP = 8
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
    elif p.design == "across_samples" and p.table is None:
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
    values["class_test"] = text.fill(PAGE, "figures.rho.class_test", p_class=_p(below[0])
                                     ) if "complex_below" in m.extra else ""  # fmt: skip
    normals = "normal_class_test" in m.extra
    drawn = [figures.rho_by_class(writer, tested, values,
                                  "figures.rho_tumour" if normals else "figures.rho")]  # fmt: skip
    if normals:
        drawn.append(_normal_rho(writer, ctx))
    out: list[Result] = [*_shifts(ctx), Group(text.word(PAGE, "groups.correlation"), tuple(drawn))]
    rna, protein = m.extra["rna"], m.extra["values"]
    annotation = m.extra["annotation"]
    rows = {g: i for i, g in enumerate(m.genes["gene"])}
    blocks = []
    tested_rows = set(m.genes.filter("tested")["gene"])
    # Tumour against normal draws both tissues, coloured by tissue; rho stays the tumours'.
    paired = "worked" in m.extra
    shown_rna, shown_protein, shown = m.extra.get("worked", (rna, protein, annotation))
    track = "tissue" if paired else p.amplified or p.plex
    for gene in [g for g in p.genes if g in tested_rows]:
        i = rows[gene]
        keep = ~np.isnan(shown_protein[i])
        frame = pl.DataFrame(
            {
                "sample": shown["sample"].filter(pl.Series(keep)),
                "rna": shown_rna[i][keep],
                "protein": shown_protein[i][keep],
                "track": shown[track].filter(pl.Series(keep)),
            }
        )
        row = m.genes.row(i, named=True)
        blocks.append(figures.worked(writer, gene, frame, {
            "rho": row["rho"], "r": row["r"], "track": track, "padj": _p(row["rho_padj"]),
            "n": int(row["n_quantified"]) if paired else frame.height,
            "unit": p.unit or text.word(PAGE, "unit_ratio")},
            caption=f"figures.worked.caption{'_paired' if paired else ''}"))  # fmt: skip
    if blocks:
        out.append(Group(text.word(PAGE, "groups.worked"), tuple(blocks)))
    named = [*(k for k in p.published.get("mertins", {}) if k.isupper()), *p.genes]
    chosen = figures.heatmap_genes(m.genes, named, PAIRED_TOP)
    if chosen:
        tracks = tuple((c, c) for c in dict.fromkeys([*p.tracks, p.plex]) if c)
        out.append(Group(text.word(PAGE, "groups.paired"), (figures.paired(
            writer, (m.genes, rna, protein), annotation, chosen, tracks),)))  # fmt: skip
    found = [_contrast(writer, made) for made in m.extra.get("folds", [])]
    if found:
        out.append(Group(text.word(PAGE, "groups.contrast"), tuple(found)))
    return out


def _normal_rho(writer: Writer, ctx: Context) -> FigureBlock:
    """F3 of the normals: their own per-gene rho by class, from the ``normal_`` columns."""
    m = ctx.measured
    own = m.genes.filter(pl.col("normal_tested")).select(
        "gene", "gene_class", rho="normal_rho", rho_padj="normal_rho_padj"
    )
    tested = figures.shown(own, ctx.options.gene_class)
    test = m.extra["normal_class_test"]
    values = {"n_samples": m.extra["n_normal"], "n_genes": tested.height, "alpha": ALPHA,
              "median": measure.median(tested["rho"]), "h": test.h, "p": _p(test.p),
              "eps": test.epsilon2,
              "n_sig": tested.filter((pl.col("rho_padj") < ALPHA)
                                     & (pl.col("rho") > 0)).height}  # fmt: skip
    return figures.rho_by_class(writer, tested, values, "figures.rho_normal",
                                "rho_by_class_normal")  # fmt: skip


def _shifts(ctx: Context) -> list[Result]:
    """One line per declared class of the paired design: its fold changes and its shift test."""
    m = ctx.measured
    if "paired_genes" not in m.extra:
        return []
    found = {name: value for name, value, _ in m.metrics}
    out: list[Result] = []
    for name in ctx.spec.proteo.classes if ctx.spec.proteo else ():
        n = m.extra["paired_genes"].filter(pl.col("gene_class") == name).height
        out.append(text.fill(PAGE, "class_shift", name=figures.class_name(name), n=n,
                             protein=found[f"class.{name}.protein_lfc.median"],
                             rna=found[f"class.{name}.rna_lfc.median"],
                             below=_p(found.get(f"class.{name}_below_other.p", np.nan)),
                             above=_p(found.get(f"class.{name}.p", np.nan))))  # fmt: skip
    return out


def _contrast(writer: Writer, made: folds.Folds) -> FigureBlock:
    """F4 for one contrast, its numbers filled into the text."""
    rho, low, high = made.rho
    c = made.contrast
    n_protein = f"{c.reference} {made.n_protein[0]}, {c.tested} {made.n_protein[1]}"
    if made.paired:
        n_protein = text.fill(PAGE, "figures.contrast.pairs", n=made.n_protein[0])
    test = "paired" if made.paired else "welch"
    values = {"tested": c.tested, "reference": c.reference, "n_genes": made.genes.height,
              "n_rna": f"{c.reference} {c.expected_n[0]}, {c.tested} {c.expected_n[1]}",
              "n_protein": n_protein,
              "test": text.word(PAGE, f"figures.contrast.test.{test}"),
              "alpha": ALPHA,
              "lfc": LFC_MIN, "ci": round(100 * concordance.LEVEL),
              "rho": np.nan if rho is None else rho, "low": np.nan if low is None else low,
              "high": np.nan if high is None else high, **made.counts()}  # fmt: skip
    return densities.contrast_density(writer, made, values)


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
    blocks: list[FigureBlock] = [densities.fold(writer, tested, values)]
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
    across = ctx.spec.proteo is not None and ctx.spec.proteo.design != "cross_dataset"
    tests = _across(writer, ctx) if across else _cross(writer, ctx)  # fmt: skip
    results: list[Result] = [*_notes(ctx), Group(text.word(PAGE, "groups.missing"), tuple(missing)),
                             *tests]  # fmt: skip
    blocks = [b for r in results for b in (r.items if isinstance(r, Group) else (r,))
              if isinstance(b, FigureBlock)]  # fmt: skip
    first = all(b.stem.endswith(f"_{w}_v1") for b, w in zip(blocks, MISSINGNESS, strict=False))
    m.metrics.append(("page.missingness_first", float(first and len(blocks) >= FIRST),
                      "the first two figures are the quantified counts and coverage"))  # fmt: skip
    writer.table(m.genes, "genes", figures.columns(*m.genes.columns))
    for what, frame, columns in m.extra.get("tables", ()):
        writer.table(frame, what, columns)
    writer.table(pl.DataFrame(m.metrics, schema=list(METRIC_COLUMNS), orient="row"), "metrics",
                 METRIC_COLUMNS)  # fmt: skip
    results.append(_sanity(writer, ctx))
    usage, budget = runtime.measure(ctx.settings.started), runtime.budget(MODULE)
    track = "A" if MODULE in runtime.TRACK_A_MODULES else "B"
    paired = ctx.spec.proteo is not None and ctx.spec.proteo.design == "tumour_normal"
    sections = text.skeleton(PAGE, "variants.tumour_normal" if paired else None)
    built = ModulePage(
        title=text.fill(PAGE, "title", dataset=ctx.spec.name), **sections,
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
