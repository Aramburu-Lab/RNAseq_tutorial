"""The M8 page of a dataset's Step-1 half: biotypes, the filter, the volcano, the neighbours.

Each group writes its tables and figures through the run's ``Writer``; the numbers come from
``figures`` and ``neighbours``, every sentence from ``content/ncrna.toml``. M8's own sanity
rules (``[ncrna.sanity]``) are decided here, on the metrics this run records.
"""

from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import polars as pl

from ctx import runtime
from ctx.de.run import Writer
from ctx.de.sanity import decide, parse
from ctx.de.stats import ALPHA, LFC_MIN
from ctx.gather import gdc as gdc_reader
from ctx.io.naming import stem
from ctx.ncrna import biotypes, figures, mirna_page
from ctx.report import text
from ctx.report.builder import (
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    write,
)
from ctx.report.views import viewed
from ctx.viz.aggregates import Aggregate, Lead
from ctx.viz.scatter import ScatterData
from ctx.viz.style import ROLES

if TYPE_CHECKING:
    from ctx.ncrna.run import Inputs

PAGE = "ncrna"
# M8 F2's presets of the filter's minimum count (M8 plan: 1 to 50; the run's own first).
FILTER_PRESETS = (1, 2, 5, 10, 20, 50)
METRIC_COLUMNS = {
    "metric": ("none", "Name the sanity rules read."),
    "value": ("varies", "The value."),
    "detail": ("none", "What it was computed on."),
}


@dataclass(frozen=True)
class Context:
    """What the page is built from: the inputs, the options and the run."""

    inputs: "Inputs"
    biotypes: tuple[str, ...]
    min_count: int
    settings: runtime.Settings
    command: str
    cards: tuple[ParameterCard, ...]


def fill(key: str, **values: Any) -> str:
    """A template of this page's content file, filled."""
    return text.fill(PAGE, key, **values)


def word(key: str) -> str:
    """A fixed text of this page's content file."""
    return text.word(PAGE, key)


def figure_text(key: str, **values: Any) -> FigureText:
    """Title, lead and reading guide of ``figures.<key>``, filled from the same values."""
    reading = tuple(line.format_map(values) for line in text.raw(PAGE, f"figures.{key}.reading"))
    return FigureText(title=fill(f"figures.{key}.title", **values),
                      lead=fill(f"figures.{key}.lead", **values), reading=reading)  # fmt: skip


def _filter_panel(ctx: Context, values: pl.DataFrame, k: int, min_count: int
                  ) -> tuple[Any, FigureText, dict[str, float]]:  # fmt: skip
    """F2 at one minimum count: the curves, their words, and each biotype's kept share, %."""
    test = figures.rank_test(values)
    data = figures.curves(values, ctx.biotypes, float(np.log2(min_count + 1)),
                          fill("figures.filter.caption", **contrast_words(ctx),
                               min_count=min_count, **test))  # fmt: skip
    kept = {b: 100 * float(np.mean(values.filter(pl.col("biotype") == b)["kept"].to_numpy()))
            for b in ctx.biotypes}  # fmt: skip
    words = figure_text("filter", protein_coding=kept.get(biotypes.CODING, float("nan")),
                        lncrna=kept.get(biotypes.LNCRNA, float("nan")), k=k,
                        min_count=min_count)  # fmt: skip
    return data, words, kept


def _filter_views(writer: Writer, ctx: Context, columns: Any, sizes: list[int]) -> Any:
    """F2 at the run's minimum count, with the plan's presets as views (M8 F2's control).

    M5's filter is the same rule at any count, so each preset is M5's own filter, computed
    here, never a stand-in; the kept shares of every preset are one table.
    """
    inp = ctx.inputs
    counts = [ctx.min_count, *(m for m in FILTER_PRESETS if m != ctx.min_count)]
    options, rows = [], []
    for m in counts:
        values, k = figures.filter_values(inp.genes, inp.matrix[:, columns], sizes, m)
        data, words, kept = _filter_panel(ctx, values, k, m)
        key = "figures.filter.preset_one" if m == 1 else "figures.filter.preset"
        options.append((fill(key, min_count=m), data, words))
        rows += [(m, b, share) for b, share in kept.items()]
    writer.table(pl.DataFrame(rows, schema=["min_count", "biotype", "kept_pct"], orient="row"),
                 "filter_presets", {
                     "min_count": ("reads", "Minimum count at the median library size."),
                     "biotype": ("none", "Gene biotype."),
                     "kept_pct": ("%", "Share of the biotype's genes M5's filter keeps."),
                 })  # fmt: skip
    block = writer.figure(options[0][1], "filter_by_biotype", options[0][2])
    return viewed(block, options, word("figures.filter.view"))


def contrast_words(ctx: Context) -> dict[str, Any]:
    """The contrast's names and group sizes, which several captions print."""
    c, used = ctx.inputs.contrast, ctx.inputs.used
    n_test = int((used["group"] == "tested").sum())
    cell = (", ".join(used["cell_type"].unique(maintain_order=True).to_list())
            if "cell_type" in used.columns else ctx.inputs.spec.name)  # fmt: skip
    return {"reference": c.reference, "tested": c.tested, "cell_type": cell,
            "n_ref": used.height - n_test, "n_test": n_test, "contrast": c.name}  # fmt: skip


def contrast_columns(ctx: Context) -> np.ndarray:
    """The contrast's samples as columns of the count matrix, reference group first."""
    position = {s: i for i, s in enumerate(ctx.inputs.sheet["sample"].to_list())}
    used = ctx.inputs.used.sort(pl.col("group") != "reference", maintain_order=True)
    return np.array([position[s] for s in used["sample"]])


def library(writer: Writer, ctx: Context) -> tuple[list[Result], dict[str, tuple[float, str]]]:
    """F1 and F2, the miRNA sentence, and the metrics they give."""
    inp = ctx.inputs
    caption = fill("figures.composition.caption", n_samples=inp.sheet.height,
                   n_genes=inp.genes.height, gtf=inp.gtf.name)  # fmt: skip
    bars, wide = figures.composition(inp.sheet, inp.genes, inp.matrix, caption)
    writer.table(wide, "sample_shares", figures.columns(wide))
    where = {"biotype": biotypes.LNCRNA}, {"biotype": biotypes.MITO}
    lead = Lead(text.raw(PAGE, "figures.composition.lead"),
                (Aggregate("lncrna", "median", "percent", where[0]),
                 Aggregate("mt_low", "min", "percent", where[1]),
                 Aggregate("mt_high", "max", "percent", where[1])),
                note=word("figures.composition.note"))  # fmt: skip
    shares = writer.figure(bars, "biotype_shares", replace(figure_text_fixed("composition"),
                                                           follows=lead))  # fmt: skip
    selected = "library_selection" in wide.columns
    poly_a = wide.filter(pl.col("library_selection").is_in(biotypes.POLY_A)) if selected else wide
    top = float(poly_a["mirna_percent"].to_numpy().max()) if poly_a.height else float("nan")
    which = "poly-A" if selected else "every library, its selection not on the sheet"
    columns = contrast_columns(ctx)
    used = ctx.inputs.used
    sizes = [int((used["group"] == "reference").sum()), int((used["group"] == "tested").sum())]
    values, k = figures.filter_values(inp.genes, inp.matrix[:, columns], sizes, ctx.min_count)
    writer.table(values, "filter_values", figures.columns(values))
    curve = _filter_views(writer, ctx, columns, sizes)
    name = inp.contrast.name
    metrics = {
        "mirna.max_percent": (top, f"largest miRNA-gene share of counted reads, {which}, %"),
        "mirna.poly_a_samples": (float(poly_a.height), f"samples read: {which}"),
        "filter.kept": (float(values["kept"].sum()), f"genes the filter keeps at "
                        f"{ctx.min_count} reads, samples of {name}"),
        "m5.tested": (float(inp.results["padj"].is_not_null().sum()), f"genes M5 tested in {name}"),
    }  # fmt: skip
    note = fill("notes.mirna", pct=top)
    return [Group(word("groups.library"), (shares, note, curve))], metrics


def figure_text_fixed(key: str) -> FigureText:
    """Title and reading of a figure whose lead follows its rows."""
    return FigureText(title=word(f"figures.{key}.title"),
                      reading=tuple(text.raw(PAGE, f"figures.{key}.reading")))  # fmt: skip


def de(writer: Writer, ctx: Context) -> tuple[Result, pl.DataFrame]:
    """F3 and the calls table: BH over every tested gene against BH within each biotype."""
    frame = figures.by_biotype(ctx.inputs.results, ctx.inputs.genes, ctx.biotypes)
    writer.table(frame, "results_by_biotype", figures.columns(frame))
    words = contrast_words(ctx)
    n = {b: frame.filter(pl.col("biotype") == b).height for b in ctx.biotypes}
    strata = ", ".join(fill("figures.volcano.stratum", biotype=b, n=n[b]) for b in ctx.biotypes)
    caption = fill("figures.volcano.caption", **words, n_tested=frame.height, strata=strata)
    data = figures.volcano(frame, ctx.biotypes, words, caption)
    called = {b: frame.filter(pl.col("level") == b).height for b in biotypes.TESTED}
    shown = figure_text("volcano", **words, alpha=ALPHA, lfc_min=LFC_MIN,
                        called_protein_coding=called[biotypes.CODING],
                        called_lncrna=called[biotypes.LNCRNA])  # fmt: skip
    block = writer.figure(data, "volcano_by_biotype", shown)
    table = ResultTable(word("figures.volcano.table_title"),
                        text.words(PAGE, "figures.volcano.columns"),
                        tuple(figures.calls_table(frame, ctx.biotypes)),
                        numeric=(1, 2, 3, 4, 5))  # fmt: skip
    return Group(word("groups.de"), (block, table)), frame


def pair_live(shown: pl.DataFrame, a: np.ndarray, b: np.ndarray, strata: np.ndarray,
              members: list[list[int]]) -> dict[str, Any]:  # fmt: skip
    """The pair select's settings: each pair's values and words, the samples' strata."""
    pairs, options = [], []
    for i, r in enumerate(shown.iter_rows(named=True)):
        words = {"lncrna": r["lncrna_name"], "neighbour": r["neighbour_name"],
                 "kind": figures.kind_words(r["kind"]), "gap": f"{r['gap']:,}",
                 "q": f"{r['q']:.2g}", "q_within": f"{r['q_within']:.2g}"}  # fmt: skip
        pairs.append({"x": a[i].tolist(), "y": b[i].tolist(), "words": words})
        options.append(fill("figures.neighbour.option", **words))
    first = {**pairs[0]["words"], "rho": f"{shown['rho'][0]:.2f}",
             "rho_within": f"{shown['rho_within'][0]:.2f}"}  # fmt: skip
    template = text.raw(PAGE, "figures.neighbour.summary")
    return {
        "kind": "pairs", "options": options, "pairs": pairs, "strata": strata.tolist(),
        "members": members, "none": word("figures.neighbour.none"),
        "summary_template": template, "summary": template.format_map(first),
        "lead": {"template": text.raw(PAGE, "figures.neighbour.lead"), "fixed": {},
                 "start": first},
    }  # fmt: skip


def lncrna(writer: Writer, ctx: Context) -> Result:
    """F4 with its select, and the pairs by kind against the background."""
    inp, cols = ctx.inputs, contrast_columns(ctx)
    sub = inp.matrix[:, cols].astype(float)
    logcpm = np.log2(sub / sub.sum(axis=0) * 1e6 + 1)
    used = inp.used.sort(pl.col("group") != "reference", maintain_order=True)
    facets = [c for c in ("cell_type", "condition") if c in used.columns]
    strata = used.select(pl.concat_str(facets, separator="|"))[:, 0].to_numpy()
    pairs = figures.pair_table(inp.genes, inp.results, logcpm, strata)
    writer.table(pairs, "pairs", figures.columns(pairs))
    assert inp.spec.ncrna is not None  # noqa: S101 - resolve() refuses a file without [ncrna]
    shown = figures.shown_pairs(pairs, inp.spec.ncrna.pairs)
    row = {g: i for i, g in enumerate(inp.genes["gene_id"].to_list())}
    a = logcpm[[row[g] for g in shown["lncrna_id"]]]
    b = logcpm[[row[g] for g in shown["neighbour_id"]]]
    frame = used.select("sample", *facets).with_columns(
        lncrna=pl.Series(a[0]), neighbour=pl.Series(b[0])
    )
    colours = {inp.contrast.reference: ROLES["reference"], inp.contrast.tested: ROLES["comparison"]}
    members = [[i for i, c in enumerate(frame["condition"]) if c == level] for level in colours]
    live = pair_live(shown, a, b, strata, members)
    cis = pairs.filter(pl.col("kind") != "trans")
    caption = fill("figures.neighbour.caption", n_samples=frame.height, n_pairs=cis.height,
                   n_background=pairs.height - cis.height, contrast=inp.contrast.name)  # fmt: skip
    data = ScatterData(
        frame=frame, columns=figures.columns(frame), x="lncrna", y="neighbour",
        x_label=word("figures.neighbour.x"), y_label=word("figures.neighbour.y"),
        colour="condition", colours=colours, caption=caption, hover=("sample", *facets[:-1]),
        live=live,
    )  # fmt: skip
    reading = tuple(line.format(pairs=shown.height)
                    for line in text.raw(PAGE, "figures.neighbour.reading"))  # fmt: skip
    words = FigureText(title=word("figures.neighbour.title"),
                       lead=text.raw(PAGE, "figures.neighbour.lead").format_map(
                           live["lead"]["start"]), reading=reading)  # fmt: skip
    block = writer.figure(data, "lncrna_neighbour", words)
    table = ResultTable(fill("figures.neighbour.table_title", alpha=ALPHA),
                        text.words(PAGE, "figures.neighbour.columns"),
                        tuple(figures.kinds_table(pairs)), numeric=(1, 2, 3, 4))  # fmt: skip
    return Group(word("groups.lncrna"), (block, table))


def sanity(writer: Writer, ctx: Context, metrics: dict[str, float],
           reasons: dict[str, str]) -> Result | None:  # fmt: skip
    """M8's own rules decided on this run's metrics; the verdicts as a table and a group.

    A condition whose metric this run could not measure is NOT MEASURED, its reason printed
    where its value would be.
    """
    spec = ctx.inputs.spec.ncrna
    if spec is None or not spec.sanity:
        return None
    rows = decide(metrics, parse(spec.sanity), reasons)
    frame = pl.DataFrame(rows, schema={"rule": pl.String, "text": pl.String,
                                       "condition": pl.String, "left": pl.Float64,
                                       "right": pl.Float64, "verdict": pl.String,
                                       "source": pl.String, "rule_verdict": pl.String,
                                       "reason": pl.String})  # fmt: skip
    from ctx.de.sanity import VERDICT_COLUMNS  # noqa: PLC0415 - the one data dictionary

    writer.table(frame.select(list(VERDICT_COLUMNS)), "sanity", VERDICT_COLUMNS)
    cells = tuple((r["rule"], r["text"], r["condition"],
                   r["reason"] if r["left"] is None else f"{r['left']:.4g}", r["verdict"],
                   r["source"])
                  for r in rows)  # fmt: skip
    table = ResultTable("", text.words(PAGE, "sanity.columns"), cells, numeric=(3,))
    return Group(
        word("groups.sanity"),
        (table,),
        lead=fill("sanity.lead", registered=spec.sanity["registered_on"]),
    )


def keys(ctx: Context, frame: pl.DataFrame, top: float) -> tuple[KeyNumber, ...]:
    """Genes tested, lncRNAs tested, genes called within biotype, the largest miRNA share."""
    genes = ctx.inputs.genes
    in_groups = genes.filter(pl.col("biotype").is_in(list(ctx.biotypes))).height
    lnc = frame.filter(pl.col("biotype") == biotypes.LNCRNA).height
    all_lnc = genes.filter(pl.col("biotype") == biotypes.LNCRNA).height
    called = frame.filter(pl.col("call") != figures.NOT_CALLED).height
    return (
        KeyNumber(word("keys.tested"), f"{frame.height:,}",
                  fill("keys.tested_detail", n_genes=in_groups)),
        KeyNumber(word("keys.lncrna"), f"{lnc:,}",
                  fill("keys.lncrna_detail", pct=100 * lnc / all_lnc if all_lnc else 0.0)),
        KeyNumber(word("keys.called"), f"{called:,}",
                  fill("keys.called_detail", alpha=ALPHA, lfc_min=LFC_MIN)),
        KeyNumber(word("keys.mirna"), fill("keys.mirna_value", pct=top),
                  word("keys.mirna_detail")),
    )  # fmt: skip


def build(writer: Writer, ctx: Context) -> list[Path]:
    """Write every table and figure, then the page; the page first in the returned list."""
    shown, metrics = library(writer, ctx)
    tested, frame = de(writer, ctx)
    gdc = ctx.inputs.tcga is not None
    spec = ctx.inputs.spec.gdc
    lead = gdc_reader.acknowledgement(spec) if gdc and spec else (word("notes.tcga"),)
    results: list[Result] = [*lead, *shown, tested, lncrna(writer, ctx)]
    reasons: dict[str, str] = {}
    if gdc:
        found, measured, reasons = mirna_page.build(writer, ctx)
        results += found
        metrics |= measured
    metric_rows = [(k, float(v), detail) for k, (v, detail) in metrics.items()]
    writer.table(pl.DataFrame(metric_rows, schema={"metric": pl.String, "value": pl.Float64,
                                                   "detail": pl.String}, orient="row"),
                 "metrics", METRIC_COLUMNS)  # fmt: skip
    record = {k: float(v) for k, (v, _) in metrics.items()}
    if (verdicts := sanity(writer, ctx, record, reasons)) is not None:
        results.append(verdicts)
    usage, budget = runtime.measure(ctx.settings.started), runtime.budget(PAGE)
    skeleton = text.skeleton(PAGE)
    if gdc:
        skeleton |= {"input": word("tcga.input"), "method": word("tcga.method")}
    built = ModulePage(
        title=fill("title", dataset=ctx.inputs.spec.name), **skeleton,
        parameters=ctx.cards, results=tuple(results),
        key_numbers=keys(ctx, frame, record["mirna.max_percent"]),
        reproduce=Reproduction(
            command=ctx.command, wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self, peak_rss_children=usage.peak_rss_children,
            profile=f"{ctx.settings.profile}, {ctx.settings.threads} threads",
            budget=fill("notes.budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=tuple(p.name for p in ctx.inputs.read),
        ),
    )  # fmt: skip
    page = writer.out / f"{stem(PAGE, f'{writer.dataset}_overview', 1, day=writer.day)}.html"
    return [write(built, page), *writer.files]
