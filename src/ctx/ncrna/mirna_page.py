"""The TCGA half of the M8 page: the miRNA fit, the target shifts, the pairs, the two layers.

Each group writes its tables and figures through the run's ``Writer``; the numbers come from
``mirna`` and ``targets``, every sentence from ``content/ncrna.toml``. The metrics returned
are the ones ``[ncrna.sanity]`` rules read.
"""

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

import numpy as np
import polars as pl

from ctx.de.deseq import RESULT_COLUMNS
from ctx.de.run import Writer, volcano
from ctx.de.stats import ALPHA
from ctx.ncrna import mirna, targets
from ctx.ncrna.figures import columns as figure_columns
from ctx.report import text
from ctx.report.builder import FigureText, Group, Result, ResultTable
from ctx.viz.curve import CurveData
from ctx.viz.scatter import ScatterData, needed_points
from ctx.viz.style import ALONE_SCALE, LARGEST_PLOT_PX, PLOTLY_HEIGHT_PX, ROLES, THEMES

if TYPE_CHECKING:
    from ctx.ncrna.page import Context
    from ctx.ncrna.run import Tcga

PAGE = "ncrna"
Metrics = dict[str, tuple[float, str]]
CLASS_COLOURS = {targets.VALIDATED: ROLES["comparison"], targets.PREDICTED: ROLES["reference"],
                 targets.NON_TARGET: THEMES["white"].reference}  # fmt: skip
COLUMNS = {
    "mimat": ("none", "miRBase mature accession."),
    "mirna": ("none", "miRBase name of the mature miRNA, from TargetScan's family file."),
    "gene_id": ("none", "Ensembl gene ID without its version."),
    "gene_name": ("none", "Gene symbol."),
    "source": ("none", "validated (miRTarBase strong evidence) or predicted only (TargetScan)."),
    "target": ("none", "validated, predicted only or non-target."),
    "preset": ("context++", "Cumulative weighted context++ cutoff; empty for validated targets."),
    "n": ("genes", "Genes of the class M5 tested."),
    "n_non": ("genes", "Non-target genes M5 tested: no TargetScan site, not validated."),
    "shift": ("log2 fold change", "Median log2 FC of the class minus that of the non-targets."),
    "d": ("dimensionless", "Two-sample KS statistic against the non-targets."),
    "p_ks": ("dimensionless", "KS p-value."),
    "p_mw": ("dimensionless", "Mann-Whitney p-value."),
    "q_ks": ("dimensionless", "BH of p_ks over every row."),
    "q_mw": ("dimensionless", "BH of p_mw over every row."),
    "log2FoldChange": ("log2 fold change", "M5's raw log2 fold change, tested over reference."),
    "share": ("fraction", "Share of the class's genes at or below this fold change."),
    "genes": ("none", "The called genes tested: down for a rising miRNA, up for a falling one."),
    "size": ("genes", "Targets of the set among the background."),
    "hits": ("genes", "Targets among the called genes."),
    "expected": ("genes", "Hits expected by chance."),
    "p": ("dimensionless", "Hypergeometric p, P(X >= hits)."),
    "q": ("dimensionless", "BH of p over every row."),
    "score": ("context++", "TargetScan cumulative weighted context++ score; empty if validated."),
    "evidence": ("none", "miRTarBase experiments, or TargetScan's site count."),
    **{f"{k}_{c}": (u, m.format(c=c)) for c in mirna.CONDITIONS for k, u, m in (
        ("n", "samples", "Samples of condition {c} with both libraries."),
        ("rho", "dimensionless", "Spearman rho across the {c} samples."),
        ("p", "dimensionless", "Two-sided p of rho_{c}, t approximation."),
        ("q", "dimensionless", "BH of p_{c} over the pairs."))},
    "mirna_sample": ("none", "Aliquot of the small-RNA library."),
    "rna_sample": ("none", "Aliquot of the mRNA library of the same sample barcode."),
    "condition": ("none", "tumour or normal."),
    "rho": ("dimensionless", "Spearman rho of hairpin reads against STAR counts of their genes."),
    "hairpins": ("hairpins", "Hairpins matched one to one with a GENCODE miRNA gene."),
    "mirna_value": ("log2 RPM", "log2 of 1 plus reads per million miRNA-mapped reads."),
    "target_value": ("log2 CPM", "log2 of 1 plus counts per million counted reads."),
    "case": ("none", "Patient submitter ID."),
}  # fmt: skip


def word(key: str) -> str:
    """A fixed text of this page's content file."""
    return text.word(PAGE, key)


def fill(key: str, **values: Any) -> str:
    """A template of this page's content file, filled."""
    return text.fill(PAGE, key, **values)


def described(frame: pl.DataFrame) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of a table of this module."""
    return {c: COLUMNS.get(c) or figure_columns(frame.select(c))[c] for c in frame.columns}


def _number(value: float, spec: str) -> str:
    """A number as the page prints it; empty measurements print as none."""
    return format(value, spec) if value is not None and np.isfinite(value) else word(
        "figures.mirna_pairs.none")  # fmt: skip


def fitted(writer: Writer, ctx: "Context", tcga: "Tcga") -> tuple[Result, pl.DataFrame, Metrics]:
    """F6: the paired fit of the mature miRNAs, its volcano and the featured miRNAs' rows."""
    result = mirna.fit_mature(tcga.mature, tcga.sheet, tcga.contrast, ctx.settings.threads)
    named = mirna.names(tcga.mature["mimat"].to_list(), tcga.found.names)
    frame = result.results.join(named, on="gene_id", how="left")
    writer.table(frame, "mirna_results",
                 {**RESULT_COLUMNS, "gene_name": COLUMNS["mirna"]})  # fmt: skip
    tested = frame.filter(pl.col("padj").is_not_null())
    c = tcga.contrast
    shown = FigureText(title=fill("figures.mirna_volcano.title", tested=c.tested,
                                  reference=c.reference),
                       lead=fill("figures.mirna_volcano.lead", n_tested=tested.height,
                                 up=int((tested["call"] == "up").sum()),
                                 down=int((tested["call"] == "down").sum())))  # fmt: skip
    block = writer.figure(volcano(result, named), "mirna_volcano", shown)
    assert ctx.inputs.spec.ncrna is not None  # noqa: S101 - resolve() refuses a file without it
    metrics: Metrics = {
        "mirna.tested": (float(tested.height), f"mature miRNAs with an adjusted p in {c.name}"),
        "mirna.up": (float((tested["call"] == "up").sum()), f"called up in {c.name}"),
        "mirna.down": (float((tested["call"] == "down").sum()), f"called down in {c.name}"),
    }
    rows = []
    for name in ctx.inputs.spec.ncrna.mirnas:
        hit = tested.filter(pl.col("gene_name") == name)
        if hit.is_empty():
            rows.append((name, "", "", "", "", word("figures.mirna_volcano.untested")))
            continue
        r = hit.row(0, named=True)
        metrics[f"mirna.{name}.lfc"] = (r["log2FoldChange"], f"log2 FC, {c.name}")
        metrics[f"mirna.{name}.padj"] = (r["padj"], f"BH-adjusted p, {c.name}")
        rows.append((name, r["gene_id"], f"{r['baseMean']:,.0f}", f"{r['log2FoldChange']:+.2f}",
                     f"{r['padj']:.2g}", r["call"]))  # fmt: skip
    table = ResultTable(word("figures.mirna_volcano.table_title"),
                        text.words(PAGE, "figures.mirna_volcano.columns"), tuple(rows),
                        numeric=(2, 3, 4))  # fmt: skip
    return Group(word("groups.mirna"), (block, table)), frame, metrics


def _tested(ctx: "Context") -> pl.DataFrame:
    """M5's tested genes of the featured contrast, with the columns the target tests read."""
    return ctx.inputs.results.filter(pl.col("padj").is_not_null()).select(
        "gene_id", "gene_name", "log2FoldChange", "call"
    )


def _shift_rows(shifts: pl.DataFrame) -> tuple[tuple[str, ...], ...]:
    """The tests table as the page prints it."""
    return tuple(
        (r["mirna"], r["target"],
         word("figures.shift.any") if r["preset"] is None else f"{r['preset']:g}",
         f"{int(r['n']):,}", _number(r["shift"], "+.3f"), _number(r["d"], ".3f"),
         _number(r["q_ks"], ".2g"), _number(r["q_mw"], ".2g"))
        for r in shifts.iter_rows(named=True)
    )  # fmt: skip


def shift(writer: Writer, ctx: "Context", tcga: "Tcga",
          fit: pl.DataFrame) -> tuple[Result, Metrics]:  # fmt: skip
    """F7 for the first featured miRNA, every test as a table, and the enrichment.

    ``fit`` is the miRNA fit's results with names, whose fold changes set each miRNA's
    expected direction in the enrichment.
    """
    assert ctx.inputs.spec.ncrna is not None  # noqa: S101 - resolve() refuses a file without it
    mirnas = [m for m in ctx.inputs.spec.ncrna.mirnas if m in set(tcga.found.names.values())]
    tested = _tested(ctx)
    shifts = mirna.shifts(tcga.found, mirnas, tested)
    writer.table(shifts, "target_shift", described(shifts))
    first = mirnas[0]
    curves = _thinned(mirna.curves(tcga.found, first, tested))
    row = {(r["target"], r["preset"]): r for r in shifts.filter(pl.col("mirna") == first)
           .iter_rows(named=True)}  # fmt: skip
    valid = row[(targets.VALIDATED, None)]
    pred = row[(targets.PREDICTED, targets.DEFAULT_PRESET)]
    c = ctx.inputs.contrast
    caption = fill("figures.shift.caption", n_validated=int(valid["n"]),
                   n_predicted=int(pred["n"]), n_non=int(valid["n_non"]), n_tests=shifts.height,
                   preset=targets.DEFAULT_PRESET,
                   d_validated=_number(valid["d"], ".3f"),
                   q_validated=_number(valid["q_ks"], ".2g"),
                   d_predicted=_number(pred["d"], ".3f"),
                   q_predicted=_number(pred["q_ks"], ".2g"))  # fmt: skip
    data = CurveData(
        frame=curves, columns=described(curves), series="target", x="log2FoldChange", y="share",
        x_label=fill("figures.shift.x", tested=c.tested, reference=c.reference),
        y_label=word("figures.shift.y"), colour="target",
        colours={k: v for k, v in CLASS_COLOURS.items() if k in set(curves["target"])},
        caption=caption,
    )  # fmt: skip
    reading = tuple(line.format(preset=targets.DEFAULT_PRESET)
                    for line in text.raw(PAGE, "figures.shift.reading"))  # fmt: skip
    shown = FigureText(title=fill("figures.shift.title", mirna=first),
                       lead=fill("figures.shift.lead",
                                 shift_validated=_number(valid["shift"], "+.3f"),
                                 shift_predicted=_number(pred["shift"], "+.3f")),
                       reading=reading)  # fmt: skip
    block = writer.figure(data, "target_shift", shown)
    tests = ResultTable(word("figures.shift.table_title"),
                        text.words(PAGE, "figures.shift.columns"), _shift_rows(shifts),
                        numeric=(3, 4, 5, 6, 7))  # fmt: skip
    found = dict(
        fit.filter(pl.col("padj").is_not_null()).select("gene_name", "log2FoldChange").rows()
    )
    lfc = {m: float(found.get(m, np.nan)) for m in mirnas}
    enriched = mirna.enrichment(tcga.found, lfc, tested)
    writer.table(enriched, "target_enrichment", described(enriched))
    rows = tuple((r["mirna"], r["target"], r["genes"], f"{r['size']:,}", f"{r['hits']:,}",
                  f"{r['expected']:.1f}", f"{r['p']:.2g}", f"{r['q']:.2g}")
                 for r in enriched.iter_rows(named=True))  # fmt: skip
    enrich = ResultTable(word("figures.shift.enrich_title"),
                         text.words(PAGE, "figures.shift.enrich_columns"), rows,
                         numeric=(3, 4, 5, 6, 7))  # fmt: skip
    kept = targets.derived(tcga.found, mirnas, ctx.inputs.genes)
    writer.table(kept, "targets", described(kept))
    metrics: Metrics = {f"shift.{r['mirna']}.{r['target']}.{r['preset']}": (
        r["shift"], f"median log2 FC shift; KS q {r['q_ks']:.2g}")
        for r in shifts.iter_rows(named=True)}  # fmt: skip
    return Group(word("groups.targets"), (block, tests, enrich)), metrics


def _thinned(curves: pl.DataFrame) -> pl.DataFrame:
    """Each class's curve kept to the points that draw it within a pixel at the largest plot."""
    parts = []
    for target in curves["target"].unique(maintain_order=True):
        part = curves.filter(pl.col("target") == target)
        x, y = part["log2FoldChange"].to_numpy(), part["share"].to_numpy()
        keep = needed_points(x, y, (LARGEST_PLOT_PX, PLOTLY_HEIGHT_PX["scatter"] * ALONE_SCALE))
        parts.append(part[np.asarray(keep).tolist()])
    return pl.DataFrame(pl.concat(parts))


def _pair_words(r: Mapping[str, Any]) -> dict[str, str]:
    """One pair's words for the select, the lead and the summary."""
    return {"mirna": r["mirna"], "gene": r["gene_name"], "source": r["source"],
            "q_tumour": _number(r["q_tumour"], ".2g"),
            "rho_normal": _number(r["rho_normal"], ".2f")}  # fmt: skip


def pairs(writer: Writer, ctx: "Context", tcga: "Tcga",
          both: mirna.Layers) -> tuple[Result, Metrics]:  # fmt: skip
    """F8: each featured miRNA against its targets across tumours, normals apart, with a select."""
    assert ctx.inputs.spec.ncrna is not None  # noqa: S101 - resolve() refuses a file without it
    by_name = {v: k for k, v in tcga.found.names.items()}
    wanted = mirna.target_pairs(tcga.found, ctx.inputs.spec.ncrna.mirnas, by_name, _tested(ctx))
    table = mirna.pair_table(both, wanted)
    writer.table(table, "mirna_target_pairs", described(table))
    shown = (table.filter(pl.col("q_tumour").is_finite())
             .sort("q_tumour", -pl.col("rho_tumour").abs(), "mirna", "gene_name")
             .head(ctx.inputs.spec.ncrna.pairs))  # fmt: skip
    tumours = (both.samples["condition"] == "tumour").to_numpy()
    mir_row = {m: i for i, m in enumerate(both.mimats)}
    gene_row = {g: i for i, g in enumerate(both.genes)}
    xs = both.mirna[[mir_row[m] for m in shown["mimat"]]][:, tumours]
    ys = both.rna[[gene_row[g] for g in shown["gene_id"]]][:, tumours]
    frame = both.samples.filter(pl.Series(tumours)).select("mirna_sample", "case", "condition")
    frame = frame.with_columns(mirna_value=pl.Series(xs[0]), target_value=pl.Series(ys[0]))
    options, live_pairs = [], []
    for i, r in enumerate(shown.iter_rows(named=True)):
        words = _pair_words(r)
        options.append(fill("figures.mirna_pairs.option", **words))
        live_pairs.append({"x": xs[i].tolist(), "y": ys[i].tolist(), "words": words})
    first = {**live_pairs[0]["words"], "rho": f"{shown['rho_tumour'][0]:.2f}"}
    template = text.raw(PAGE, "figures.mirna_pairs.summary")
    live = {"kind": "pairs", "label": word("figures.mirna_pairs.label"),
            "options": options, "pairs": live_pairs,
            "strata": ["tumour"] * frame.height, "members": [list(range(frame.height))],
            "none": word("figures.mirna_pairs.none"), "summary_template": template,
            "summary": template.format_map(first),
            "lead": {"template": text.raw(PAGE, "figures.mirna_pairs.lead"), "fixed": {},
                     "start": first}}  # fmt: skip
    counts = {c: int((both.samples["condition"] == c).sum()) for c in mirna.CONDITIONS}
    caption = fill("figures.mirna_pairs.caption", n_tumour=counts["tumour"],
                   n_normal=counts["normal"], n_pairs=table.height)  # fmt: skip
    data = ScatterData(
        frame=frame, columns=described(frame), x="mirna_value", y="target_value",
        x_label=word("figures.mirna_pairs.x"), y_label=word("figures.mirna_pairs.y"),
        colour="condition", colours={"tumour": ROLES["comparison"]}, caption=caption,
        hover=("mirna_sample", "case"), live=live,
    )  # fmt: skip
    reading = tuple(line.format(pairs=shown.height)
                    for line in text.raw(PAGE, "figures.mirna_pairs.reading"))  # fmt: skip
    shown_text = FigureText(title=word("figures.mirna_pairs.title"),
                            lead=text.raw(PAGE, "figures.mirna_pairs.lead").format_map(first),
                            reading=reading)  # fmt: skip
    block = writer.figure(data, "mirna_target_pairs", shown_text)
    summary = ResultTable(fill("figures.mirna_pairs.table_title", alpha=ALPHA),
                          text.words(PAGE, "figures.mirna_pairs.columns"), _pair_rows(table),
                          numeric=(2, 3, 4, 5))  # fmt: skip
    metrics: Metrics = {}
    for r in table.iter_rows(named=True):
        for key in ("rho_tumour", "p_tumour", "rho_normal"):
            metrics[f"pair.{r['mirna']}.{r['gene_name']}.{key}"] = (
                r[key],
                f"{r['source']}, {r['n_tumour']} tumours, {r['n_normal']} normals",
            )
    return Group(word("groups.pairs"), (block, summary)), metrics


def _pair_rows(table: pl.DataFrame) -> tuple[tuple[str, ...], ...]:
    """Per miRNA and source: pairs, median rho in each condition, the share called in tumours."""
    rows = []
    for (name, source), part in table.group_by("mirna", "source", maintain_order=True):
        called = 100 * float(np.mean(part["q_tumour"].fill_nan(1.0).to_numpy() < ALPHA))
        rows.append((str(name), str(source), f"{part.height:,}",
                     _number(float(np.nanmedian(part["rho_tumour"].to_numpy())), ".2f"),
                     f"{called:.0f} %",
                     _number(float(np.nanmedian(part["rho_normal"].to_numpy())),
                             ".2f")))  # fmt: skip
    return tuple(sorted(rows))


def agreement(writer: Writer, tcga: "Tcga", both: mirna.Layers) -> tuple[str, Metrics]:
    """Per sample, the STAR miRNA genes against the small-RNA library; the median rho."""
    matched = mirna.placed(tcga.loci, tcga.types)
    frame = mirna.agreement(both, tcga.hairpins, matched, tcga.counts)
    writer.table(frame, "layer_agreement", described(frame))
    rho = float(np.nanmedian(frame["rho"].to_numpy())) if frame.height else float("nan")
    hairpins = int(frame["hairpins"][0]) if frame.height else 0
    detail = f"median over {frame.height} samples, {hairpins} hairpins placed one to one"
    note = fill("tcga.agreement", hairpins=hairpins, rho=rho, n=frame.height)
    return note, {"agreement.median_rho": (rho, detail)}


def build(writer: Writer, ctx: "Context") -> tuple[list[Result], Metrics]:
    """Every TCGA group, in page order, and the metrics they measured."""
    tcga = ctx.inputs.tcga
    assert tcga is not None  # noqa: S101 - page.build calls this for a GDC dataset alone
    fit_group, fit, metrics = fitted(writer, ctx, tcga)
    both = mirna.layers(tcga.sheet, tcga.mature, tcga.counts)
    shift_group, found = shift(writer, ctx, tcga, fit)
    pair_group, paired = pairs(writer, ctx, tcga, both)
    note, agreed = agreement(writer, tcga, both)
    return [fit_group, note, shift_group, pair_group], metrics | found | paired | agreed


def missing(record: Mapping[str, float], raw: Mapping[str, Any]) -> str | None:
    """A note naming every metric a rule reads that this run did not measure, or None."""
    wanted = {c.get(k) for r in raw.get("rules", []) for c in r.get("conditions", [])
              for k in ("metric", "other")} - {None}  # fmt: skip
    absent = sorted(m for m in wanted if m not in record or not np.isfinite(record[m]))
    return fill("tcga.missing", metrics=", ".join(absent)) if absent else None
