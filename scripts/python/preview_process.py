#!/usr/bin/env python3
# /// codebase-meta
# exempt: preview builder; draws the M3 page on made-up data until ctx process is built
# ///
"""The course preview of M3 Trim, align, count (``ctx process``), drawn on made-up data.

``ctx process`` is not built yet (TODO.md, M3). This page carries the figures its plan names:
the aligner's outcome per sample and per setting of one featured sample, reads assigned at each
featureCounts ``-s``, and plain against ``-M`` gene counts as a scatter and a Bland-Altman plot.
Every number comes from a seeded generator shaped like STAR's Log.final.out and featureCounts'
summaries; the statistics on those numbers are computed. Samples, groups and the plain gene
counts are the dataset's gather outputs, so this page agrees with the other preview pages. The
words are ``content/process.toml``, which the real module will reuse. ``build_preview.py``
calls ``page``.
"""

from __future__ import annotations

import time
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from scipy import stats

from ctx import runtime
from ctx.de.run import Writer, newest
from ctx.io.naming import MODULES, stem
from ctx.qc import Fastp
from ctx.qc import registry as reg
from ctx.report import text
from ctx.report.builder import (
    FigureBlock,
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    write,
)
from ctx.viz.bar import BarData
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import DATA_COLOURS, PROJECTOR_UNSAFE, ROLES, category_colours

# The made-up dataset's name on the page, one word for every preview (content/ui.toml).
SHOWN = ("ui", "preview.dataset")

PAGE = MODULE = "process"
SEED = 3
# The planned default (--aligner star); its outcome classes are STAR's Log.final.out rows.
ALIGNER = "STAR"
OUTCOMES = ("unique", "multi", "too_many", "mismatches", "too_short", "other")
# Six outcomes fit the seven data colours without the one a projector washes out.
OUTCOME_COLOURS = dict(zip(OUTCOMES, [c for c in DATA_COLOURS if c not in PROJECTOR_UNSAFE],
                           strict=False))  # fmt: skip
STRAND = ("s0", "s1", "s2")
# The M3 plan's statistics (section 11): percentile bootstrap, 1000 resamples, seed 42.
CI, RESAMPLES, BOOT_SEED = 0.95, 1000, 42
Z = 1.96  # limits of agreement: bias plus or minus Z standard deviations (Bland-Altman)
EMPIRICAL = (2.5, 97.5)
STRAND_METRIC = "salmon_strand_bias"  # nf-core/rnaseq's thresholds, until M3 adds its own
# The planned fastp options take their defaults from M2's Fastp, as the real CLI will.
FASTP = {"--fastp-q": "q", "--fastp-u": "u", "--fastp-l": "l", "--cut-tail": "cut_tail",
         "--cut-window": "window", "--cut-mean": "mean"}  # fmt: skip


def made_up_runs(
    rng: np.random.Generator, sheet: pl.DataFrame, assigned: np.ndarray
) -> pl.DataFrame:
    """Made-up Log.final.out shares and featureCounts strand shares, one row per sample.

    Each batch shifts the unique share; every library is reverse-stranded. A sample's gather
    counts are taken as its reads assigned at ``-s 2``, which fixes its input reads.
    """
    n = sheet.height
    levels = sheet["batch"].unique(maintain_order=True).to_list()
    shift = dict(zip(levels, rng.normal(0.0, 1.5, len(levels)), strict=True))
    batch = np.array([shift[b] for b in sheet["batch"]])
    shares = {
        "multi": np.clip(rng.normal(5.0, 0.6, n) - batch / 2, 1.0, None),
        "too_many": rng.uniform(0.1, 0.4, n),
        "mismatches": rng.uniform(0.05, 0.3, n),
        "too_short": np.clip(rng.normal(5.0, 1.2, n) - batch / 2, 1.0, None),
        "other": rng.uniform(0.3, 1.0, n),
    }
    unique = 100.0 - sum(shares.values())
    genic = unique / 100 * rng.uniform(0.74, 0.82, n)  # unique reads that fall on a gene
    forward = rng.uniform(0.04, 0.08, n)
    strand = {"s0": genic * (1 - rng.uniform(0.01, 0.02, n)), "s1": genic * forward,
              "s2": genic * (1 - forward)}  # fmt: skip
    return pl.DataFrame({
        "sample": sheet["sample"],
        "group": sheet.select(pl.concat_str("cell_type", "condition", "batch", separator=", "))
        .to_series(),
        "input_reads": np.round(assigned / strand["s2"]).astype(np.int64),
        "unique": unique, **shares, **strand, "forward_share": forward,
    })  # fmt: skip


def with_multimappers(rng: np.random.Generator, plain: np.ndarray) -> np.ndarray:
    """Made-up ``-M`` counts: the plain counts plus multimapped reads, most on gene families."""
    genes = plain.shape[0]
    # Per gene, extra reads per plain read: none, a few, or a gene family's many.
    rates = np.stack([np.zeros(genes), rng.lognormal(np.log(0.01), 0.8, genes),
                      rng.lognormal(np.log(0.3), 0.8, genes)])  # fmt: skip
    kind = rng.choice(len(rates), size=genes, p=(0.45, 0.45, 0.10))
    extra: np.ndarray = rng.poisson(plain * rates[kind, np.arange(genes)][:, None])
    return np.asarray(plain + extra)


def agreement(a: np.ndarray, b: np.ndarray) -> dict[str, Any]:
    """Correlation and Bland-Altman agreement of two counts, over genes either one counts."""
    kept = (a + b) > 0
    a, b = a[kept], b[kept]
    la, lb = np.log2(a + 1), np.log2(b + 1)
    d = lb - la
    n, bias, sd = a.size, float(d.mean()), float(d.std(ddof=1))
    half = float(stats.t.ppf((1 + CI) / 2, n - 1)) * sd / np.sqrt(n)
    pearson = stats.pearsonr(la, lb)
    r_low, r_high = pearson.confidence_interval(confidence_level=CI)
    d_low, d_high = np.percentile(d, EMPIRICAL)
    return {
        "n_genes": n, "n_identical": int((a == b).sum()),
        "spearman": float(stats.spearmanr(a, b).statistic), "spearman_low": None,
        "spearman_high": None, "pearson": float(pearson.statistic), "pearson_low": float(r_low),
        "pearson_high": float(r_high), "bias": bias, "bias_low": bias - half,
        "bias_high": bias + half, "sd": sd, "loa_low": bias - Z * sd, "loa_high": bias + Z * sd,
        "d_low": float(d_low), "d_high": float(d_high), "fold": 2**bias,
    }  # fmt: skip


def concordance(plain: np.ndarray, multi: np.ndarray, samples: list[str]) -> pl.DataFrame:
    """The gene sums over samples, with Spearman's bootstrap interval; then each sample."""
    a, b = plain.sum(axis=1), multi.sum(axis=1)
    kept = (a + b) > 0
    boot = stats.bootstrap((a[kept], b[kept]), lambda x, y: stats.spearmanr(x, y).statistic,
                           paired=True, vectorized=False, n_resamples=RESAMPLES,
                           confidence_level=CI, method="percentile",
                           rng=np.random.default_rng(BOOT_SEED))  # fmt: skip
    ci = boot.confidence_interval
    rows = [{"unit": "all", **agreement(a, b), "spearman_low": ci.low, "spearman_high": ci.high}]
    rows += [{"unit": s, **agreement(plain[:, i], multi[:, i])} for i, s in enumerate(samples)]
    return pl.DataFrame(rows).with_columns(comparison=pl.lit("plain_vs_M"))


def _columns(key: str, frame: pl.DataFrame) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column, from a content table with ``units`` and ``columns``."""
    entry = text.raw(PAGE, key)
    return {c: (entry["units"].get(c, "none"), entry["columns"][c]) for c in frame.columns}


def _text(key: str, **values: Any) -> FigureText:
    """Title, lead and reading guide of ``figures.<key>``, the lead filled from values."""
    return FigureText(title=text.word(PAGE, f"figures.{key}.title"),
                      lead=text.fill(PAGE, f"figures.{key}.lead", **values),
                      reading=text.words(PAGE, f"figures.{key}.reading"))  # fmt: skip


def _bars(key: str, frame: pl.DataFrame, caption: str, **layout: Any) -> BarData:
    """Bars of ``figures.<key>``: category first column, value ``percent``, series second."""
    category, series = frame.columns[:2]
    return BarData(
        frame=frame, columns=_columns(f"figures.{key}", frame), category=category,
        value="percent", category_label=text.word(PAGE, f"figures.{key}.category"),
        value_label=text.word(PAGE, f"figures.{key}.value"), series=series, caption=caption,
        **layout,
    )  # fmt: skip


def mapping_figures(writer: Writer, runs: pl.DataFrame) -> list[FigureBlock]:
    """The outcome bars of every sample, then of the featured sample under each setting."""
    names = {o: text.word(PAGE, f"outcomes.{o}") for o in OUTCOMES}
    long = runs.unpivot(OUTCOMES, index=["sample", "group", "input_reads"],
                        variable_name="outcome", value_name="percent").with_columns(
        reads=(pl.col("percent") / 100 * pl.col("input_reads")).round().cast(pl.Int64))  # fmt: skip
    long = long.select("sample", "outcome", "percent", "reads", "input_reads", "group")
    reads = int(runs["input_reads"].median())  # type: ignore[arg-type]
    low = runs.sort("unique").row(0, named=True)
    caption = text.fill(PAGE, "figures.mapping.caption", aligner=ALIGNER, n=runs.height,
                        reads=reads)  # fmt: skip
    bars = _bars("mapping", long, caption, colours=OUTCOME_COLOURS, names=names, stacked=True,
                 hover=("reads", "input_reads", "group"), width="double")  # fmt: skip
    lead = {"median": float(runs["unique"].median()), "low": low["unique"],  # type: ignore[arg-type]
            "low_sample": low["sample"]}  # fmt: skip
    blocks = [writer.figure(bars, "mapping", _text("mapping", **lead))]
    featured = runs.row(0, named=True)
    shares = settings_shares(featured)
    sample = featured["sample"]
    frame = pl.DataFrame([
        {"setting": text.fill(PAGE, f"settings.{s}", sample=sample), "outcome": o,
         "percent": p, "reads": round(p / 100 * featured["input_reads"])}
        for s, row in shares.items() for o, p in row.items()
    ])  # fmt: skip
    caption = text.fill(PAGE, "figures.settings.caption", sample=sample,
                        reads=featured["input_reads"])  # fmt: skip
    bars = _bars("settings", frame, caption, colours=OUTCOME_COLOURS, names=names, stacked=True,
                 horizontal=True, hover=("reads",), width="double")  # fmt: skip
    lead = {"sample": sample, "main": shares["main"]["unique"], "cut": shares["cut"]["unique"]}
    blocks.append(writer.figure(bars, "settings", _text("settings", **lead)))
    return blocks


def settings_shares(run: dict[str, Any]) -> dict[str, dict[str, float]]:
    """Made-up outcome shares of one sample as run, on half the genome, untrimmed, and at limit 1.

    Limit 1 is derived as the plan derives it: multimapped reads become too many loci.
    """
    main = {o: float(run[o]) for o in OUTCOMES}
    # About half the reads come from genes on the removed half; a few align to paralogues.
    cut = main | {"unique": main["unique"] * 0.52, "multi": main["multi"] * 0.6,
                  "mismatches": main["mismatches"] + 0.4, "other": main["other"] + 0.5}  # fmt: skip
    cut["too_short"] = 100.0 - sum(v for k, v in cut.items() if k != "too_short")
    # Adapter tails that fastp would cut leave reads too short to align, or mismatched.
    untrimmed = main | {"unique": main["unique"] - 1.4, "mismatches": main["mismatches"] + 0.3,
                        "too_short": main["too_short"] + 1.1}  # fmt: skip
    limit = main | {"multi": 0.0, "too_many": main["too_many"] + main["multi"]}
    return {"main": main, "cut": cut, "untrimmed": untrimmed, "limit": limit}


def strand_figure(writer: Writer, runs: pl.DataFrame) -> tuple[FigureBlock, str, int]:
    """Reads assigned at each ``-s`` per sample; also the majority call and how many share it."""
    metric = reg.load(reg.resolve(None)).metrics[STRAND_METRIC]
    calls = runs.with_columns(call=pl.col("forward_share").map_elements(
        metric.strandedness, return_dtype=pl.String))  # fmt: skip
    frame = calls.unpivot(STRAND, index=["sample", "input_reads", "forward_share", "call"],
                          variable_name="setting", value_name="share").with_columns(
        percent=pl.col("share") * 100,
        reads=(pl.col("share") * pl.col("input_reads")).round().cast(pl.Int64),
    ).select("sample", "setting", "percent", "reads", "forward_share", "call")  # fmt: skip
    call, n_call = Counter(calls["call"]).most_common(1)[0]
    share = runs["forward_share"]
    caption = text.fill(PAGE, "figures.strandedness.caption", n=runs.height,
                        share=float(share.median()),  # type: ignore[arg-type]
                        stranded=metric.stranded, unstranded=metric.unstranded)  # fmt: skip
    bars = _bars("strandedness", frame, caption, colours=category_colours(list(STRAND)),
                 names=text.raw(PAGE, "figures.strandedness.names"),
                 hover=("reads", "forward_share", "call"), width="double")  # fmt: skip
    lead = _text("strandedness", low=float(share.min()), high=float(share.max()),  # type: ignore[arg-type]
                 n_call=n_call, n=runs.height, call=call)  # fmt: skip
    return writer.figure(bars, "strandedness", lead), call, n_call


def counter_figures(writer: Writer, genes: pl.Series, plain: np.ndarray, multi: np.ndarray,
                    table: pl.DataFrame) -> list[FigureBlock]:  # fmt: skip
    """Plain against ``-M`` gene sums as a scatter, then as a Bland-Altman agreement plot."""
    a, b = plain.sum(axis=1), multi.sum(axis=1)
    frame = pl.DataFrame({"gene_id": genes, "plain": a, "multi": b}).filter(
        pl.col("plain") + pl.col("multi") > 0).with_columns(
        log2_plain=(pl.col("plain") + 1).log(2), log2_multi=(pl.col("multi") + 1).log(2),
        keep=pl.col("plain") != pl.col("multi"))  # fmt: skip
    frame = frame.with_columns(agreement=pl.when("keep").then(pl.lit("differs"))
                               .otherwise(pl.lit("same")))  # fmt: skip
    colours = {"same": None, "differs": ROLES["comparison"]}
    names = {k: text.word(PAGE, f"figures.counters.{k}") for k in colours}
    every, per_sample = table.row(0, named=True), table.filter(pl.col("unit") != "all")["spearman"]
    values = {**every, "n": every["n_genes"], "samples": plain.shape[1], "ci": round(CI * 100),
              "rho": every["spearman"], "rho_low": every["spearman_low"],
              "rho_high": every["spearman_high"], "rho_med": per_sample.median(),
              "rho_min": per_sample.min(), "rho_max": per_sample.max(), "r": every["pearson"],
              "r_low": every["pearson_low"], "r_high": every["pearson_high"],
              "n_same": every["n_identical"],
              "n_diff": every["n_genes"] - every["n_identical"]}  # fmt: skip
    span = (float(frame["log2_plain"].min()), float(frame["log2_multi"].max()))  # type: ignore[arg-type]
    scatter = ScatterData(
        frame=frame, columns=_columns("figures.counters", frame), x="log2_plain",
        y="log2_multi", x_label=text.word(PAGE, "figures.counters.x"),
        y_label=text.word(PAGE, "figures.counters.y"), colour="agreement", colours=colours,
        caption=text.fill(PAGE, "figures.counters.caption", **values),
        hover=("gene_id", "plain", "multi"), keep="keep", names=names,
        lines=(Line(text.word(PAGE, "figures.counters.identity"), span, span),),
    )  # fmt: skip
    blocks = [writer.figure(scatter, "counters_multimap", _text("counters", **values))]
    blocks.append(writer.figure(*agreement_plot(frame, values, colours, names)))
    return blocks


def agreement_plot(frame: pl.DataFrame, values: dict[str, Any], colours: dict[str, str | None],
                   names: dict[str, str]) -> tuple[ScatterData, str, FigureText]:  # fmt: skip
    """The Bland-Altman plot of the counter scatter's genes, with its stem and words."""
    ba = frame.select("gene_id", mean_log2=(pl.col("log2_plain") + pl.col("log2_multi")) / 2,
                      diff_log2=pl.col("log2_multi") - pl.col("log2_plain"),
                      agreement="agreement", keep="keep")  # fmt: skip
    span = (float(ba["mean_log2"].min()), float(ba["mean_log2"].max()))  # type: ignore[arg-type]
    bias, limits = (text.word(PAGE, f"figures.agreement.{k}") for k in ("bias", "limits"))
    outside = ~ba["diff_log2"].is_between(values["loa_low"], values["loa_high"])
    data = ScatterData(
        frame=ba, columns=_columns("figures.agreement", ba), x="mean_log2", y="diff_log2",
        x_label=text.word(PAGE, "figures.agreement.x"),
        y_label=text.word(PAGE, "figures.agreement.y"), colour="agreement", colours=colours,
        caption=text.fill(PAGE, "figures.agreement.caption", **values), hover=("gene_id",),
        keep="keep", names=names,
        lines=(Line(bias, span, (values["bias"],) * 2),
               *(Line(limits, span, (values[k],) * 2, dash=True) for k in ("loa_low", "loa_high"))),
    )  # fmt: skip
    return data, "agreement_multimap", _text("agreement", **values, n_out=int(outside.sum()))


def cards() -> tuple[ParameterCard, ...]:
    """The options ``ctx process`` will take, in the plan's order, fastp's defaults from M2."""
    fastp = Fastp()
    shown = {True: "on", False: "off"}
    return tuple(
        ParameterCard(flag, card["meaning"],
                      card.get("default") or shown.get(v := getattr(fastp, FASTP[flag]), str(v)),
                      card["effect"])
        for flag, card in text.raw(PAGE, "planned").items()
    )  # fmt: skip


def page(root: Path, dataset: str, day: date) -> Path:
    """Write the M3 preview page and its figures under ``root/outputs/03_process``.

    Args:
        root: A ctx root holding the dataset's gather outputs.
        dataset: The dataset slug.
        day: The date every output name carries.

    Returns:
        The page.
    """
    started = time.perf_counter()
    gathered = root / "outputs" / "01_gather"
    inputs = [newest(gathered, dataset, what) for what in ("samplesheet", "counts")]
    sheet = pl.read_parquet(inputs[0]).unique("sample", keep="first", maintain_order=True)
    counts = pl.read_parquet(inputs[1])
    samples = sheet["sample"].to_list()
    plain = counts.select(samples).to_numpy().astype(np.int64)
    rng = np.random.default_rng(SEED)
    runs = made_up_runs(rng, sheet, plain.sum(axis=0))
    multi = with_multimappers(rng, plain)
    releases = {"input:samplesheet": inputs[0].name, "input:counts": inputs[1].name,
                "made_up": f"numpy default_rng({SEED})"}  # fmt: skip
    writer = Writer(root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}", dataset, day,
                    releases)  # fmt: skip
    table = concordance(plain, multi, samples)
    writer.table(table, "concordance", _columns("tables.concordance", table))
    strand, call, n_call = strand_figure(writer, runs)
    counters = counter_figures(writer, counts["gene_id"], plain, multi, table)
    every = table.row(0, named=True)
    results = (
        text.word("ui", "preview.note"),
        Group(text.word(PAGE, "groups.align"), tuple(mapping_figures(writer, runs))),
        Group(text.word(PAGE, "groups.strand"), (strand,)),
        Group(text.word(PAGE, "groups.counters"), tuple(counters)),
    )
    keys = (
        KeyNumber(text.word(PAGE, "keys.samples"), f"{runs.height:,}"),
        KeyNumber(text.word(PAGE, "keys.unique"),
                  f"{float(runs['unique'].median()):.1f} %"),  # type: ignore[arg-type]
        KeyNumber(text.word(PAGE, "keys.strand"), call,
                  text.fill(PAGE, "keys.strand_detail", n_call=n_call, n=runs.height)),
        KeyNumber(text.word(PAGE, "keys.differ"), f"{every['n_genes'] - every['n_identical']:,}",
                  text.fill(PAGE, "keys.differ_detail", n=every["n_genes"])),
    )  # fmt: skip
    usage, budget = runtime.measure(started), runtime.budget(MODULE)
    command = f"ctx process --dataset {dataset} --dataset-file {dataset}.toml --genome genome.fa.gz"
    built = ModulePage(
        title=text.fill(PAGE, "title", dataset=text.word(*SHOWN)), **text.skeleton(PAGE),
        parameters=cards(), results=results, key_numbers=keys,
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            # page() starts no child process; the build's own ctx children are not this page's.
            peak_rss_children=0, profile=text.word(PAGE, "profile"),
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=tuple(p.name for p in inputs),
        ),
    )  # fmt: skip
    return write(built, writer.out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html")
