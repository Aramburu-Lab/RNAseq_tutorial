#!/usr/bin/env python3
# /// codebase-meta
# exempt: preview builder; draws the M9 Gene cards page on made-up data until ctx genecard exists
# ///
r"""The M9 Gene cards preview page, on made-up data (scripts/python/build_preview.py).

Every value comes from one seeded numpy Generator: a made-up gene in the dataset's samples,
across made-up GTEx tissues, HPA immune cell types and exons. Only sample names, cell types,
conditions and batches come from the made-up dataset's gather sheet, so the preview pages show
the same samples. The text is ``content/genecard_preview.toml``; tau is ``ctx genecard``'s own.

Usage:
    pixi run python scripts/python/build_preview.py --root <dir> --only genecard
"""

from __future__ import annotations

from dataclasses import replace
from datetime import date
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import polars as pl

from ctx import runtime
from ctx.de.run import Writer
from ctx.genecard.summary import tau
from ctx.io.naming import MODULES, stem
from ctx.report import text
from ctx.report.builder import (
    FigureBlock,
    FigureText,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    ResultTable,
    write,
)
from ctx.report.entry import write_entry
from ctx.viz.annotate import Caption, wrap
from ctx.viz.bar import BarData
from ctx.viz.box import box_data
from ctx.viz.heatmap import HeatmapData, diverging_levels
from ctx.viz.strip import Range, StripData
from ctx.viz.style import DATA_COLOURS, PROJECTOR_UNSAFE, category_colours

# The made-up dataset's name on the page, one word for every preview (content/ui.toml).
SHOWN = ("ui", "preview.dataset")

MODULE, PAGE = "genecard", "genecard_preview"
SEED = 9
# Not an HGNC symbol: no real gene is drawn with made-up values.
GENE = "MADEUP1"
# Our evidence, log2 CPM: first cell type's mean, the step down per later cell type, the case
# shift planted in the first cell type alone, and the sample-to-sample SD.
BASE, CELL_STEP, SHIFT, SAMPLE_SD = 5.0, -2.0, 1.8, 0.5
# Made-up GTEx: median TPM, and the use of the alternative exon, per tissue.
TISSUES = {"Whole blood": (85.0, 1.0), "Spleen": (40.0, 0.9), "Lung": (22.0, 0.7),
           "Liver": (6.0, 0.5), "Kidney cortex": (5.0, 0.6), "Skin": (4.0, 0.5),
           "Brain cortex": (0.8, 0.15), "Skeletal muscle": (0.5, 0.2)}  # fmt: skip
DONORS, MALE_SHARE, PER_TISSUE, DONOR_SD = 400, 0.66, (60, 160), 0.6
# Made-up HPA immune cells: cell type, lineage, pooled nTPM. The 18 names are HPA's as recalled,
# not yet checked against its files; the external-data spike reads them.
IMMUNE = (
    ("neutrophil", "granulocytes", 420.0), ("eosinophil", "granulocytes", 150.0),
    ("basophil", "granulocytes", 90.0), ("classical monocyte", "monocytes", 260.0),
    ("intermediate monocyte", "monocytes", 180.0), ("non-classical monocyte", "monocytes", 120.0),
    ("myeloid DC", "dendritic cells", 70.0), ("plasmacytoid DC", "dendritic cells", 35.0),
    ("naive B-cell", "B cells", 12.0), ("memory B-cell", "B cells", 10.0),
    ("naive CD4 T-cell", "T cells", 6.0), ("memory CD4 T-cell", "T cells", 8.0),
    ("T-reg", "T cells", 7.0), ("naive CD8 T-cell", "T cells", 5.0),
    ("memory CD8 T-cell", "T cells", 9.0), ("MAIT T-cell", "T cells", 11.0),
    ("gdT-cell", "T cells", 10.0), ("NK-cell", "NK cells", 14.0),
)  # fmt: skip
IMMUNE_SD = 0.15
# Six lineages, six colours: every data colour but the yellow a projector washes out.
LINEAGE_COLOURS = tuple(c for c in DATA_COLOURS if c not in PROJECTOR_UNSAFE)
# Made-up collapsed gene model: exon lengths in nt, the alternative exon's index, the tissue
# whose degraded RNA raises exons toward the poly(A) tail, and that rise per exon.
EXON_NT = (120, 180, 150, 90, 200, 160, 140, 170, 600)
ALT_EXON, DEGRADED, BIAS_STEP, EXON_SD = 3, "Liver", 0.12, 0.05


class MadeUpCaption(Caption):
    """The box caption with the content file's made-up note, in the PDF as on the page."""

    def text(self, width: int = 0) -> str:
        """The note and the test sentence, wrapped to ``width``; 0 leaves one line."""
        return wrap(text.fill(PAGE, "figures.ours.caption", stats=super().text()), width)


def _columns(key: str) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of one figure's table, from the content file."""
    return {c: (u, m) for c, (u, m) in text.raw(PAGE, f"figures.{key}.columns").items()}


def _words(key: str, **values: Any) -> FigureText:
    """Title, lead and reading guide of one figure, filled from its computed numbers."""
    return FigureText(title=text.fill(PAGE, f"figures.{key}.title", gene=GENE, **values),
                      lead=text.fill(PAGE, f"figures.{key}.lead", **values),
                      reading=text.words(PAGE, f"figures.{key}.reading"))  # fmt: skip


def ours(writer: Writer, rng: np.random.Generator, sheet: pl.DataFrame) -> FigureBlock:
    """The gene in the dataset's own samples: box per cell type and condition, Mann-Whitney."""
    cells = sheet["cell_type"].unique(maintain_order=True).to_list()
    reference, comparison = sheet["condition"].unique(maintain_order=True).to_list()
    index = sheet["cell_type"].replace_strict({c: i for i, c in enumerate(cells)}).to_numpy()
    case = (sheet["condition"] == comparison).to_numpy()
    value = BASE + CELL_STEP * index + SHIFT * (case & (index == 0))
    frame = sheet.select(feature="cell_type", group="condition", sample="sample", batch="batch")
    frame = frame.with_columns(value=value + rng.normal(0.0, SAMPLE_SD, sheet.height))
    box = box_data(frame, columns=_columns("ours"), feature_label=text.word(PAGE, "figures.ours.x"),
                   value_label=text.word(PAGE, "figures.ours.y"), groups=(reference, comparison),
                   hover=("sample", "batch"))  # fmt: skip
    box = replace(box, caption=MadeUpCaption(box.caption.groups, box.caption.tests))
    medians = (frame.group_by("feature", "group").agg(pl.col("value").median())
               .pivot("group", index="feature", values="value")
               .with_columns(shift=pl.col(comparison) - pl.col(reference)))  # fmt: skip
    top = medians.sort(pl.col("shift").abs(), descending=True).row(0, named=True)
    words = _words("ours", comparison=comparison, reference=reference, shift=top["shift"],
                   cell_type=top["feature"])  # fmt: skip
    return writer.figure(box, "own_evidence", words)


def tissues(writer: Writer, rng: np.random.Generator) -> FigureBlock:
    """Made-up GTEx per-sample TPM by tissue, with each tissue's median and IQR; no test."""
    male = rng.random(DONORS) < MALE_SHARE
    rows = []
    for tissue, (median, _) in TISSUES.items():
        donors = rng.choice(DONORS, int(rng.integers(*PER_TISSUE)), replace=False)
        tpm = median * np.exp(rng.normal(0.0, DONOR_SD, donors.size))
        rows += [(tissue, f"D{d:04d}", "male" if male[d] else "female", float(t))
                 for d, t in zip(donors, tpm, strict=True)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=list(_columns("tissues")), orient="row")
    stats = frame.group_by("tissue", maintain_order=True).agg(
        n=pl.len(), median=pl.col("tpm").median(), q1=pl.col("tpm").quantile(0.25, "linear"),
        q3=pl.col("tpm").quantile(0.75, "linear"))  # fmt: skip
    ranges = tuple(Range(t, m, lo, hi) for t, _, m, lo, hi in stats.rows())
    by_median = stats.sort("median", descending=True)
    caption = text.fill(PAGE, "figures.tissues.caption", samples=frame.height, tissues=stats.height,
                        least=stats["n"].min(), most=stats["n"].max())  # fmt: skip
    strip = StripData(frame=frame, columns=_columns("tissues"), category="tissue", value="tpm",
                      category_label=text.word(PAGE, "figures.tissues.x"),
                      value_label=text.word(PAGE, "figures.tissues.y"), colour="sex",
                      colours=category_colours(["female", "male"]), caption=caption,
                      hover=("donor",), ranges=ranges,
                      range_label=text.word(PAGE, "figures.tissues.range"), log=True,
                      width="double")  # fmt: skip
    words = _words(
        "tissues",
        top=by_median["tissue"][0],
        top_median=by_median["median"][0],
        bottom=by_median["tissue"][-1],
        bottom_median=by_median["median"][-1],
    )
    return writer.figure(strip, "normal_tissues", words)


def immune(writer: Writer, rng: np.random.Generator) -> tuple[FigureBlock, float]:
    """Made-up HPA immune cell nTPM, one pooled bar per cell type, with our tau; no test."""
    names, lineages, base = zip(*IMMUNE, strict=True)
    ntpm = np.round(np.array(base) * np.exp(rng.normal(0.0, IMMUNE_SD, len(base))), 1)
    frame = pl.DataFrame({"cell_type": names, "lineage": lineages, "ntpm": ntpm})
    specificity, top = float(tau(ntpm) or 0.0), int(np.argmax(ntpm))
    caption = text.fill(PAGE, "figures.immune.caption", cells=len(names), tau=specificity)
    bars = BarData(frame=frame, columns=_columns("immune"), category="cell_type", value="ntpm",
                   category_label=text.word(PAGE, "figures.immune.x"),
                   value_label=text.word(PAGE, "figures.immune.y"), series="lineage",
                   colours=dict(zip(dict.fromkeys(lineages), LINEAGE_COLOURS, strict=True)),
                   caption=caption, horizontal=True, stacked=True, width="double")  # fmt: skip
    words = _words("immune", top=names[top], top_ntpm=ntpm[top], tau=specificity, cells=len(names))
    return writer.figure(bars, "immune_cells", words), specificity


def exons(writer: Writer, rng: np.random.Generator) -> tuple[FigureBlock, str]:
    """Made-up exon shares per tissue against each exon's median, with CV; no test.

    Returns the figure and the adjacent exon pair whose larger CV is smallest, for qPCR.
    """
    names, tissue_names = [f"E{i + 1}" for i in range(len(EXON_NT))], list(TISSUES)
    reads = np.outer(EXON_NT, np.ones(len(TISSUES)))
    reads[ALT_EXON] *= [use for _, use in TISSUES.values()]
    reads[:, tissue_names.index(DEGRADED)] *= 1.0 + BIAS_STEP * np.arange(len(EXON_NT))
    reads *= np.exp(rng.normal(0.0, EXON_SD, reads.shape))
    share = reads / reads.sum(axis=0)
    ratio = np.log2(share / np.median(share, axis=1, keepdims=True))
    cv = share.std(axis=1, ddof=1) / share.mean(axis=1)
    levels, swatches = diverging_levels(ratio.ravel().tolist())
    frame = pl.DataFrame({"exon": np.repeat(names, len(TISSUES)),
                          "tissue": tissue_names * len(names), "share": share.ravel(),
                          "log2_ratio": ratio.ravel(), "level": levels})  # fmt: skip
    first = int(np.argmin(np.maximum(cv[:-1], cv[1:])))
    caption = text.fill(PAGE, "figures.exons.caption", exons=len(names), tissues=len(TISSUES))
    heat = HeatmapData(frame=frame, columns=_columns("exons"), row="exon", column="tissue",
                       value="log2_ratio", level="level",
                       row_label=text.word(PAGE, "figures.exons.y"),
                       column_label=text.word(PAGE, "figures.exons.x"), swatches=swatches,
                       caption=caption, row_order=tuple(names), column_order=tuple(tissue_names),
                       row_stats=dict(zip(names, cv.tolist(), strict=True)),
                       row_stats_label=text.word(PAGE, "figures.exons.stat"),
                       row_stats_format=".2f")  # fmt: skip
    pair = {"first": names[first], "second": names[first + 1]}
    words = _words("exons", **pair, uneven=names[int(np.argmax(cv))], uneven_cv=float(cv.max()))
    return writer.figure(heat, "exon_use", words), text.fill(PAGE, "keys.pair_value", **pair)


def page(root: Path, dataset: str, day: date) -> Path:
    """Write the M9 preview page and its figures under root; return the page."""
    started = perf_counter()
    gathered = root / "outputs" / "01_gather"
    sheet_file = sorted(gathered.glob(f"*_gather_{dataset}_samplesheet_v1.parquet"))[-1]
    sheet = pl.read_parquet(sheet_file).unique("sample", keep="first", maintain_order=True)
    out = root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    writer = Writer(out, dataset, day, {"made_up_seed": str(SEED),
                                        "input:samplesheet": sheet_file.name})  # fmt: skip
    rng = np.random.default_rng(SEED)
    own, normal = ours(writer, rng, sheet), tissues(writer, rng)
    (cells, specificity), (heat, pair) = immune(writer, rng), exons(writer, rng)
    group = {k: text.word(PAGE, f"groups.{k}") for k in ("ours", "reference", "exons")}
    lead = {k: text.word(PAGE, f"groups.{k}_lead") for k in group}
    sources = ResultTable(text.word(PAGE, "sources.title"), text.words(PAGE, "sources.columns"),
                          tuple(tuple(r) for r in text.raw(PAGE, "sources.rows")))  # fmt: skip
    results = (text.word("ui", "preview.note"),
               Group(group["ours"], (own,), lead=lead["ours"]),
               Group(group["reference"], (normal, cells), lead=lead["reference"]),
               Group(group["exons"], (heat,), lead=lead["exons"]), sources)  # fmt: skip
    keys = (KeyNumber(text.word(PAGE, "keys.samples"), f"{sheet.height:,}"),
            KeyNumber(text.word(PAGE, "keys.tissues"), f"{len(TISSUES):,}"),
            KeyNumber(text.word(PAGE, "keys.tau"), f"{specificity:.2f}"),
            KeyNumber(text.word(PAGE, "keys.pair"), pair))  # fmt: skip
    cards = tuple(ParameterCard(flag, c["meaning"], c["default"], c["effect"])
                  for flag, c in text.raw(PAGE, "planned").items())  # fmt: skip
    usage, budget = runtime.measure(started), runtime.budget(MODULE)
    built = ModulePage(
        title=text.fill(PAGE, "title", dataset=text.word(*SHOWN)), **text.skeleton(PAGE),
        parameters=cards,
        results=results, key_numbers=keys,
        reproduce=Reproduction(
            command=text.fill(PAGE, "command", gene=GENE, dataset=dataset),
            wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children, profile=text.word(PAGE, "profile"),
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(), inputs=(sheet_file.name,),
        ),
    )  # fmt: skip
    written = write(built, out / f"{stem(MODULE, f'{dataset}_overview', 1, day=day)}.html")
    # A kit refuses this folder without its entry page (D-072), as for a real root.
    write_entry(out, MODULE)
    return written
