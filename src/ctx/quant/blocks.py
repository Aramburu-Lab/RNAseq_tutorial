"""What the M4 page shows of each figure: its frame, test, caption and lead, then its block.

F1 and F3 read files every nf-core/rnaseq outdir carries; F2, F4, F6 and F7 read copy Q1
through ctx's port of tximport, which ``ctx.quant.checks`` has held equal to nf-core's tables.
Every word is ``content/quant.toml``'s.
"""

from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import polars as pl

from ctx.concordance import rank_correlation
from ctx.de.deseq import select_samples
from ctx.de.filters import expressed
from ctx.de.run import Writer
from ctx.quant import CountsFrom, outdir
from ctx.quant import figures as fig
from ctx.quant import tables as tab
from ctx.quant.tximport import GeneTables
from ctx.report import text
from ctx.report.builder import FigureBlock, FigureText, ResultTable
from ctx.viz.aggregates import Aggregate, Lead
from ctx.viz.annotate import Caption, compare_groups

PAGE = "quant"
# Samples each read-length batch needs within a cell type before the two are compared.
MIN_PER_BATCH = 3
# A fold change differing by this much (log2) is always drawn on the page: the genes F6's lead
# counts, two-fold or more.
FAR = 1.0


@dataclass(frozen=True)
class Shown:
    """The sheet as the figures name it.

    Attributes:
        sheet: sample, cell_type, condition, batch, donor, one row per sample in sheet order.
        groups: sample, then each sheet column the bars filter by.
        batches: (reference, tested) read-length batches, in sheet order; empty if not two.
    """

    sheet: pl.DataFrame
    groups: pl.DataFrame
    batches: tuple[str, ...]

    @property
    def cells(self) -> list[str]:
        """Cell types in sheet order."""
        return self.sheet["cell_type"].unique(maintain_order=True).to_list()


def figure_text(key: str, lead: str = "", **values: Any) -> FigureText:
    """Title, lead and reading guide of ``figures.<key>``."""
    reading = tuple(line.format_map(values) for line in text.words(PAGE, f"figures.{key}.reading"))
    return FigureText(title=text.word(PAGE, f"figures.{key}.title"), lead=lead, reading=reading)


def batch_tests(frame: pl.DataFrame, value: str, shown: Shown) -> Caption | None:
    """Mann-Whitney by read-length batch within each cell type with enough of both, BH."""
    if len(shown.batches) != 2:  # noqa: PLR2004 - the test compares two batches
        return None
    sized = frame.group_by("cell_type", "batch").len()
    ready = sized.filter(pl.col("len") >= MIN_PER_BATCH).group_by("cell_type").len()
    cells = [c for c in shown.cells if c in set(ready.filter(pl.col("len") == 2)["cell_type"])]  # noqa: PLR2004
    if not cells:
        return None
    kept = frame.filter(pl.col("cell_type").is_in(cells))
    return compare_groups(kept, feature="cell_type", group="batch", value=value,
                          groups=(shown.batches[0], shown.batches[1]))  # fmt: skip


def test_table(found: Caption, key: str) -> ResultTable:
    """One row per cell type tested: n per batch, U, p, BH p and rank-biserial r."""
    rows = tuple((t.feature, str(t.n_reference), str(t.n_comparison), f"{t.u:g}", f"{t.p:.2g}",
                  f"{t.p_adj:.2g}", f"{t.r:.2f}") for t in found.tests)  # fmt: skip
    return ResultTable(text.word(PAGE, f"tables.{key}.title"),
                       tuple(text.words(PAGE, "tables.tests.columns")), rows,
                       numeric=(1, 2, 3, 4, 5, 6))  # fmt: skip


def tests_frame(found: Caption) -> pl.DataFrame:
    """A batch test as a table: one row per cell type."""
    return pl.DataFrame([{"cell_type": t.feature, "n_reference": t.n_reference,
                          "n_tested": t.n_comparison, "u": t.u, "p": t.p, "p_bh": t.p_adj,
                          "r": t.r} for t in found.tests])  # fmt: skip


def caption(key: str, found: Caption | None, **values: Any) -> str:
    """A figure's caption, with its batch tests, or the variant that says none could run."""
    if found is None:
        return text.fill(PAGE, f"figures.{key}.caption_untested", **values)
    return text.fill(PAGE, f"figures.{key}.caption", **values, **fig.test_caption(found))


def fates_figure(writer: Writer, fates: pl.DataFrame, shown: Shown
                 ) -> tuple[FigureBlock, Caption | None]:  # fmt: skip
    """F1, its lead following the filter: the median mapped share and the lowest."""
    frame = tab.fates_long(fates, shown.groups)
    rate = fates.select("sample", mapped_pct=100 * pl.col("mapped") / pl.col("processed")).join(
        shown.sheet, on="sample")  # fmt: skip
    found = batch_tests(rate, "mapped_pct", shown)
    data = fig.fate_bars(frame, shown.groups.columns[1:],
                         caption("fates", found, n=fates.height))  # fmt: skip
    mapped = {"fate": "mapped"}
    lead = Lead(text.raw(PAGE, "figures.fates.lead"),
                (Aggregate("median", "median", "percent", mapped),
                 Aggregate("low", "min", "percent", mapped, row=("low_sample", "sample"))),
                note=text.word(PAGE, "figures.fates.note"))  # fmt: skip
    return writer.figure(data, "fates", replace(figure_text("fates"), follows=lead)), found


def correlation_figure(writer: Writer, matrix: np.ndarray, shown: Shown) -> FigureBlock:
    """F3: Spearman between samples; cell-type labels permuted within each donor."""
    rho, genes = tab.spearman(matrix)
    sheet = shown.sheet
    strata = [d or s for d, s in zip(sheet["donor"].fill_null(""), sheet["sample"], strict=True)]
    found = tab.within_minus_between(rho, sheet["cell_type"].to_list(), strata)
    if found is None:
        words = text.fill(PAGE, "figures.correlation.caption_untested", n=sheet.height,
                          genes=genes, min_reads=tab.MIN_READS)  # fmt: skip
        lead = text.word(PAGE, "figures.correlation.lead_untested")
    else:
        words = text.fill(PAGE, "figures.correlation.caption", n=sheet.height, genes=genes,
                          diff=found.difference, p=found.p, perms=found.permutations,
                          strata=found.strata, min_reads=tab.MIN_READS)  # fmt: skip
        lead = text.fill(PAGE, "figures.correlation.lead", within=found.within,
                         between=found.between)  # fmt: skip
    data = fig.correlation_heatmap(rho, sheet, words)
    return writer.figure(data, "correlation", figure_text("correlation", lead))


def detection_figure(writer: Writer, port: GeneTables, shown: Shown
                     ) -> tuple[FigureBlock, pl.DataFrame]:  # fmt: skip
    """F2: genes detected as each sample's estimated counts are thinned; tested at the common."""
    frame, common = tab.thinned(port.counts, port.samples)
    frame = frame.join(shown.sheet.select("sample", "cell_type", "condition", "batch"),
                       on="sample", maintain_order="left")  # fmt: skip
    at = frame.filter("common")
    found = batch_tests(at, "genes", shown)
    words = caption("detection", found, n=at.height, depth=f"{common:,.0f}",
                    min_reads=tab.MIN_READS)  # fmt: skip
    data = fig.detection_curves(frame.drop("common"), shown.cells, common, words)
    genes = at["genes"]
    lead = text.fill(PAGE, "figures.detection.lead", depth=f"{common:,.0f}",
                     median=f"{tab.stat(genes, 'median'):,.0f}",
                     low=f"{tab.stat(genes, 'min'):,.0f}",
                     high=f"{tab.stat(genes, 'max'):,.0f}")  # fmt: skip
    block = writer.figure(data, "detection", figure_text("detection", lead,
                                                         min_reads=tab.MIN_READS))  # fmt: skip
    return block, frame


def units_figure(writer: Writer, port: GeneTables, fates: pl.DataFrame,
                 names: dict[str, str]) -> tuple[FigureBlock, pl.DataFrame]:  # fmt: skip
    """F4: the median-depth sample's genes, reads against TPM; descriptive."""
    sample = fates.sort("mapped", "sample")["sample"][fates.height // 2]
    j = port.samples.index(sample)
    whole = pl.DataFrame({"gene_id": port.genes, "tpm": port.tpm[:, j],
                          "counts": port.counts[:, j], "length": port.lengths[:, j]})  # fmt: skip
    whole = whole.filter(pl.col("counts") >= tab.MIN_READS).with_columns(
        gene_name=pl.col("gene_id").replace(names))  # fmt: skip
    whole, classes = fig.length_classes(whole)
    per_tpm = dict(whole.group_by("length_class").agg((pl.col("counts") / pl.col("tpm"))
                                                       .median()).iter_rows())  # fmt: skip
    drawn = whole.select("gene_id", "gene_name", "tpm", "counts", "length", "length_class",
                         keep=pl.lit(value=False))  # fmt: skip
    data = fig.units_scatter(drawn, classes, text.fill(PAGE, "figures.units.caption",
                                                       sample=sample, total=whole.height,
                                                       min_reads=tab.MIN_READS))  # fmt: skip
    ratio = per_tpm.get(classes["long"], float("nan")) / per_tpm.get(classes["short"], float("nan"))
    lead = text.fill(PAGE, "figures.units.lead", sample=sample, ratio=ratio)
    return writer.figure(data, "units", figure_text("units", lead)), whole


def contrast_samples(inputs: Any) -> pl.DataFrame:
    """The featured contrast's samples, as M5 picks them, with ``group``."""
    used, _ = select_samples(inputs.samples, inputs.contrast, inputs.spec.order)
    return used.select("sample", "group")


def estimates_figure(writer: Writer, port: GeneTables, used: pl.DataFrame, genes: pl.DataFrame,
                     name: str) -> tuple[FigureBlock, pl.DataFrame]:  # fmt: skip
    """F6: the contrast's fold change from estimated and from length-scaled counts, per gene."""
    columns = [port.samples.index(s) for s in used["sample"]]
    tested = (used["group"] == "tested").to_numpy()
    scaled = port.counts_from(CountsFrom.LENGTH_SCALED_TPM)[:, columns]
    estimated = port.counts_from(CountsFrom.NO)[:, columns]
    keep = expressed(scaled, [int((~tested).sum()), int(tested.sum())])
    frame = pl.DataFrame({
        "gene_id": np.asarray(port.genes)[keep],
        "estimated": tab.fold_change(estimated, tested)[keep],
        "length_scaled": tab.fold_change(scaled, tested)[keep],
    }).with_columns(difference=pl.col("estimated") - pl.col("length_scaled")).join(
        genes.select("gene_id", "gene_name", "gene_type"), on="gene_id", how="left")  # fmt: skip
    frame = fig.type_names(frame)
    rho = float(rank_correlation(frame["estimated"].to_numpy(), frame["length_scaled"].to_numpy()))
    frame = frame.with_columns(keep=pl.col("difference").abs() >= FAR)
    far = int(frame["keep"].sum())
    data = fig.estimates_scatter(frame, text.fill(
        PAGE, "figures.estimates.caption", contrast=name, n_ref=int((~tested).sum()),
        n_test=int(tested.sum()), genes=frame.height, rho=rho))  # fmt: skip
    lead = text.fill(PAGE, "figures.estimates.lead", far=far, genes=frame.height, rho=rho)
    return writer.figure(data, "estimates", figure_text("estimates", lead)), frame


def isoforms_figure(writer: Writer, folder: Any, samples: Sequence[str], genes: pl.DataFrame
                    ) -> tuple[FigureBlock, pl.DataFrame] | None:  # fmt: skip
    """F7: isoforms expressed against isoforms annotated, coding against lncRNA."""
    half = len(samples) / 2
    expressed_tx = (pl.scan_csv(folder / outdir.MERGED.format("transcript_tpm"), separator="\t",
                                schema_overrides=dict.fromkeys(samples, pl.Float64))
                    .select("tx", expressed=pl.sum_horizontal([pl.col(s) >= tab.MIN_TPM
                                                               for s in samples]) >= half)
                    .collect())  # fmt: skip
    per_gene = fig.type_names(tab.isoforms(expressed_tx, outdir.tx2gene(folder), genes))
    names = [text.word(PAGE, f"gene_types.{t}") for t in tab.GENE_TYPES]
    typed = per_gene.filter(pl.col("gene_type").is_in(names))
    # Genes with no isoform expressed share 0 %: on Step-1 that is 12,904 of 18,253 coding
    # genes and the median of both types, so the share is read where a gene is expressed.
    compared = typed.filter(pl.col("expressed") > 0).with_columns(feature=pl.lit(
        text.word(PAGE, "figures.isoforms.feature")), share_pct=100 * pl.col("share"))  # fmt: skip
    if compared["gene_type"].n_unique() < 2:  # noqa: PLR2004 - two types compared
        return None
    found = compare_groups(compared, feature="feature", group="gene_type", value="share_pct",
                           groups=(names[0], names[1]))  # fmt: skip
    bins = compared.with_columns(annotated_bin=pl.col("annotated").map_elements(
        lambda a: text.words(PAGE, "figures.isoforms.bins")[tab.isoform_bin(a)],
        return_dtype=pl.String))  # fmt: skip
    order = {b: i for i, b in enumerate(text.words(PAGE, "figures.isoforms.bins"))}
    bars = (bins.group_by("annotated_bin", "gene_type")
            .agg(share=pl.col("share_pct").median(), n_genes=pl.len())
            .sort(pl.col("annotated_bin").replace_strict(order), "gene_type"))  # fmt: skip
    medians = dict(compared.group_by("gene_type").agg(pl.col("share_pct").median()).iter_rows())
    sizes = dict(typed.group_by("gene_type").len().iter_rows())
    shown = dict(compared.group_by("gene_type").len().iter_rows())
    data = fig.isoform_bars(bars, caption("isoforms", found, n=compared.height,
                                          total=typed.height, tpm=tab.MIN_TPM))  # fmt: skip
    lead = text.fill(PAGE, "figures.isoforms.lead", coding=medians[names[0]],
                     lncrna=medians[names[1]], coding_n=shown[names[0]],
                     coding_total=sizes[names[0]], lncrna_n=shown[names[1]],
                     lncrna_total=sizes[names[1]])  # fmt: skip
    return writer.figure(data, "isoforms", figure_text("isoforms", lead, tpm=tab.MIN_TPM)), per_gene
