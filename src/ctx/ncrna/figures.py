"""The data contracts of M8's figures, from the run's inputs; the words come from content.

F1 each sample's counted reads by biotype; F2 the value M5's filter tests, by biotype; F3 M5's
featured contrast by biotype, BH within each; F4 a lncRNA against its nearest protein-coding
gene, the pair chosen on the page. Each builder returns the figure's data and its numbers.
"""

import base64
import math
from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import mannwhitneyu

from ctx.de.filters import TOLERANCE, expressed, min_sample_size
from ctx.de.stats import ALPHA, LFC_MIN, P_FLOOR, adjust_within, classify
from ctx.function.live import bitset
from ctx.ncrna import biotypes, neighbours
from ctx.report import text
from ctx.viz.bar import BarData
from ctx.viz.curve import CurveData
from ctx.viz.scatter import CUTOFF_PRESETS, MIN_CUTOFF, Line, ScatterData, needed_points
from ctx.viz.style import ALONE_SCALE, LARGEST_PLOT_PX, PLOTLY_HEIGHT_PX

PAGE = "ncrna"
NOT_CALLED = "ns"
# The sheet columns F1 carries where the sheet has them, Step-1's all three.
FACETS = ("cell_type", "condition", "library_selection")
# One meaning per column, for every table's and figure's README.
COLUMNS = {
    "sample": ("none", "Sample, from the gather sample sheet."),
    "cell_type": ("none", "Cell type, from the sample sheet."),
    "condition": ("none", "Condition, from the sample sheet."),
    "library_selection": ("none", "ENA library_selection, from the sample sheet."),
    "biotype": ("none", "GENCODE gene type, grouped (ctx.ncrna.biotypes)."),
    "percent": ("%", "Share of the sample's counted reads."),
    "mirna_percent": ("%", "Share of the sample's counted reads on miRNA genes, hairpin loci."),
    **{
        f"{g}_percent": ("%", f"Share of the sample's counted reads on {g} genes.")
        for g in biotypes.GROUPS
    },
    "gene_id": ("none", "GENCODE gene ID, versioned."),
    "gene_name": ("none", "GENCODE gene name."),
    "gene_type": ("none", "GENCODE gene type."),
    "value": (
        "log2 reads",
        "log2 of 1 plus the reads, at the median library size, that the "
        "filter's sample count reaches.",
    ),
    "kept": ("none", "True where the filter, at this run's minimum count, keeps the gene."),
    "share": ("fraction", "Share of the biotype's genes at or below this value."),
    "baseMean": ("normalised counts", "M5: mean of the normalised counts."),
    "log2FoldChange": ("log2 fold change", "M5: raw log2 fold change, tested over reference."),
    "pvalue": ("dimensionless", "M5: Wald test p-value."),
    "padj": ("dimensionless", "M5: Benjamini-Hochberg over every tested gene, pooled."),
    "padj_within": ("dimensionless", "Benjamini-Hochberg within the gene's biotype."),
    "drawn": ("none", "True where the page draws the gene at the build's setting."),
    "call_pooled": ("none", "up, down or ns by the call rule on padj."),
    "call": ("none", "up, down or ns by the call rule on padj_within."),
    "neg_log10_p": (
        "-log10 p",
        f"Minus log10 of the raw p-value, drawn at {P_FLOOR:g} where it is smaller.",
    ),
    "level": ("none", "The biotype of a gene called on padj_within; ns otherwise."),
    "keep": ("none", "True where the page draws the gene whatever the thinning: a call."),
    "lncrna": ("log2 CPM", "The lncRNA, log2 of 1 plus counts per million."),
    "neighbour": ("log2 CPM", "Its nearest protein-coding gene, log2 of 1 plus CPM."),
    "lncrna_id": ("none", "The lncRNA's GENCODE gene ID."),
    "lncrna_name": ("none", "The lncRNA's GENCODE gene name."),
    "neighbour_id": ("none", "The protein-coding gene's GENCODE gene ID."),
    "neighbour_name": ("none", "The protein-coding gene's GENCODE gene name."),
    "gap": ("nt", "Nucleotides between the two genes; 0 when they touch or overlap."),
    "kind": ("none", "antisense or same_strand overlap, near, far, or trans (background)."),
    "rho": ("dimensionless", "Spearman rho over the contrast's samples, log2 CPM."),
    "p": ("dimensionless", "Two-sided p of rho, t approximation."),
    "q": ("dimensionless", "BH of p over the pairs of its set: cis or background."),
    "rho_within": ("dimensionless", "Spearman rho of the residuals from each condition's mean."),
    "p_within": ("dimensionless", "Two-sided p of rho_within, t approximation."),
    "q_within": ("dimensionless", "BH of p_within over the pairs of its set."),
}


def columns(frame: pl.DataFrame) -> dict[str, tuple[str, str]]:
    """(unit, meaning) of every column of a table."""
    return {c: COLUMNS[c] for c in frame.columns}


def word(key: str) -> str:
    """A fixed text of this page's content file."""
    return text.word(PAGE, key)


def composition(sheet: pl.DataFrame, genes: pl.DataFrame, matrix: np.ndarray,
                caption: str) -> tuple[BarData, pl.DataFrame]:  # fmt: skip
    """F1: each sample's counted reads by biotype, filterable by biotype, cell type and condition.

    Returns:
        The figure's data, and one row per sample with every group's share and the miRNA
        genes' share, which rule N1 reads.
    """
    groups = biotypes.shares(matrix, genes["biotype"].to_numpy(), biotypes.GROUPS)
    (mirna,) = biotypes.shares(matrix, genes["gene_type"].to_numpy(), (biotypes.MIRNA,))
    # A GDC sheet has neither cell type nor library selection: its bars split by condition.
    facets = [c for c in FACETS if c in sheet.columns]
    shown = tuple(c for c in facets if c != "library_selection")
    wide = sheet.select("sample", *facets).with_columns(
        *(pl.Series(f"{g}_percent", groups[i]) for i, g in enumerate(biotypes.GROUPS)),
        mirna_percent=pl.Series(mirna),
    )
    frame = (wide.unpivot([f"{g}_percent" for g in biotypes.GROUPS],
                          index=["sample", *shown], variable_name="biotype",
                          value_name="percent")
             .with_columns(pl.col("biotype").str.strip_suffix("_percent"))
             .sort(*shown, "sample", maintain_order=True))  # fmt: skip
    data = BarData(
        frame=frame, columns=columns(frame), category="sample", value="percent",
        category_label=word("figures.composition.x"), value_label=word("figures.composition.y"),
        series="biotype", colours=biotypes.COLOURS, caption=caption, stacked=True,
        hover=shown, filters=(*shown, "biotype"),
        summary_median=False,
    )  # fmt: skip
    return data, wide


def filter_values(genes: pl.DataFrame, sub: np.ndarray, sizes: list[int],
                  min_count: int) -> tuple[pl.DataFrame, int]:  # fmt: skip
    """The value M5's filter tests, per gene, and whether it passes at ``min_count``.

    The value is the count, at the median library size, that the filter's k samples reach:
    the k-th largest of the gene's scaled counts. On that scale the cutoff is ``min_count``
    itself, whatever the depth.

    Args:
        genes: The annotation's record of each count row, with ``biotype``.
        sub: Counts of the contrast's samples, genes by samples.
        sizes: Samples per compared group.
        min_count: The cutoff, in reads at the median library size.

    Returns:
        gene_id, gene_name, biotype, value (log2 of 1 plus reads) and kept; and k.
    """
    k = math.ceil(min_sample_size(sizes) - TOLERANCE)
    library = sub.sum(axis=0)
    at_median = np.sort(sub / library * float(np.median(library)), axis=1)[:, -k]
    frame = genes.select("gene_id", "gene_name", "biotype").with_columns(
        value=pl.Series(np.log2(at_median + 1)), kept=pl.Series(expressed(sub, sizes, min_count))
    )
    return frame, k


def curves(values: pl.DataFrame, chosen: Sequence[str], cutoff: float, caption: str) -> CurveData:
    """F2: each chosen biotype's cumulative share of genes against the filter's value.

    Each line keeps the points that draw it to within a pixel at the largest plot
    (``needed_points``); the per-gene values are their own table.
    """
    parts = []
    for b in chosen:
        x = np.sort(values.filter(pl.col("biotype") == b)["value"].to_numpy())
        if not x.size:
            continue
        y = np.arange(1, x.size + 1) / x.size
        kept = needed_points(x, y, (LARGEST_PLOT_PX, PLOTLY_HEIGHT_PX["scatter"] * ALONE_SCALE))
        parts.append(pl.DataFrame({"biotype": b, "value": x[kept], "share": y[kept]}))
    frame = pl.concat(parts)
    return CurveData(
        frame=frame, columns=columns(frame), series="biotype", x="value", y="share",
        x_label=word("figures.filter.x"), y_label=word("figures.filter.y"), colour="biotype",
        colours={b: biotypes.COLOURS[b] for b in chosen}, caption=caption,
        thresholds=(Line(word("figures.filter.cutoff"), (cutoff, cutoff), (0.0, 1.0),
                         dash=True),),
    )  # fmt: skip


def rank_test(values: pl.DataFrame) -> dict[str, float]:
    """Mann-Whitney U of the filter's value, lncRNAs against protein-coding genes; r too.

    r is the rank-biserial correlation, 2U / (n1 n2) - 1, positive when lncRNAs rank higher.
    """
    lnc, coding = (values.filter(pl.col("biotype") == b)["value"].to_numpy()
                   for b in (biotypes.LNCRNA, biotypes.CODING))  # fmt: skip
    u, p = mannwhitneyu(lnc, coding)
    return {"p": float(p), "r": float(2 * u / (lnc.size * coding.size) - 1),
            "n_lncrna": lnc.size, "n_coding": coding.size}  # fmt: skip


def by_biotype(results: pl.DataFrame, genes: pl.DataFrame, chosen: Sequence[str]) -> pl.DataFrame:
    """M5's tested genes of the chosen biotypes, BH within each, called both ways.

    The tested genes are those with M5's adjusted p; M5's own ``padj`` is BH over all of them,
    pooled, and stays beside the stratified one.
    """
    tested = (results.filter(pl.col("padj").is_not_null())
              .join(genes.select("gene_id", "biotype"), on="gene_id", how="left")
              .filter(pl.col("biotype").is_in(list(chosen))))  # fmt: skip
    within = adjust_within(tested["pvalue"].to_numpy(), tested["biotype"].to_numpy(), "padj")
    return (
        tested.select("gene_id", "gene_name", "biotype", "baseMean", "log2FoldChange", "pvalue",
                      "padj")
        .with_columns(padj_within=pl.Series(within))
        .with_columns(call_pooled=classify(pl.col("padj"), pl.col("log2FoldChange")),
                      call=classify(pl.col("padj_within"), pl.col("log2FoldChange")),
                      neg_log10_p=-pl.max_horizontal(pl.col("pvalue"), pl.lit(P_FLOOR)).log10())
        .with_columns(level=pl.when(pl.col("call") == NOT_CALLED).then(pl.lit(NOT_CALLED))
                      .otherwise(pl.col("biotype")), keep=pl.col("call") != NOT_CALLED)
    )  # fmt: skip


def volcano(frame: pl.DataFrame, chosen: Sequence[str], words: Mapping[str, Any],
            caption: str) -> ScatterData:  # fmt: skip
    """F3: the featured contrast's volcano, a called gene coloured by its biotype.

    M5's cut-offs move live (Rytis, TODO.md M8): ctxCalls calls every tested gene of the
    chosen biotypes on a p adjusted within its biotype, and the page colours the called ones by
    biotype (``_page_strata.js.j2``). Every tested gene travels in the control's rows, since a
    looser cut-off calls genes the build's thinning left out.
    """
    top = float(frame["neg_log10_p"].to_numpy().max(initial=1.0)) * 1.05
    name = word("figures.volcano.threshold")
    data = ScatterData(
        frame=frame, columns=columns(frame), x="log2FoldChange", y="neg_log10_p",
        x_label=text.fill(PAGE, "figures.volcano.x", **words),
        y_label=word("figures.volcano.y"), colour="level",
        colours={**{b: biotypes.COLOURS[b] for b in chosen}, NOT_CALLED: None},
        names={NOT_CALLED: word("figures.volcano.ns")}, caption=caption,
        hover=("gene_name", "padj_within", "padj"), keep="keep",
        lines=tuple(Line(name, (x, x), (0.0, top), dash=True) for x in (-LFC_MIN, LFC_MIN)),
    )  # fmt: skip
    return replace(data, live=strata(data, chosen))


# The bases the page offers: BH within the biotype (Rytis, TODO.md M8), then the raw p. The
# other adjustments within each biotype would add 3 columns of float64 per gene, and put
# TCGA-BRCA's page at 14,660,097 B with them, 340 kB under the cap.
STRATA_BASES = ("padj_within", "pvalue")


def strata_summary(frame: pl.DataFrame, chosen: Sequence[str], basis: str, alpha: float,
                   lfc_min: float) -> str:  # fmt: skip
    """The live line at one setting: each biotype's calls, as the page fills it again."""
    calls = frame.select("biotype", call=classify(pl.col(basis), pl.col("log2FoldChange"),
                                                  alpha, lfc_min))  # fmt: skip
    entries = []
    for b in chosen:
        mine = calls.filter(pl.col("biotype") == b)["call"]
        up, down = int((mine == "up").sum()), int((mine == "down").sum())
        entries.append(text.fill(PAGE, "figures.volcano.entry", biotype=b, n=up + down, up=up,
                                 down=down))  # fmt: skip
    return text.fill(
        PAGE, "figures.volcano.summary", basis=word(f"figures.volcano.bases.{basis}"),
        cutoff=f"{alpha:g}", fc=f"{lfc_min:.2f}", counts="; ".join(entries),
    )  # fmt: skip


def packed(rows: pl.DataFrame, chosen: Sequence[str]) -> dict[str, Any]:
    """F3's rows as the page receives them, packed.

    Names joined, biotypes by index, the numbers the call rule reads as float64 and the hover's
    pooled padj as float32, in base64. As JSON, Step-1's 12,246 rows took 2,177,184 B, and
    TCGA-BRCA's 19,650 would have put its page over the cap.
    """

    def encoded(column: str, kind: str) -> str:
        return base64.b64encode(rows[column].to_numpy().astype(kind).tobytes()).decode("ascii")

    return {
        "n": rows.height, "names": "\n".join(rows["gene_name"].fill_null("").to_list()),
        "biotype": [chosen.index(b) for b in rows["biotype"]],
        "f64": {c: encoded(c, "<f8") for c in ("log2FoldChange", *STRATA_BASES)},
        "f32": {"padj": encoded("padj", "<f4")},
        "drawn": bitset([i for i, d in enumerate(rows["drawn"]) if d], rows.height),
    }  # fmt: skip


def strata(data: ScatterData, chosen: Sequence[str]) -> dict[str, Any]:
    """F3's live cut-offs: the rows a cut-off can call, the bases, the words the page fills.

    A gene whose raw p is at or above the loosest cut-off is called at no setting, as every
    adjusted p is at least the raw one: it travels only where the build draws it.
    """
    frame = data.frame.with_columns(drawn=data.thinning.page)
    callable_ = pl.col("drawn") | (pl.col("pvalue") < max(CUTOFF_PRESETS))
    carried = frame.filter(callable_)
    traces = {level: i for i, level in enumerate(data.colours)}
    # Moved, the files draw the genes no cut-off calls and the page leaves out, as the ns
    # trace's points (export.settings' thinned layer: float32 x, y pairs by trace index).
    left = frame.filter(~callable_).select("log2FoldChange", "neg_log10_p").to_numpy()
    never = {str(traces[NOT_CALLED]):
             base64.b64encode(left.astype("<f4").tobytes()).decode("ascii")}  # fmt: skip
    return {
        "kind": "strata", "lfc": "log2FoldChange", "y": "neg_log10_p", "gate": "padj_within",
        "stratum": "biotype", "levels": list(chosen), "ns": NOT_CALLED, "traces": traces,
        "lines": [len(data.colours), len(data.colours) + 1], "hover": list(data.hover),
        "bases": {b: word(f"figures.volcano.bases.{b}") for b in STRATA_BASES},
        "alpha": ALPHA, "lfc_min": LFC_MIN, "max_cutoff": max(CUTOFF_PRESETS),
        "min_cutoff": MIN_CUTOFF, "presets": [f"{c:g}" for c in CUTOFF_PRESETS],
        "p_floor": P_FLOOR, "rows": packed(carried, chosen), "never": never,
        "summary": text.raw(PAGE, "figures.volcano.summary"),
        "entry": text.raw(PAGE, "figures.volcano.entry"),
        "every": word("figures.volcano.every"),
        "default_summary": strata_summary(frame, chosen, "padj_within", ALPHA, LFC_MIN),
    }  # fmt: skip


def calls_table(frame: pl.DataFrame, chosen: Sequence[str]) -> list[tuple[str, ...]]:
    """Per biotype: tested, called on pooled BH, called within, up and down within."""
    rows: list[tuple[str, ...]] = []
    for b in chosen:
        rows_b = frame.filter(pl.col("biotype") == b)
        rows.append((b, f"{rows_b.height:,}", f"{(rows_b['call_pooled'] != NOT_CALLED).sum():,}",
                     f"{(rows_b['call'] != NOT_CALLED).sum():,}",
                     f"{(rows_b['call'] == 'up').sum():,}",
                     f"{(rows_b['call'] == 'down').sum():,}"))  # fmt: skip
    return rows


def pair_table(genes: pl.DataFrame, results: pl.DataFrame, logcpm: np.ndarray,
               strata: np.ndarray) -> pl.DataFrame:  # fmt: skip
    """Every tested lncRNA with its nearest tested protein-coding gene, and the background.

    Args:
        genes: The annotation's record of each count row, with ``biotype``.
        results: M5's results of the contrast; its tested genes are the candidates.
        logcpm: log2 of 1 plus CPM, count rows by the contrast's samples.
        strata: Each sample's stratum, whose mean the residuals remove.

    Returns:
        One row per pair, cis pairs first, then the background (kind ``trans``), each set
        with its own BH.
    """
    tested = results.filter(pl.col("padj").is_not_null())["gene_id"]
    loci = genes.filter(pl.col("gene_id").is_in(tested.implode()))
    lnc = loci.filter(pl.col("biotype") == biotypes.LNCRNA)
    coding = loci.filter(pl.col("biotype") == biotypes.CODING)
    cis = neighbours.nearest(lnc, coding)
    row = {g: i for i, g in enumerate(genes["gene_id"].to_list())}
    name = dict(genes.select("gene_id", "gene_name").rows())
    parts = []
    for pairs in (cis, neighbours.background(cis, coding, lnc)):
        a = logcpm[[row[g] for g in pairs["lncrna"]]]
        b = logcpm[[row[g] for g in pairs["neighbour"]]]
        stats = neighbours.correlate(a, b, strata)
        parts.append(pairs.select(
            lncrna_id="lncrna", lncrna_name=pl.col("lncrna").replace_strict(name),
            neighbour_id="neighbour", neighbour_name=pl.col("neighbour").replace_strict(name),
            gap="gap", kind="kind",
        ).with_columns(**{k: pl.Series(v) for k, v in stats.items()}))  # fmt: skip
    return pl.concat(parts)


def kinds_table(pairs: pl.DataFrame) -> list[tuple[str, ...]]:
    """Per kind, background last: pairs, median rho on residuals, share called both ways."""
    rows: list[tuple[str, ...]] = []
    for kind in (*neighbours.KINDS, neighbours.BACKGROUND):
        found = pairs.filter(pl.col("kind") == kind)
        if not found.height:
            continue
        rho = found["rho_within"].to_numpy()
        called = (f"{100 * np.mean(found[q].to_numpy() < ALPHA):.0f} %" for q in ("q", "q_within"))
        rows.append((kind_words(kind), f"{found.height:,}", f"{np.nanmedian(rho):.2f}", *called))
    return rows


def kind_words(kind: str) -> str:
    """A pair kind in page words."""
    return text.fill(PAGE, f"kinds.{kind}", kb=neighbours.NEAR_NT / 1000)


def shown_pairs(pairs: pl.DataFrame, n: int) -> pl.DataFrame:
    """The cis pairs the select offers: smallest q on residuals, then largest |rho|."""
    cis = pl.col("kind") != neighbours.BACKGROUND
    return (pairs.filter(cis & pl.col("q_within").is_finite())
            .sort("q_within", -pl.col("rho_within").abs(), "lncrna_id").head(n))  # fmt: skip
