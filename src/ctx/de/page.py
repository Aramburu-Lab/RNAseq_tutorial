"""M5 page copy: the words around the numbers, each sentence filled from the run's own fits.

Nothing here names a cohort. A self-check tied to a feature of the design (a batch term, a
paired contrast) is asked only when the dataset file declares that feature.
"""

import math
import re
from collections.abc import Mapping, Sequence

from ctx.dataset import Contrast, Dataset
from ctx.de.deseq import TOP_GENES, Fit
from ctx.de.stats import ADJUSTMENTS, ALPHA, LFC_MIN
from ctx.report.builder import FigureText, SelfCheck
from ctx.viz.annotate import Permanova

QUESTION = "Which genes change between two groups, and how sure can we be?"
CONCEPT = (
    "DESeq2 asks, gene by gene, whether read counts differ between two groups by more than "
    "they vary within each group. The design names everything else that differs between "
    "samples, such as the individual or the sequencing batch, so that it is not mistaken for "
    "the group effect."
)
TERMS = (
    ("Size factor", "A per-sample scaling for sequencing depth, so counts compare across samples."),
    ("Dispersion", "How much a gene varies between samples beyond counting noise."),
    ("log2 fold change", "log2 of tested over reference: 1 is twice as high, -1 half as high."),
    ("Wald test", "Asks whether a gene's log2 fold change differs from 0, given its variation."),
    ("padj", "The p-value adjusted over all tested genes; DESeq2 uses BH."),
    ("BH, BY", "Benjamini-Hochberg and Benjamini-Yekutieli: they control the false discovery "
     "rate when the p-values are exact, BY under any dependence."),
    ("Holm, Bonferroni", "They control the chance of even one false call; stricter."),
    ("Called", f"padj below {ALPHA} and a log2 fold change of at least {LFC_MIN:g} either way."),
    ("Shrunk fold change", "An estimate pulled toward 0 where the data are thin; ranks genes."),
    ("PERMANOVA", "Tests whether groups of samples sit apart; R² is the share of variation "
     "they explain."),
)  # fmt: skip
PRACTICE = (
    "Check the p-value histogram before the volcano. Put a known batch in the design rather "
    "than dropping samples. Rank genes by the shrunk fold change; the test uses the raw one."
)
MISTAKES = (
    "Blaming batch for the gap between a small within-batch contrast and a large pooled one: "
    "the smaller contrast also has less power.",
    "Counting several samples from one individual as independent replicates: the test "
    "becomes overconfident.",
    "Reading the volcano's height as the call: it is the raw p-value, and the call uses the "
    "adjusted one.",
    # Measured in D-018 on tests/synth.py design_counts, seeds 0-199; not this run's numbers.
    "Reading padj < 0.05 as at most 5 % false genes: BH bounds the average share, and only "
    "if the p-values of unchanged genes are exact. DESeq2 takes each gene's estimated "
    "dispersion as known, so they run small. On synthetic data with 6 samples per group, "
    "8.4 % of the genes at padj < 0.05 were false on average, against the 4.5 % BH "
    "promised; the call's fold-change filter brought it to 3.2 %.",
)
NA_CHECK = SelfCheck(
    "Why can a gene with a large fold change have padj = NA?",
    "Either one sample is an outlier (Cook's distance), or the gene is so weakly expressed "
    "that independent filtering removed it before the adjustment.",
)
BATCH_ANSWER = (
    "When one group sits mostly in one batch, a design without batch reads the batch effect "
    "as a group effect. With batch in the design, DESeq2 estimates it from the groups present "
    "in both batches and removes it."
)
PAIRED_CHECK = SelfCheck(
    "Why does a paired design usually call more genes than an unpaired one on the same "
    "individuals?",
    "Each individual is compared with themself, so the differences between individuals leave "
    "the error term.",
)
# The words of each adjustment on the page, from the one list of methods; the two long
# names are abbreviated as the field does, and TERMS spells them out.
SHORT = {"Benjamini-Hochberg": "BH", "Benjamini-Yekutieli": "BY"}
BASES = {column: f"{SHORT.get(m, m)}-adjusted p" for column, m in ADJUSTMENTS.items()}
RAW_WORDS = "raw p (no correction)"
VOLCANO_READING = (
    "Each dot is one of {n_tested} {unit} tested in both groups; right of 0 is higher in "
    "{tested}, left is higher in {reference}.",
    "Up in {tested}: log2 FC at least {fc} ({fold_up} times as high, +{pct_up} %) and {line_rule}.",
    "Down in {tested}: log2 FC at most -{fc} ({fold_down} times as high, -{pct_down} %) under "
    "the same p-value rule. Grey dots miss at least one of the two conditions.",
)
VOLCANO_SUMMARY = (
    "Up in {tested}: {n_up} {unit}. Down: {n_down} {unit}. Rule: log2 FC at or beyond "
    "\u00b1{fc} and {basis} below {cutoff}."
)
LINE_RULE = "above the dashed line at -log10 p = {y_line}, where {basis} falls below {cutoff}"
NO_LINE_RULE = "{basis} below {cutoff}; no {unit} were tested, so no line is drawn"
# What one point of M5's volcano is.
UNIT = "genes"


def _terms(design: str) -> set[str]:
    """The terms of a design formula, e.g. ``{"batch", "group"}`` for ``~batch + group``."""
    return set(re.findall(r"[A-Za-z_]\w*", design))


def think(spec: Dataset, fits: Mapping[str, Fit]) -> tuple[SelfCheck, ...]:
    """The self-checks this dataset's design can answer, with its own names and counts."""
    checks = []
    assert spec.de is not None  # noqa: S101 - M5 runs only with [de]
    for names in spec.de.lessons.values():
        batch = [
            (a, b)
            for i, a in enumerate(names)
            for b in names[i + 1 :]
            if set(fits[a].samples["sample"]) == set(fits[b].samples["sample"])
            and _terms(fits[a].contrast.design) ^ _terms(fits[b].contrast.design) == {"batch"}
        ]
        if batch:
            a, b = batch[0]
            n = fits[a].samples.height
            checks.append(
                SelfCheck(
                    f"{a} and {b} fit the same {n} samples; only the design differs. Why can "
                    "they call different genes?",
                    BATCH_ANSWER,
                )
            )
            break
    if any(c.paired for c in spec.contrasts):
        checks.append(PAIRED_CHECK)
    return (*checks, NA_CHECK)


def _subset(filters: Mapping[str, Sequence[str]]) -> str:
    """The samples a figure uses, in words."""
    return "; ".join(f"{k.replace('_', ' ')} {', '.join(v)}" for k, v in filters.items())


def pca_text(name: str, filters: Mapping[str, Sequence[str]], factor: str, symbol: str | None,
             share: Sequence[float], test: tuple[str, Permanova] | None) -> FigureText:  # fmt: skip
    """Title, lead and reading guide of one PCA; ``test`` is the tested factor and its test.

    The colour and the test can name different factors: a factor with one level colours the
    points but cannot be tested, and the next factor is.
    """
    shown = factor.replace("_", " ")
    reading = [
        f"Each dot is one sample; samples close together have similar expression over the "
        f"{TOP_GENES} most variable genes.",
        f"Colour marks the {shown}"
        + (f", shape the {symbol.replace('_', ' ')}" if symbol else "")
        + f". PC1 holds {share[0]:.0%} and PC2 {share[1]:.0%} of those genes' variation.",
    ]
    lead = ""
    if test is not None:
        tested, result = test[0].replace("_", " "), test[1]
        p = "no p-value possible" if math.isnan(result.p) else f"p = {result.p:.3g}"
        reading.append(
            f"The {tested} accounts for {result.r2:.0%} of the variation between samples "
            f"(PERMANOVA R², {p}): the larger, the further apart the groups sit."
        )
        lead = f"The {tested} accounts for {result.r2:.0%} of the variation between samples."
    subset = _subset(filters)
    return FigureText(
        title=f"Samples{': ' + subset if subset else ''}, coloured by {shown}",
        lead=lead,
        reading=tuple(reading),
        tab=name,
    )


def volcano_text(contrast: Contrast) -> FigureText:
    """Title and lead of the volcano; its reading guide comes live from the figure."""
    return FigureText(
        title=f"Volcano: {contrast.tested} over {contrast.reference}",
        tab="Volcano",
    )


def ma_text(contrast: Contrast, shrunk: bool = True) -> FigureText:
    """Title, lead and reading guide of the MA plot; the raw estimate alone without shrinkage."""
    if not shrunk:
        return FigureText(
            title=f"Fold change against expression: {contrast.tested} over {contrast.reference}",
            lead="Shrinkage was off for this run (--no-lfc-shrink): the raw estimate alone.",
            reading=(
                "Each gene is one point: its raw log2 fold change against its mean expression.",
                "Weakly expressed genes spread widest, because few reads make a noisy fold "
                "change. Shrinkage, on by default, would pull them toward 0.",
                "The test uses this raw estimate either way.",
            ),
            tab="MA plot",
        )
    return FigureText(
        title=f"Fold change against expression: {contrast.tested} over {contrast.reference}",
        lead="Shrinkage tames the noisy fold changes of weakly expressed genes.",
        reading=(
            "Each gene appears twice against its mean expression: its raw log2 fold change and "
            "its shrunk one, as the legend marks.",
            "Shrinkage pulls weakly expressed genes, which have few reads to go on, toward 0; "
            "strongly expressed genes barely move.",
            "The test uses the raw estimate; rank genes by the shrunk one.",
        ),
        tab="MA plot",
    )


def phist_text(
    name: str, tested: int, low: int, flat: float, first: int, width: float
) -> FigureText:
    """Title, lead and reading guide of the p-value histogram.

    ``low`` counts p below ALPHA and ``first`` the first bar; they agree only when the bin
    width equals ALPHA, so each sentence names its own.
    """
    ratio = f"{first / flat:.1f} times" if flat else "far above"
    return FigureText(
        title=f"p-value histogram: {name}",
        lead=f"{low:,} genes have p < {ALPHA}, against a flat level of {flat:,.1f} per bin.",
        reading=(
            f"Each bar counts the genes whose p-value falls in a {width:g}-wide bin; "
            f"{tested:,} genes in all.",
            f"The dashed line is the flat level, {flat:,.1f} genes per bin: how p-values spread "
            "where nothing changes.",
            f"The first bar holds {first:,} genes, {ratio} the flat level. A peak at 0 over a "
            "flat floor is signal; a hill anywhere else is a model problem.",
        ),
        tab="p-value histogram",
    )


def dispersion_text(name: str) -> FigureText:
    """Title, lead and reading guide of the dispersion plot."""
    return FigureText(
        title=f"Dispersion: {name}",
        lead="Each gene's final dispersion is pulled toward the trend fitted across genes.",
        reading=(
            "Each gene appears twice against its mean expression: its own (gene-wise) "
            "dispersion and the final one the test uses.",
            "The line is the trend DESeq2 fits across genes; it pulls each final value toward "
            "it, more where the gene's own estimate is less certain.",
            "A cloud far above the trend points to outlier samples or a missing design term.",
        ),
        tab="Dispersion",
    )


def calls_text(title: str, names: Sequence[str], fits: Mapping[str, Fit]) -> FigureText:
    """Title, lead and reading guide of one lesson's call counts."""
    counts = {n: int((fits[n].results["call"] != "ns").sum()) for n in names}
    most = max(counts, key=lambda n: counts[n])
    return FigureText(
        title=f"Genes called: {title.lower()} lesson",
        lead=f"{most} calls the most genes: {counts[most]:,}.",
        reading=(
            f"Each bar is one contrast; its length counts the genes called up and down (padj "
            f"below {ALPHA}, log2 FC at or beyond ±{LFC_MIN:g}).",
            "Compare contrasts only at equal n: a smaller contrast calls fewer genes even when "
            "nothing else differs.",
            "The caption gives each contrast's n and design; the difference between designs "
            "on the same samples is the lesson.",
        ),
    )


def lfc_text(x: str, y: str, rho: float) -> FigureText:
    """Title, lead and reading guide of one contrast's fold changes against another's."""
    return FigureText(
        title=f"{x} against {y}",
        lead=f"The two contrasts agree with Spearman rho = {rho:.2f}.",
        reading=(
            f"Each dot is one gene tested in both; it sits at its raw log2 fold change in {x} "
            f"(across) and in {y} (up).",
            "Dots on the solid diagonal changed the same way in both; the tighter the cloud "
            "around it, the more the two contrasts agree.",
            "Highlighted dots are called in at least one of the two.",
        ),
    )


def box_text(contrast: Contrast, n_genes: int) -> FigureText:
    """Title, lead and reading guide of the named-gene boxes."""
    return FigureText(
        title=f"Named genes: {contrast.tested} against {contrast.reference}",
        reading=(
            "Each box is one gene in one group: the line is the median sample, the box holds "
            "the middle half, and each dot is a sample.",
            "Expression is log2 of the count divided by the sample's size factor, plus 1.",
            f"Above each gene: the Mann-Whitney p, BH-adjusted over the {n_genes} plotted "
            f"genes, and r, the rank-biserial effect, positive when {contrast.tested} is higher.",
        ),
    )
