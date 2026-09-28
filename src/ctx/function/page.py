"""M6 page copy: the words around the enrichment figures, filled from the run's own results."""

from ctx.dataset import Contrast
from ctx.de.stats import ALPHA, LFC_MIN
from ctx.function import lite
from ctx.report.builder import FigureText, SelfCheck

QUESTION = "Which pathways move together, beyond single genes?"
CONCEPT = (
    "Over-representation analysis (ORA) asks whether a pathway holds more of the called genes "
    "than chance predicts, counted against the genes that could have been called. Gene set "
    "enrichment analysis (GSEA) needs no cutoff: it ranks every gene with a Wald statistic and "
    "asks whether a pathway's genes gather at one end of the list."
)
TERMS = (
    ("Gene set", "A named list of genes that share a pathway or a function."),
    ("Background", "The genes that could have been called: the ones DESeq2 tested."),
    ("NES", "Normalised enrichment score: positive when the set gathers among genes that "
     "went up."),
    ("ORA FDR", "Benjamini-Hochberg over the hypergeometric p-values of the sets tested."),
    ("GSEA FDR", "The false discovery rate over all sets tested, from gene-set permutations."),
    ("Wald statistic", "A gene's fold change divided by its standard error; it carries the "
     "sign."),
)  # fmt: skip
PRACTICE = (
    "Use the tested genes as the background, rank by a signed statistic, and read each NES "
    "with its FDR. A pathway can move with no single gene called."
)
UNION_CHECK = SelfCheck(
    "Why does the gene-set union background make the same hits look more significant?",
    "It adds genes that could never have been called, such as ones not expressed in the cells, "
    "so the expected overlap shrinks and ordinary overlaps look rare.",
)
THINK = (
    SelfCheck(
        "Why rank by the Wald statistic rather than by the fold change?",
        "The statistic weighs the change by its certainty; a large fold change on few counts "
        "would otherwise sit at the top.",
    ),
)
MISTAKES = (
    "Using the whole genome, or every gene in the gene sets, as the ORA background: only "
    "tested genes could have been called.",
    "Reading the smallest GSEA p-value as exact: it cannot fall below one over the null "
    "scores of the same sign plus one, so it is a bound.",
    "Ranking genes for GSEA by padj: it has no sign, and many genes share its smallest values.",
)
CALLS = f"padj below {ALPHA}, log2 FC at or beyond ±{LFC_MIN:g}"


def fdr_words(fdr: float) -> str:
    """An FDR as the page prints it; 0 means no null score was as extreme, not certainty."""
    return "0 (no null score as extreme)" if fdr == 0 else f"{fdr:.2g}"


def ora_text(contrast: Contrast, background: int, passing: int) -> FigureText:
    """Title, lead and reading guide of the ORA bars."""
    return FigureText(
        title=f"Over-representation: {contrast.tested} over {contrast.reference}",
        lead=f"{passing} gene sets pass FDR < {ALPHA} on the up or the down calls.",
        reading=(
            "Each bar is one gene set; its length is -log10 of its FDR, so a longer bar is "
            "stronger evidence of over-representation.",
            f"Colour says which calls were counted: up or down in {contrast.tested} ({CALLS}).",
            f"A bar past the dashed line has FDR below {ALPHA}: more calls in the set than the "
            f"{background:,} tested genes predict.",
        ),
        tab="ORA",
    )


def gsea_text(contrast: Contrast, passing: int, tested: int) -> FigureText:
    """Title, lead and reading guide of the GSEA bars."""
    return FigureText(
        title=f"Gene set enrichment: {contrast.tested} over {contrast.reference}",
        lead=f"{passing} of {tested} tested gene sets pass FDR < {ALPHA}.",
        reading=(
            "Each bar is one gene set, the sets with the smallest FDR on top; its length is the "
            "NES.",
            f"Right of 0, the set's genes gather among genes higher in {contrast.tested}; left "
            f"of 0, among genes higher in {contrast.reference}.",
            f"Highlighted bars pass FDR < {ALPHA}; the others are shown for context.",
        ),
        tab="GSEA",
    )


def nes_text(gene_set: str, contrasts: int) -> FigureText:
    """Title, lead and reading guide of one set's NES across contrasts."""
    return FigureText(
        title=f"{gene_set} in every contrast",
        lead=f"The same set, ranked and permuted separately in each of {contrasts} contrasts.",
        reading=(
            "Each bar is one contrast; its length is the set's NES there.",
            "Right of 0, the set's genes sit among genes higher in that contrast's tested "
            "group; left of 0, lower.",
            f"Highlighted bars pass FDR < {ALPHA}.",
        ),
    )


def running_text(contrast: Contrast, result: lite.Enriched, ranked: int) -> FigureText:
    """Title, lead and reading guide of one running enrichment score."""
    return FigureText(
        title=f"Running score: {result.name} in {contrast.name}",
        lead=f"ES = {result.es:.2f}, NES = {result.nes:.2f}, FDR = {fdr_words(result.fdr)}.",
        reading=(
            f"The {ranked:,} genes are ranked from most up to most down in {contrast.tested}; "
            "the curve rises at each set member and falls between them.",
            f"Dots mark the set's {result.size} members; the curve's furthest point from 0 is "
            "the enrichment score.",
            f"A peak near the start means the members gather among genes higher in "
            f"{contrast.tested}; a trough near the end, among genes higher in "
            f"{contrast.reference}.",
        ),
    )


def think(tested: int, union: int) -> tuple[SelfCheck, ...]:
    """The self-checks; the background question only where the two backgrounds differ."""
    return (UNION_CHECK, *THINK) if union > tested else THINK


def trap_text(tested: int, union: int, more: int, shown: int) -> FigureText:
    """Title, lead and reading guide of the background lesson, from the bars drawn.

    ``more`` counts the plotted sets that look more significant against the union.
    """
    lead = (
        "Here both backgrounds hold the same genes, so the bars match."
        if union == tested
        else f"Against {union:,} genes instead of {tested:,}, {more} of {shown} sets look more "
        "significant."
    )
    return FigureText(
        title="The background decides the p-value",
        lead=lead,
        reading=(
            "Each gene set appears twice: tested against the genes DESeq2 tested, and against "
            "every gene in any set plus the tested ones, gseapy's default.",
            "A longer bar looks more significant; genes that could never be called make the "
            "same overlap look rarer.",
            "Use the tested genes: only they could have been called.",
        ),
    )
