"""M6 page copy: the words around the enrichment figures, filled from the run's own results.

The sentences moved so far are in ``content/function.toml``; the rest move with P2S.4.
"""

from ctx.dataset import Contrast
from ctx.de.stats import ALPHA, LFC_MIN
from ctx.function import lite
from ctx.report.builder import FigureText, SelfCheck
from ctx.report.text import fill, word

PAGE = "function"
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


def ora_text(contrast: Contrast, background: int, up: int, down: int) -> FigureText:
    """Title, lead and reading guide of the ORA bars, which draw only the sets that pass.

    ``up`` and ``down`` count the sets passing on each direction's calls; the lead states
    both, so a direction with no bar says it has none.
    """
    return FigureText(
        title=f"Over-representation: {contrast.tested} over {contrast.reference}",
        lead=fill(PAGE, "figures.ora.lead", up=up, down=down, alpha=ALPHA),
        reading=(
            "Each bar is one gene set; its length is -log10 of its FDR, so a longer bar is "
            "stronger evidence of over-representation.",
            f"Colour says which calls were counted: up or down in {contrast.tested} ({CALLS}).",
            fill(PAGE, "figures.ora.passing", alpha=ALPHA, background=background),
        ),
        tab="ORA",
    )


def gsea_text(contrast: Contrast, passing: int, tested: int, drawn: tuple[int, int]) -> FigureText:
    """Title, lead and reading guide of the GSEA bars.

    ``drawn`` is (bars that pass, bars): the last line names what the highlight separates
    among the bars drawn. On Step-1 all 15 drawn passed, and "the others" pointed at none.
    """
    shown, bars = drawn
    which = "all" if shown == bars else "no" if shown == 0 else "some"
    return FigureText(
        title=f"Gene set enrichment: {contrast.tested} over {contrast.reference}",
        lead=f"{passing} of {tested} tested gene sets pass FDR < {ALPHA}.",
        reading=(
            "Each bar is one gene set, the sets with the smallest FDR on top; its length is the "
            "NES.",
            f"Right of 0, the set's genes gather among genes higher in {contrast.tested}; left "
            f"of 0, among genes higher in {contrast.reference}.",
            fill(PAGE, f"figures.gsea.{which}", alpha=ALPHA),
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


def union_direction(more: int, less: int) -> str:
    """Which way the gene-set union moved most plotted sets: "more" or "less" significant.

    Against the union a set's share of the background grows by its own untested members and
    shrinks by the others', so the direction is measured, never assumed: on Step-1 ndn_gcsf
    6 of 7 plotted sets looked less significant. "mixed" when as many moved each way.
    """
    return "more" if more > less else "less" if less > more else "mixed"


def think(direction: str | None) -> tuple[SelfCheck, ...]:
    """The self-checks, with the background question where the two backgrounds differ.

    ``direction`` is ``union_direction``'s, and None where both hold the same genes.
    """
    if direction is None:
        return THINK
    key = f"questions.union_{direction}"
    return (SelfCheck(word(PAGE, f"{key}.prompt"), word(PAGE, f"{key}.answer")), *THINK)


def trap_text(tested: int, union: int, moved: tuple[int, int], shown: int) -> FigureText:
    """Title, lead and reading guide of the background lesson, from the bars drawn.

    ``moved`` counts the plotted sets that look more and less significant against the union;
    the second reading line explains the direction most of them took.
    """
    more, less = moved
    lead = (
        word(PAGE, "figures.background.lead_same")
        if union == tested
        else fill(PAGE, "figures.background.lead_moved", union=union, tested=tested, more=more,
                  less=less, shown=shown)
    )  # fmt: skip
    return FigureText(
        title="The background decides the p-value",
        lead=lead,
        reading=(
            "Each gene set appears twice: tested against the genes DESeq2 tested, and against "
            "every gene in any set plus the tested ones, gseapy's default.",
            word(PAGE, f"figures.background.{union_direction(more, less)}"),
            "Use the tested genes: only they could have been called.",
        ),
    )
