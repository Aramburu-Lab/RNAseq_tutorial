"""M6 page copy: the words around the enrichment figures, from ``content/function.toml``.

Code names a key and fills it from the run's own results; it writes no sentence. A contrast
reads by its plain label and a gene set in words (Rytis 2026-10-01).
"""

from typing import Any

from ctx.de.page import Words
from ctx.de.stats import ALPHA, LFC_MIN
from ctx.function import lite, live
from ctx.report import text
from ctx.report.builder import FigureText, SelfCheck

PAGE = "function"
CALLS = text.fill(PAGE, "figures.ora.calls", alpha=ALPHA, lfc_min=f"{LFC_MIN:g}")


def figure(name: str, key: str, **values: Any) -> str:
    """One figure's word at ``figures.<name>.<key>``, filled."""
    return text.fill(PAGE, f"figures.{name}.{key}", **values)


def set_words(name: str) -> str:
    """A gene set in words: "Interferon alpha response, MSigDB Hallmark".

    The identifier's name after its collection prefix, in sentence case, then the collection
    (``[gene_sets.collections]``); an identifier with no listed prefix reads as itself in words.
    """
    prefix, _, rest = name.partition("_")
    collections = text.raw(PAGE, "gene_sets.collections")
    plain = (rest if prefix in collections else name).replace("_", " ").lower()
    plain = plain[:1].upper() + plain[1:]
    if prefix not in collections:
        return plain
    return text.fill(PAGE, "gene_sets.named", name=plain, collection=collections[prefix])


def fdr_words(fdr: float) -> str:
    """An FDR as the page prints it; 0 means no null score was as extreme, not certainty."""
    return text.word(PAGE, "fdr.zero") if fdr == 0 else f"{fdr:.2g}"


def ora_text(said: Words, background: int, up: int, down: int, failing: int = 0,
             alpha: float = ALPHA) -> FigureText:  # fmt: skip
    """Title, lead and reading guide of the ORA bars (``live.bars``), at FDR cutoff ``alpha``.

    ``up`` and ``down`` count the sets passing on each direction's calls; the lead states
    both, so a direction with no passing set says it has none. ``failing`` counts the grey
    bars, which the last reading line names when there are any (PK.3).
    """
    last = (figure("ora", "failing", alpha=alpha, failing=failing, top=live.TOP) if failing
            else figure("ora", "passing", alpha=alpha, background=background))  # fmt: skip
    return FigureText(
        title=figure("ora", "title", label=said.label),
        lead=figure("ora", "lead", up=up, down=down, alpha=alpha),
        reading=(*text.lines(PAGE, "figures.ora.reading", tested=said.tested, calls=CALLS), last),
        tab=figure("ora", "tab"),
    )


def gsea_text(said: Words, passing: int, tested: int, drawn: tuple[int, int],
              alpha: float = ALPHA) -> FigureText:  # fmt: skip
    """Title, lead and reading guide of the GSEA bars, marked at FDR cutoff ``alpha``.

    ``drawn`` is (bars that pass, bars): the last line names what the highlight separates
    among the bars drawn. On Step-1 all 15 drawn passed, and "the others" pointed at none.
    """
    shown, bars = drawn
    which = "all" if shown == bars else "no" if shown == 0 else "some"
    return FigureText(
        title=figure("gsea", "title", label=said.label),
        lead=figure("gsea", "lead", passing=passing, tested=tested, alpha=alpha),
        reading=(
            *text.lines(PAGE, "figures.gsea.reading", tested=said.tested, reference=said.reference),
            figure("gsea", which, alpha=alpha),
        ),  # fmt: skip
        tab=figure("gsea", "tab"),
    )


def nes_text(gene_set: str, contrasts: int, alpha: float = ALPHA) -> FigureText:
    """Title, lead and reading guide of one set's NES across contrasts, marked at ``alpha``."""
    shown = set_words(gene_set)
    return FigureText(
        title=figure("nes", "title", gene_set=shown),
        lead=figure("nes", "lead", contrasts=contrasts),
        reading=text.lines(PAGE, "figures.nes.reading", alpha=alpha),
    )


def running_text(said: Words, result: lite.Enriched, ranked: int) -> FigureText:
    """Title, lead and reading guide of one running enrichment score."""
    return FigureText(
        title=figure("running", "title", gene_set=set_words(result.name), label=said.label),
        lead=figure("running", "lead", es=result.es, nes=result.nes, fdr=fdr_words(result.fdr)),
        reading=text.lines(
            PAGE,
            "figures.running.reading",
            ranked=ranked,
            tested=said.tested,
            reference=said.reference,
            size=result.size,
        ),  # fmt: skip
    )


def union_direction(more: int, less: int) -> str:
    """Which way the gene-set union moved most plotted sets: "more" or "less" significant.

    Against the union a set's share of the background grows by its own untested members and
    shrinks by the others', so the direction is measured, never assumed: on Step-1 ndn_gcsf
    6 of 7 plotted sets looked less significant. "mixed" when as many moved each way.
    """
    return "more" if more > less else "less" if less > more else "mixed"


def _check(key: str) -> SelfCheck:
    """One self-check from ``[questions.<key>]``."""
    return SelfCheck(text.word(PAGE, f"questions.{key}.prompt"),
                     text.word(PAGE, f"questions.{key}.answer"))  # fmt: skip


def questions(direction: str | None) -> tuple[SelfCheck, ...]:
    """The self-checks, with the background question where the two backgrounds differ.

    ``direction`` is ``union_direction``'s, and None where both hold the same genes.
    """
    rank = _check("rank")
    return (rank,) if direction is None else (_check(f"union_{direction}"), rank)


def trap_text(tested: int, union: int, moved: tuple[int, int], shown: int,
              passing: tuple[int, int]) -> FigureText:  # fmt: skip
    """Title, lead and reading guide of the background lesson, from the bars drawn.

    ``moved`` counts the plotted sets that look more and less significant against the wider
    background, ``union`` its size; ``passing`` the sets passing FDR < ALPHA against the tested
    genes and against it. The second reading line explains the direction most of them took.
    """
    more, less = moved
    lead = (
        figure("background", "lead_same")
        if union == tested
        else figure("background", "lead_moved", union=union, tested=tested, more=more,
                    less=less, shown=shown)
    )  # fmt: skip
    lead = " ".join([lead, figure("background", "lead_passing", right=passing[0],
                                  wrong=passing[1], alpha=ALPHA)])  # fmt: skip
    return FigureText(
        title=figure("background", "title"),
        lead=lead,
        reading=(
            figure("background", "reading_first"),
            figure("background", union_direction(more, less)),
            figure("background", "reading_last"),
        ),  # fmt: skip
    )
