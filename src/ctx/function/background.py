"""M6's background lesson with every background a reader may pick (WO7, D-051, D-073).

The same calls are tested by ``lite.ora`` against the genes DESeq2 tested and against a wider
background: the gene-set union with the tested genes, gseapy's default when given none (the
build's), every protein-coding gene symbol of the annotation, or every gene symbol of it. The
page looks each one up; nothing is recomputed there. On Step-1 ndn_gcsf's up calls the sets
that pass run 16 against the tested genes, 6 against the union, 19 against protein-coding and
44 against every annotated gene (figure controls plan), so a wrong background moves p either way.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import polars as pl

from ctx.de.page import Words
from ctx.de.stats import ALPHA
from ctx.function import lite
from ctx.function.page import figure, set_words
from ctx.report import text
from ctx.viz.bar import BarData
from ctx.viz.style import category_colours

PAGE = "function"
# The backgrounds beside the tested genes, the build's first; the calls; the sets drawn.
BACKGROUNDS = ("union", "protein_coding", "annotated")
DIRECTIONS = ("up", "down")
SHOWN = ("top", "all")
# The build draws the 7 sets of smallest p against the tested genes (TOP_SETS // 2).
TOP = 7
CODING = "protein_coding"


def backgrounds(tested: Sequence[str], sets: Mapping[str, Sequence[str]],
                types: pl.DataFrame | None) -> dict[str, list[str]]:  # fmt: skip
    """Each wider background's symbols; the annotation's two where its gene types were read.

    Args:
        tested: The tested genes' symbols.
        sets: The gene sets.
        types: The annotation's gene records (gene_name, gene_type), or None.
    """
    found = {"union": sorted({g for gs in sets.values() for g in gs} | set(tested))}
    if types is not None:
        named = types.filter(pl.col("gene_name").is_not_null() & (pl.col("gene_name") != ""))
        found[CODING] = sorted(set(named.filter(pl.col("gene_type") == CODING)["gene_name"]))
        found["annotated"] = sorted(set(named["gene_name"]))
    return found


@dataclass(frozen=True)
class Tested:
    """One direction's calls against the tested genes and against one wider background."""

    direction: str
    background: str
    size: int
    hits: int
    right: Mapping[str, lite.Overrepresented]
    wrong: Mapping[str, lite.Overrepresented]

    def shown(self, how: str) -> list[str]:
        """The sets drawn: tested in both, by p against the tested genes, the top ``TOP``.

        A tie keeps the gene-set file's order, as the build's bars always took it.
        """
        both = sorted((s for s in self.right if s in self.wrong), key=lambda s: self.right[s].p)
        return both[:TOP] if how == "top" else both

    def moved(self, how: str) -> tuple[int, int, int]:
        """Sets drawn that look more and less significant against the wider background; all.

        Counted on FDR, the bars' length, so the words match what the bars show.
        """
        drawn = self.shown(how)
        more = sum(self.wrong[s].fdr < self.right[s].fdr for s in drawn)
        less = sum(self.wrong[s].fdr > self.right[s].fdr for s in drawn)
        return more, less, len(drawn)

    def passing(self) -> tuple[int, int]:
        """Sets passing FDR < ALPHA against the tested genes and against the wider one."""
        return (sum(r.fdr < ALPHA for r in self.right.values()),
                sum(r.fdr < ALPHA for r in self.wrong.values()))  # fmt: skip


def tested_against(hits: Sequence[str], tested: Sequence[str], wide: Sequence[str],
                   sets: Mapping[str, Sequence[str]], direction: str,
                   background: str) -> Tested:  # fmt: skip
    """``lite.ora`` of the same hits against the tested genes and against ``wide``."""
    right = {r.name: r for r in lite.ora(hits, tested, sets)}
    wrong = {r.name: r for r in lite.ora(hits, wide, sets)}
    return Tested(direction, background, len(wide), len(hits), right, wrong)


def words(background: str) -> str:
    """A background as the bars' legend names it."""
    return figure("background", background)


def option(background: str, direction: str, how: str) -> str:
    """How the live line names one view."""
    return figure("background", "option", background=words(background),
                  direction=figure("background", f"direction_{direction}"),
                  shown=figure("background", f"shown_{how}", top=TOP))  # fmt: skip


def axes(kept: Sequence[str]) -> list[tuple[str, list[str]]]:
    """The view's selects: the wider background, the calls and the sets drawn, in words."""
    return [(figure("background", "axis_background"), [words(b) for b in kept]),
            (figure("background", "axis_direction"),
             [figure("background", f"direction_{d}") for d in DIRECTIONS]),
            (figure("background", "axis_shown"),
             [figure("background", f"shown_{h}", top=TOP) for h in SHOWN])]  # fmt: skip


def bars(said: Words, calls: int, tested: int, found: Tested, how: str) -> BarData:
    """One view's bars: each set drawn twice, against the tested genes and the wider one.

    The length is -log10 FDR with the dashed line at ``ALPHA``, as on the ORA bars: BH gives each
    background its own p cutoff, so no line on p could mark the sets that pass (Rytis 2026-10-26).
    """
    drawn = found.shown(how)
    labels = (figure("background", "tested_genes"), words(found.background))
    rows = [(set_words(s), label, float(-np.log10(r[s].fdr)))
            for s in drawn
            for label, r in zip(labels, (found.right, found.wrong), strict=True)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["gene_set", "background", "neg_log10_fdr"], orient="row")
    columns = {
        "gene_set": ("none", "Gene set in words."),
        "background": ("none", "ORA background."),
        "neg_log10_fdr": (
            "-log10 FDR",
            f"BH FDR over the sets of the {found.direction} calls' hypergeometric p.",
        ),
    }
    return BarData(
        frame=frame, columns=columns, category="gene_set", value="neg_log10_fdr",
        category_label=figure("background", "x"), value_label=figure("background", "y"),
        series="background",
        colours=dict(zip(labels, category_colours(list(labels)).values(), strict=True)),
        caption=text.fill(PAGE, "figures.background.caption", label=said.label, calls=calls,
                          direction=figure("background", f"direction_{found.direction}"),
                          hits=found.hits, tested=tested, union=found.size, alpha=ALPHA,
                          wider=figure("background", f"wider_{found.background}")),
        horizontal=True, width="double",
        hline=(-np.log10(ALPHA), figure("ora", "line", alpha=ALPHA)),
    )  # fmt: skip


TABLE = {
    "direction": ("none", "The calls tested: up or down."),
    "background": ("none", "The wider background: union, protein_coding or annotated."),
    "background_size": ("symbols", "Gene symbols in the wider background."),
    "gene_set": ("none", "Gene set."),
    "p_tested": ("dimensionless", "Hypergeometric p against the tested genes."),
    "fdr_tested": ("dimensionless", "BH FDR over the sets, against the tested genes."),
    "p_background": ("dimensionless", "Hypergeometric p against the wider background."),
    "fdr_background": ("dimensionless", "BH FDR over the sets, against the wider background."),
}  # fmt: skip


def table(found: Sequence[Tested]) -> pl.DataFrame:
    """Every p and FDR the views draw or count, one row per direction, background and set."""
    rows = [(t.direction, t.background, t.size, s, t.right[s].p, t.right[s].fdr,
             t.wrong[s].p, t.wrong[s].fdr)
            for t in found for s in t.shown("all")]  # fmt: skip
    return pl.DataFrame(rows, schema=list(TABLE), orient="row")
