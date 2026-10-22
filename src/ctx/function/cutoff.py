"""M6's FDR cutoff (P2S.5, D-109): a page control on each figure that marks gene sets by FDR.

Three figures mark sets by FDR: the over-representation bars, whose page reruns ``ctxOra`` at
the cutoff, the GSEA bars and one set's NES in every contrast. The background lesson colours
its bars by background and keeps the build's FDR in its lead. One rule: a set passes at an FDR
below the cutoff, strictly, as the build marks it at ``ALPHA``; the presets are the volcano's,
``scatter.CUTOFF_PRESETS``, ``ALPHA`` among them and selected. ``templates/_fdr.js.j2`` is the
page twin: ``ctxFdr.mark`` splits the drawn bars at another cutoff as ``BarData`` draws them
there, and ``ctxFdr.fill`` fills the words, each a template of the cutoff from ``content/``.
"""

import math
from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import Any

import polars as pl

from ctx.de.stats import ALPHA
from ctx.report import export, marks, text
from ctx.report.builder import FigureBlock
from ctx.viz.bar import BarData
from ctx.viz.scatter import CUTOFF_PRESETS

PAGE = "function"
# As JavaScript's String() writes each preset, the select's values.
PRESETS = tuple(f"{c:g}" for c in CUTOFF_PRESETS)
# What a word holds in place of the cutoff, which the page fills (ctxFdr.fill).
HOLE = "{alpha}"


def significance(fdr: pl.Expr, alpha: float = ALPHA) -> pl.Expr:
    """A set's mark on the GSEA and NES bars: the passing words below ``alpha``, else failing."""
    passing = text.fill(PAGE, "significance.passing", alpha=alpha)
    return pl.when(fdr < alpha).then(pl.lit(passing)).otherwise(
        pl.lit(text.word(PAGE, "significance.failing")))  # fmt: skip


def _levels(bars: BarData) -> list[dict[str, Any]]:
    """Each series' trace style, passing then failing, the cutoff in the passing words a hole."""
    built = text.fill(PAGE, "significance.passing", alpha=ALPHA)
    hole = text.fill(PAGE, "significance.passing", alpha=HOLE)
    styles = [bars.trace_style(level) for level in bars.colours]
    for style in styles:
        style["name"] = style["name"].replace(built, hole)
        style["hovertemplate"] = style["hovertemplate"].replace(built, hole)
    return styles


def marking(
    bars: BarData,
    reading: tuple[int, Mapping[str, str]],
    fdrs: Sequence[Sequence[float | None]] = (),
    lead: str = "",
) -> dict[str, Any]:
    """The cutoff's settings for bars marked by ``significance``: ``ctxFdr.mark`` reads them.

    Args:
        bars: The figure as built, its hover carrying ``fdr``.
        reading: The reading line the cutoff rewrites, by position, and its words, filled but
            for ``HOLE``, when ``all``, ``no`` or ``some`` of the bars drawn pass.
        fdrs: Every tested set's FDR, per view of the figure in view order, which the lead counts.
        lead: The lead's key in ``figures``, filled with the count passing of those tested.
    """
    index, words = reading
    return {
        "presets": PRESETS, "alpha": ALPHA, "at": list(bars.hover).index("fdr"),
        "horizontal": bars.horizontal, "levels": _levels(bars),
        "label": text.word(PAGE, "cutoff.label"), "moved": text.raw(PAGE, "cutoff.moved"),
        "noscript": text.word(PAGE, "cutoff.noscript"),
        "reading": {"index": index, **{w: words[w] for w in ("all", "no", "some")}},
        "fdrs": [list(f) for f in fdrs],
        "lead": text.raw(PAGE, f"figures.{lead}") if lead else "",
    }  # fmt: skip


def gsea(bars: BarData, tables: Sequence[pl.DataFrame]) -> dict[str, Any]:
    """The GSEA bars' cutoff: the lead counts each view's tested sets, the last line the bars.

    ``tables`` are the GSEA tables in the view's order, the build's first.
    """
    words = {w: text.fill(PAGE, f"figures.gsea.{w}", alpha=HOLE) for w in ("all", "no", "some")}
    # NaN has no JSON; None fails every cutoff on the page, as NaN fails < in the build.
    fdrs = [[None if f is None or math.isnan(f) else f for f in t["fdr"].to_list()] for t in tables]
    lines = len(text.raw(PAGE, "figures.gsea.reading"))
    return marking(bars, (lines, words), fdrs, "gsea.lead")


def nes(bars: BarData) -> dict[str, Any]:
    """One set's NES bars' cutoff: their last reading line names the highlight's cutoff."""
    lines = text.raw(PAGE, "figures.nes.reading")
    return marking(bars, (len(lines) - 1, dict.fromkeys(("all", "no", "some"), lines[-1])))


def attach(block: FigureBlock, bars: BarData, spec: Mapping[str, Any]) -> FigureBlock:
    """``block`` with the cutoff: beside a view's select, or as the figure's own control.

    Its own control carries the export settings a moved control needs (D-033), as
    ``save.figure_block`` gives a contract's controls.
    """
    if block.controls is not None and block.controls.get("kind") == "view":
        return replace(block, controls={**block.controls, "fdr": dict(spec)})
    if block.controls is not None:
        raise ValueError(f"{block.stem}: the cutoff joins a view or stands alone")
    settings = export.settings(bars.inches(), bars.caption, block.stem, bars.export_layers())
    return replace(block, controls={"kind": "fdr", "fdr": dict(spec), "export": settings})


def ora(caption: str, background: int, reading: int) -> dict[str, Any]:
    """The live ORA's cutoff settings: its words, each a template of the cutoff.

    Args:
        caption: The figure's caption, filled but for ``HOLE``.
        background: The build's tested gene symbols, which the passing reading line counts.
        reading: The position of the reading line that names the dashed line.
    """

    def ora_word(key: str) -> str:
        return str(text.raw(PAGE, f"figures.ora.{key}"))

    return {
        "presets": PRESETS, "alpha": ALPHA, "label": text.word(PAGE, "cutoff.label"),
        "line": ora_word("line"), "caption": caption,
        "export_caption": marks.plain(export.file_caption(caption)),
        "leads": {"each": ora_word("lead"), **text.raw(PAGE, "live.leads")},
        "reading": {"index": reading, "passing": ora_word("passing"),
                    "failing": ora_word("failing"), "background": background},
    }  # fmt: skip
