"""M6's running enrichment score for any set the contrast tested (WO7, D-051, D-012).

The page draws a set's running sum with ``ctxRunning`` (``templates/_running.js.j2``), the twin
of ``lite._running`` tested against it under node; NES, p and FDR, which need the permutations,
are looked up, as are the words. The page carries the ranking's weights in rank order once, as
float64, and each set as its members' ranks. Its ES as the page computes it is checked against
the GSEA table's, and the live line says where they differ. A second select draws another
contrast's ranking and sets, each carried as the first is; changing it keeps the set shown where
that contrast tested it (``keep``, whose page twin is ``ctxRunning.keep``).
"""

import base64
from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import Any

import numpy as np
import polars as pl

from ctx.de.page import Words
from ctx.de.run import Writer
from ctx.function import lite
from ctx.function.page import fdr_words, figure, running_text, set_words
from ctx.report import export, marks, text
from ctx.report.builder import FigureBlock
from ctx.viz.scatter import Line, ScatterData
from ctx.viz.style import ROLES

# One contrast's words, M5 results and GSEA results by set.
Contrasted = tuple[Words, pl.DataFrame, Mapping[str, lite.Enriched]]
# Above a float64 running sum's drift over a ranking (numpy sums the weights pairwise, the page
# from the first rank), below the third decimal any ES is printed to.
ES_TOLERANCE = 1e-9


def ranking(results: pl.DataFrame) -> tuple[list[str], np.ndarray]:
    """Symbols and Wald statistics, one gene per symbol (highest baseMean)."""
    rows = results.filter(pl.col("stat").is_not_null())
    return lite.collapse(
        rows["gene_id"].to_list(),
        rows["gene_name"].fill_null("").to_list(),
        rows["stat"].to_list(),
        rows["baseMean"].to_list(),
    )


def caption(said: Words, result: lite.Enriched, ranked: int) -> str:
    """The running score's caption: the set's statistics and the permutations behind them."""
    return figure("running", "caption", gene_set=set_words(result.name), label=said.label,
                  es=result.es, nes=result.nes, p=result.p, fdr=fdr_words(result.fdr),
                  size=result.size, ranked=ranked, permutations=lite.PERMUTATIONS,
                  seed=lite.SEED)  # fmt: skip


def weights(results: pl.DataFrame) -> np.ndarray:
    """The weight of each rank of the contrast's ranking, as ``lite.prerank`` walks it."""
    return lite.ranking(*ranking(results))[1]


def order(enriched: Mapping[str, lite.Enriched], shown: str) -> list[str]:
    """A contrast's tested sets as the select lists them, ``shown`` first where it was tested.

    The rest follow by their names in words.
    """
    return [*([shown] if shown in enriched else []),
            *sorted((s for s in enriched if s != shown), key=set_words)]  # fmt: skip


def keep(names: Sequence[str], current: str, build: str) -> int:
    """The set a contrast change shows, as an index into that contrast's ``names``.

    The one shown where that contrast tested it, else the dataset's running set, else its
    first. ``ctxRunning.keep`` is its page twin.
    """
    for wanted in (current, build):
        if wanted in names:
            return list(names).index(wanted)
    return 0


def _ranking(said: Words, results: pl.DataFrame, enriched: Mapping[str, lite.Enriched],
             shown: str) -> dict[str, Any]:  # fmt: skip
    """One contrast's ranking in base64 float64 and its tested sets, each with its own words."""
    weight = weights(results)
    sets = []
    for name in order(enriched, shown):
        hit = enriched[name]
        words = running_text(said, hit, len(weight))
        written = caption(said, hit, len(weight))
        sets.append({"label": set_words(name), "name": name, "title": words.title,
                     "lead": str(marks.html(words.lead)), "caption": str(marks.html(written)),
                     "export_caption": marks.plain(written), "es": hit.es,
                     "positions": list(hit.positions)})  # fmt: skip
    encoded = base64.b64encode(np.asarray(weight, dtype="<f8").tobytes()).decode("ascii")
    return {"label": said.label, "sets": sets, "weights": encoded}


def panel(said: Words, results: pl.DataFrame, enriched: Mapping[str, lite.Enriched],
          shown: str, exported: Mapping[str, Any],
          others: Sequence[Contrasted] = ()
          ) -> dict[str, Any]:  # fmt: skip
    """The running score's control: the ranking and every tested set, the build's first.

    Args:
        said: The contrast in words.
        results: Its M5 results.
        enriched: Its GSEA results by set.
        shown: The set the figure draws as built.
        exported: The figure's export settings; another set's export changes the caption.
        others: Every other contrast's words, results and GSEA results, in page order; each
            becomes an option of the contrast select.

    Returns:
        JSON-ready settings for the page script; ``contrasts`` holds the build's contrast
        first, whose ranking and sets are the top level's.
    """
    built = _ranking(said, results, enriched, shown)
    rest = [_ranking(w, r, e, shown) for w, r, e in others if e]
    return {"kind": "running", "label": figure("running", "view"), "sets": built["sets"],
            "weights": built["weights"], "tolerance": ES_TOLERANCE,
            "contrast_label": figure("running", "contrast"),
            "contrasts": [{"label": said.label}, *rest], "build_set": shown,
            "member": figure("running", "member"), "line": figure("running", "line"),
            "moved": text.raw("ui", "view.moved"), "build": set_words(shown),
            "build_contrast": said.label,
            "differs": text.raw("function", "figures.running.differs"),
            "export": dict(exported)}  # fmt: skip


def scatter(said: Words, result: lite.Enriched) -> ScatterData:
    """The running enrichment score, with the set's genes as ticks along the ranking."""
    ranks = tuple(float(i) for i in range(len(result.running)))
    # A rank is a whole number: as a float, the hover printed 12,345 as 1.23e+4.
    frame = pl.DataFrame(
        {
            "rank": pl.Series([int(p) for p in result.positions], dtype=pl.Int64),
            "running": [result.running[p] for p in result.positions],
            "mark": [figure("running", "member")] * result.size,
        }
    )
    columns = {
        "rank": ("rank", "Position in the ranking, from the most positive statistic."),
        "running": ("dimensionless", "Running enrichment score at this member."),
        "mark": ("none", "Set member."),
    }
    member = figure("running", "member")
    return ScatterData(
        frame=frame,
        columns=columns,
        x="rank",
        y="running",
        x_label=figure("running", "x"),
        y_label=figure("running", "y"),
        colour="mark",
        colours={member: ROLES["comparison"]},
        names={member: member},
        caption=caption(said, result, len(ranks)),
        lines=(Line(figure("running", "line"), ranks, result.running),),
        width="double",
    )


def block(writer: Writer, name: str, featured: Contrasted, shown: str,
          others: Sequence[Contrasted] = ()) -> FigureBlock:  # fmt: skip
    """The running score of the dataset's set, and of any set any contrast tested (WO7).

    The page draws another set's running sum through ctxRunning, the twin of lite._running,
    and looks its NES, p and FDR up: the student sees what an enrichment score is for any
    pathway, an early peak, a late trough or a flat wander (figure controls plan). ``featured``
    is the dataset's running contrast; ``others`` each other contrast, in page order.
    """
    said, results, enriched = featured
    hit = enriched[shown]
    data = scatter(said, hit)
    built = writer.figure(data, f"running_{name}",
                          running_text(said, hit, len(hit.running)))  # fmt: skip
    exported = export.settings(data.inches(), data.caption, built.stem, data.export_layers())
    if len(enriched) < 2 and not others:  # noqa: PLR2004 - a choice needs two
        return built
    return replace(built, controls=panel(said, results, enriched, shown, exported,
                                         others))  # fmt: skip
