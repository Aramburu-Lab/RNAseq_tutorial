"""M6's running enrichment score for any set the contrast tested (WO7, D-051, D-012).

The page draws a set's running sum with ``ctxRunning`` (``templates/_running.js.j2``), the twin
of ``lite._running`` tested against it under node; NES, p and FDR, which need the permutations,
are looked up, as are the words. The page carries the ranking's weights in rank order once, as
float64, and each set as its members' ranks. Its ES as the page computes it is checked against
the GSEA table's, and the live line says where they differ.
"""

import base64
from collections.abc import Mapping
from typing import Any

import numpy as np
import polars as pl

from ctx.de.page import Words
from ctx.function import lite
from ctx.function.page import fdr_words, figure, running_text, set_words
from ctx.report import marks, text

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


def panel(said: Words, results: pl.DataFrame, enriched: Mapping[str, lite.Enriched],
          shown: str, exported: Mapping[str, Any]) -> dict[str, Any]:  # fmt: skip
    """The running score's control: the ranking and every tested set, the build's first.

    Args:
        said: The contrast in words.
        results: Its M5 results.
        enriched: Its GSEA results by set.
        shown: The set the figure draws as built.
        exported: The figure's export settings; another set's export changes the caption.

    Returns:
        JSON-ready settings for the page script.
    """
    weight = weights(results)
    order = [shown, *sorted((s for s in enriched if s != shown), key=set_words)]
    sets = []
    for name in order:
        hit = enriched[name]
        words = running_text(said, hit, len(weight))
        written = caption(said, hit, len(weight))
        sets.append({"label": set_words(name), "title": words.title,
                     "lead": str(marks.html(words.lead)), "caption": str(marks.html(written)),
                     "export_caption": marks.plain(written), "es": hit.es,
                     "positions": list(hit.positions)})  # fmt: skip
    encoded = base64.b64encode(np.asarray(weight, dtype="<f8").tobytes()).decode("ascii")
    return {"kind": "running", "label": figure("running", "view"), "sets": sets,
            "weights": encoded, "tolerance": ES_TOLERANCE,
            "member": figure("running", "member"), "line": figure("running", "line"),
            "moved": text.raw("ui", "view.moved"), "build": set_words(shown),
            "differs": text.raw("function", "figures.running.differs"),
            "export": dict(exported)}  # fmt: skip
