"""M9's list from a link (E.1a): the M5 evidence of every gene a contrast tested, for the page.

A list from M5's picker can name any tested gene, and only the build's card list has cards, so
the page carries every tested gene's raw log2 FC, padj and D-012 call in each contrast. The
numbers are rounded to 3 and 2 significant digits, as the gene selection plan measured them
(1,983,413 B over 10 contrasts on Step-1), and the call travels as M5 made it: a padj rounded
at 0.05 cannot recall it. The block is gzipped (mtime 0) and base64-encoded, as M10's panel is.
"""

import base64
import gzip
import json
import math
from collections.abc import Mapping, Sequence
from typing import Any

import polars as pl

from ctx.report import genes
from ctx.viz.heatmap import diverging_levels

# Significant digits each number travels at.
LFC_DIGITS, PADJ_DIGITS = 3, 2
# A call as one character: up, down, tested and not called, not tested.
CALL_CODES = {"up": "u", "down": "d", "ns": "n"}
UNTESTED = "-"
CALLED = ("up", "down")
# The evidence view's steps, most negative first; a step travels as its index.
STEPS = tuple(s.level for s in diverging_levels([])[1])
# Table rows per page of the list's table, as the M5 picker's.
PAGE_ROWS = 25
# Genes the live heatmap draws (E.1d), as the evidence view does (CONVENTIONS: heatmaps cap
# at 100 genes); the rest stay in the table.
HEAT_ROWS = 100


def rounded(value: float | None, digits: int) -> float | None:
    """``value`` at ``digits`` significant digits; None for a missing or non-finite one."""
    if value is None or not math.isfinite(value):
        return None
    return float(f"{value:.{digits}g}")


def evidence(results: Mapping[str, pl.DataFrame], table: genes.GeneTable) -> dict[str, Any]:
    """Every gene a contrast tested, ascending by number, and its evidence in each contrast.

    Args:
        results: Each contrast's M5 table as ``cards.results`` reads it, in page order.
        table: Gather's genes, for each gene's symbol.

    Returns:
        ``numbers`` and ``symbols``, one per gene; ``lfc`` and ``padj``, one list per contrast
        aligned to the genes, None where M5 has no value; ``calls``, one string per contrast
        with a character per gene (``CALL_CODES``, ``UNTESTED`` where padj is missing); and
        ``steps``, the same shape: a called gene's step of the evidence view's scale as a
        digit, from the exact log2 FC, since the rounded one can cross a step's edge, and
        ``UNTESTED`` for the rest.
    """
    tested = sorted({n for frame in results.values()
                     for n in frame.filter(pl.col("padj").is_not_null())["number"]})  # fmt: skip
    at = {n: i for i, n in enumerate(tested)}
    out: dict[str, Any] = {
        "numbers": tested,
        "symbols": [table.symbols[table.row[n]] if n in table.row else "" for n in tested],
        "lfc": [], "padj": [], "calls": [], "steps": [],
    }  # fmt: skip
    for frame in results.values():
        lfc: list[float | None] = [None] * len(tested)
        padj: list[float | None] = [None] * len(tested)
        calls, steps = [UNTESTED] * len(tested), [UNTESTED] * len(tested)
        rows = frame.select("number", "log2FoldChange", "padj", "call").rows()
        levels, _ = diverging_levels([fc for _, fc, _, _ in rows])
        for (n, fc, p, call), level in zip(rows, levels, strict=True):
            if (i := at.get(n)) is None:
                continue
            lfc[i], padj[i] = rounded(fc, LFC_DIGITS), rounded(p, PADJ_DIGITS)
            if p is not None:
                calls[i] = CALL_CODES[call]
            if p is not None and call in CALLED and level is not None:
                steps[i] = str(STEPS.index(level))
        out["lfc"].append(lfc)
        out["padj"].append(padj)
        out["calls"].append("".join(calls))
        out["steps"].append("".join(steps))
    return out


def block(data: Mapping[str, Any]) -> str:
    """The evidence as the page carries it: compact JSON, gzip with no timestamp, base64."""
    body = json.dumps(data, separators=(",", ":"), allow_nan=False).encode()
    return base64.b64encode(gzip.compress(body, mtime=0)).decode("ascii")


def unblock(encoded: str) -> dict[str, Any]:
    """The evidence back from ``block``: what the page inflates, for the tests."""
    data: dict[str, Any] = json.loads(gzip.decompress(base64.b64decode(encoded)))
    return data


def controls(data: Mapping[str, Any], contrasts: Sequence[tuple[str, str]], featured: str,
             cards: Sequence[tuple[int, str, str]], command: str, words: Mapping[str, Any]
             ) -> dict[str, Any]:  # fmt: skip
    """The panel's settings for the page script (``listing`` in ``_page_m9.js.j2``).

    Args:
        data: ``evidence``'s output.
        contrasts: (name, words) of each contrast, in ``data``'s order.
        featured: The contrast the table opens on.
        cards: (number, anchor, symbol) of each gene with a card on the page.
        command: The command that builds cards, to which the page adds ``--genes``.
        words: The panel's words, the picker's for names it cannot read and the selection's
            refusals: ``list``, ``picker`` and ``errors``; ``calls``, a call code's words; and
            ``heat``, the live heatmap's steps, colours and legend (``figures.live_scale``).
    """
    return {"kind": "genelist", "block": block(data), "contrasts": [list(c) for c in contrasts],
            "featured": featured, "cards": [list(c) for c in cards], "command": command,
            "page_rows": PAGE_ROWS, "words": {**words["list"], "calls": words["calls"]},
            "picker": words["picker"], "errors": words["errors"],
            "heat": {**words["heat"], "rows": HEAT_ROWS}}  # fmt: skip
