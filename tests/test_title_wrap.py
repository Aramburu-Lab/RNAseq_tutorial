"""A vertical axis title and a flow's labels end a line only where the next word would cross.

An average glyph width (``style.AVERAGE_CHAR_EM``) wrapped both (D14): M3's "Difference d
(log2(count + 1))", set 130.2 pt wide against the 145.8 pt its axis gives a title, broke after
"d", and a flow label broke Step-1's batch level "Reads ≤ 76 nt (min 35 nt)" at "(min". Each
line is measured here as matplotlib sets it in the PDF, never with ``ctx.viz.metrics``, the
measure under test, and a parenthesised span is found by its own pattern, not by ``words``.
"""

import dataclasses
import re

import polars as pl
import pytest
from matplotlib.patches import Rectangle
from synth import demo_linked
from test_caption_fit import _set_width
from test_page_script import EXPORT, _node

from ctx.viz.annotate import fit_title
from ctx.viz.metrics import table, words, wrap_to
from ctx.viz.style import (
    AXIS_TITLE_SHARE,
    FIGURE_WIDTH,
    FLOW,
    SIZES,
    THEMES,
    apply_matplotlib,
    figure_inches,
)

# Every axis title of more than one word on the Step-1 pages built on 2026-10-07: the y titles
# fit_title wraps in the files and the x titles the page wraps.
TITLES = (
    "Expression (log2 reads)", "Genes called", "Mean expression (normalised reads)",
    "Dispersion (dimensionless)", "log2 FC, ndn_hsct_early", "log2 FC, mono_hsct_early",
    "log2 FC, ndn_pdac_naive", "log2 FC, ndn_pdac_adjusted", "log2 fold change",
    "PC1 (VST units)", "PC2 (VST units)", "log2 fold change, gcsf over steady_state",
    "-log10 raw p", "Significance (-log10 p)", "Gene set", "NES (dimensionless)",
    "Significance (-log10 FDR)", "Enrichment score (dimensionless)", "Cell type and condition",
    "Run, by batch and depth", "Reads (million reads)", "Mean of both (log2(count + 1))",
    "Difference d (log2(count + 1))", "Count with -O (log2(count + 1))",
    "Count with -O, split (log2(count + 1))", "Plain count (log2(count + 1))",
    "Count with -M (log2(count + 1))", "Either strand (log2(count + 1))",
    "Called strand (log2(count + 1))", "Reads (%)", "Assigned reads (%)",
    "PC1, 51% of variance (VST units)", "PC2, 19% of variance (VST units)", "Raw (fragments)",
    "Useful (fragments)", "XIST (log2 CPM)", "chrY genes (log2 CPM)", "Library size (reads)",
    "Genes detected", "Globin (fraction)", "Mitochondrial transcripts (fraction)",
    "Mitochondrial rRNA (fraction)", "Distance from its group (VST units)",
    "Agreement with its group (rho)", "Nuclear rRNA (fraction)",
    "Reads compatible with the library type (fraction)", "Reads Salmon assigned (%)",
    "Reads kept by Trim Galore (%)", "Bases trimmed for quality (%)", "Reads with adapter (%)",
    "Processed reads", "Share of the top gene (fraction)", "Position in read (nt)",
    "Mean quality (Phred)", "fastp -q (Phred)", "fastp -l (nt)", "Reads kept (%)",
    "Reads after fastp", "Reads passing fastp (%)", "Duplicated reads (%)", "GC content (%)",
    "Q30 bases (%)", "Read length (nt)",
    # M3's title before 2026-10-07: two spans side by side.
    "Count with -O, split (--fraction) (log2(count + 1))",
)  # fmt: skip
# A parenthesised span, nested at most once, as every title above holds its unit.
SPAN = re.compile(r"\((?:[^()]|\([^()]*\))*\)")
# The files' heights (93.5 mm for a scatter, strip, box, curve or contour, 85 mm for bars),
# and two shorter ones, the shortest giving a unit wider than its line.
HEIGHTS = {f"{kind}_x{scale:g}": figure_inches(FIGURE_WIDTH, kind)[1] * scale
           for kind, scale in (("scatter", 1.0), ("bar", 1.0), ("scatter", 0.6),
                               ("scatter", 0.3))}  # fmt: skip
# Step-1's read-length batches, as the QC page's parallel sets name them.
BATCHES = {"short": "Reads ≤ 76 nt (min 35 nt)", "long": "Reads > 76 nt (max 500 nt)"}


def _check_lines(lines: list[str], room: float, size: float) -> None:
    """Each line fits ``room`` as matplotlib sets it, or is one word wider than any line; the
    next line's first word would cross it; no parenthesised span is broken."""
    text = " ".join(lines)
    assert all(any(span in line for line in lines) for span in SPAN.findall(text)), lines
    for line in lines:
        assert _set_width(line, size) <= room or len(words(line)) == 1, (line, room)
    for line, after in zip(lines, lines[1:], strict=False):
        assert _set_width(f"{line} {words(after)[0]}", size) > room, (line, after, room)


@pytest.mark.parametrize("height", HEIGHTS.values(), ids=HEIGHTS.keys())
def test_a_title_line_ends_only_where_its_next_word_would_cross_the_axis(height):
    """Every real title at every height, in the room fit_title gives it, AXIS_TITLE_SHARE of
    the figure's height. Mutations caught: the average glyph width restored in fit_title,
    which ends lines early; the title split at every space, which broke a unit."""
    apply_matplotlib(THEMES["white"])
    room, wrapped = height * 72 * AXIS_TITLE_SHARE, 0
    for title in TITLES:
        lines = fit_title(title, height).split("\n")
        assert " ".join(lines) == title
        _check_lines(lines, room, SIZES["axis_label"])
        wrapped += len(lines) > 1
    assert wrapped >= 3


def test_m3s_difference_title_takes_the_one_line_it_fits():
    """The audit's case, its numbers stated apart from the code: at a scatter's height the axis
    gives a 10 pt title 145.8 pt, and the title is set 130.2 pt wide. Mutation caught: the
    average glyph width, which broke it after "d", at 26 characters."""
    apply_matplotlib(THEMES["white"])
    _, height = figure_inches(FIGURE_WIDTH, "scatter")
    title = "Difference d (log2(count + 1))"
    assert height * 72 * AXIS_TITLE_SHARE == pytest.approx(145.8, abs=0.1)
    assert _set_width(title, SIZES["axis_label"]) == pytest.approx(130.2, abs=0.1)
    assert fit_title(title, height) == title


EDGES = (
    "", "  a \t b\n", "open (a b", "close) a b", "a)(b c) d", "x (a (b c) d) e",
    f"a{chr(0xA0)}b (c{chr(0xA0)}d)", f"{BATCHES['short']} 111", f"{BATCHES['long']} 99",
)  # fmt: skip


@pytest.mark.parametrize(
    ("text", "held"),
    [
        ("Difference d (log2(count + 1))", ["Difference", "d", "(log2(count + 1))"]),
        (TITLES[-1], ["Count", "with", "-O,", "split", "(--fraction)", "(log2(count + 1))"]),
        (EDGES[7], ["Reads", "≤", "76", "nt", "(min 35 nt)", "111"]),
        (EDGES[1], ["a", "b"]),
        (EDGES[0], []),
        (EDGES[2], ["open", "(a b"]),
        (EDGES[3], ["close)", "a", "b"]),
        (EDGES[4], ["a)(b", "c)", "d"]),
        (EDGES[5], ["x", "(a (b c) d)", "e"]),
    ],
)
def test_words_hold_a_parenthesised_span_as_one(text, held):
    """A span to its closing parenthesis is one word, a nested one included; an unclosed one
    runs to the end, and a closing one with no opening holds nothing. Mutation caught: words
    splitting at every space."""
    assert words(text) == held


TWIN = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.titles.map((text) => {
  const held = ctxExport.words(text);
  return [held, input.rooms.map((room) => ctxExport.wrap(held, (line) =>
    ctxExport.width(line, input.size, input.metrics) <= room))];
})));
"""


def test_the_page_splits_and_wraps_a_title_where_the_files_do():
    """ctxExport.words and ctxExport.wrap under node, with the page's numbers for the face
    (metrics.table), against words and wrap_to on every title above and the edge cases, at the
    rooms of HEIGHTS: the page's titles break where the files' titles break. JavaScript's \\s
    and Python's str.isspace differ only on U+FEFF, U+0085 and U+001C to U+001F (every code
    point compared, 2026-10-07), as ctxExport.wrap and wrap_to already do; no title holds one.
    Mutations caught: words splitting at every space; a span held only while its depth is
    above 1."""
    titles = [*TITLES, *EDGES]
    rooms = [height * 72 * AXIS_TITLE_SHARE for height in HEIGHTS.values()]
    size = SIZES["axis_label"]
    payload = {"titles": titles, "rooms": rooms, "size": size,
               "metrics": table(frozenset("".join(titles)))}  # fmt: skip
    out = _node(payload, TWIN, EXPORT)
    assert out == [[words(t), [wrap_to(words(t), room, size) for room in rooms]] for t in titles]


def _flow_labels(fig) -> list[tuple[list[str], float, float]]:
    """Each node label's lines, its room as drawn, from its own x to a pad before the next
    stage's bars or the axes' right edge, in pt, and its size."""
    (ax,) = fig.axes
    bars = sorted({p.get_x() for p in ax.patches if isinstance(p, Rectangle)})
    axes_pt = ax.get_position().width * fig.get_figwidth() * 72
    labels = []
    for text in ax.texts:
        if text.get_va() != "center":  # a stage's title hangs from the top
            continue
        x = text.get_position()[0]
        end = min((b for b in bars if b > x), default=1.0)
        labels.append((text.get_text().split("\n"), (end - FLOW["pad"] - x) * axes_pt,
                       text.get_fontsize()))  # fmt: skip
    return labels


def test_a_flow_label_line_ends_only_where_its_next_word_would_cross_the_next_stage():
    """Step-1's batch levels in the parallel sets' flow: every line fits before the next
    stage's bars as matplotlib sets it, the next word would cross them, "(min 35 nt)" stays
    whole. Mutations caught: the average glyph width restored in FlowData._labels with
    textwrap, which ended a line at "(min" as Step-1's file of 2026-10-07 does; the same
    estimate over held words, which ended "Reads ≤ 76" a word early."""
    apply_matplotlib(THEMES["white"])
    demo = demo_linked()
    batches = demo.frame.with_columns(pl.col("batch").replace(BATCHES))
    labels = _flow_labels(dataclasses.replace(demo, frame=batches).static(THEMES["white"]))
    for lines, room, size in labels:
        _check_lines(lines, room, size)
    levels = {" ".join(lines).rsplit(" ", 1)[0] for lines, _, _ in labels}
    assert set(BATCHES.values()) <= levels
    assert sum(len(lines) > 1 for lines, _, _ in labels) >= len(BATCHES)
