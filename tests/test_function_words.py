"""M6's bars and words follow what is drawn (item 2, D9): the ORA bars only the sets that pass,
the GSEA highlight line, the background's direction, and the calls against their symbols.

A planted ORA case: 3,000 tested genes, G0001 to G1200 called up (the last two sharing one
symbol, as Step-1's CMAHP and RRN3P3 each name two called genes), G1201 to G1500 down. UP_A
and UP_B hold up calls far above chance; ten null sets and every down call hold none.
"""

from datetime import date
from pathlib import Path

import polars as pl
import pytest

from ctx.dataset import load
from ctx.de.page import words as de_words
from ctx.de.run import Writer
from ctx.de.stats import ALPHA
from ctx.function import background
from ctx.function import page as words
from ctx.function import run as m6
from ctx.report import text
from ctx.report.builder import FigureBlock

ROOT = Path(__file__).resolve().parents[1]
SPEC = load(ROOT / "tests" / "fixtures" / "datasets" / "de_synthetic.toml")
FEATURED = SPEC.de.featured
SAID = de_words(next(c for c in SPEC.contrasts if c.name == FEATURED), SPEC.columns)
N, UP, DOWN = 3000, 1200, 300


def _results() -> pl.DataFrame:
    names = [f"G{i:04d}" for i in range(1, N + 1)]
    names[UP - 1] = names[UP - 2]
    return pl.DataFrame({
        "gene_id": [f"g{i:04d}" for i in range(1, N + 1)], "gene_name": names,
        "padj": [0.5] * N, "call": ["up"] * UP + ["down"] * DOWN + ["ns"] * (N - UP - DOWN),
    })  # fmt: skip


def _genes(first: int, n: int, prefix: str = "G") -> tuple[str, ...]:
    return tuple(f"{prefix}{i:04d}" for i in range(first, first + n))


def _sets(untested_a: int = 0, untested_b: int = 0, outside: int = 0) -> dict:
    """UP_A: 30 up calls and 10 not called, plus ``untested_a`` members never tested; UP_B: 25
    and 5, plus ``untested_b``; ten null sets of 30 genes not called; with ``outside``, one set
    of untested genes too large to test, which only the gene-set union counts."""
    sets = {"UP_A": _genes(1, 30) + _genes(1501, 10) + _genes(1, untested_a, "U"),
            "UP_B": _genes(31, 25) + _genes(1511, 5) + _genes(1, untested_b, "V")}  # fmt: skip
    sets |= {f"NULL_{j:02d}": _genes(1516 + 30 * j, 30) for j in range(10)}
    return sets | ({"OUTSIDE": _genes(1, outside, "W")} if outside else {})


def _gsea(passing: int, failing: int) -> pl.DataFrame:
    rows = [(f"GS_{i:02d}", 50, 0.5, 2.0, 0.001, 0.001 * (i + 1)) for i in range(passing)]
    rows += [(f"GS_{passing + i:02d}", 50, 0.1, 0.5, 0.3, 0.2 + 0.1 * i) for i in range(failing)]
    return pl.DataFrame(rows, schema=list(m6.GSEA_COLUMNS), orient="row")


def _enrichment(sets: dict, gsea: pl.DataFrame) -> m6._Enrichment:
    return m6._Enrichment({FEATURED: _results()}, sets, {FEATURED: gsea}, {FEATURED: {}},
                          Path("planted.gmt"))  # fmt: skip


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    """The page's results where the union lowers the plotted sets, and its self-checks."""
    writer = Writer(tmp_path_factory.mktemp("m6") / "06_function", "desynth",
                    date(2026, 10, 7), {})  # fmt: skip
    found, _, checks = m6._page(writer, SPEC, _enrichment(_sets(400, 300), _gsea(16, 4)))
    return found, checks


def test_ora_draws_only_the_sets_that_pass_and_says_how_many(planted):
    """Step-1's down bars sat at FDR 1, 5 of 7 (a length of 1.6e-13): only 2 of 48 down sets
    passed. Every bar drawn passes; a direction with none draws none, and the lead counts it.
    Mutation: the top 7 per direction drawn again."""
    found, _ = planted
    ora, _ = found[1].items
    assert isinstance(ora, FigureBlock)
    assert {(r["gene_set"], r["direction"]) for r in ora.rows} == {("UP_A", "up"), ("UP_B", "up")}
    assert all(r["fdr"] < ALPHA for r in ora.rows)
    assert ora.text.lead == text.fill("function", "figures.ora.lead", up=2, down=0, alpha=ALPHA)


def test_no_ora_figure_where_no_set_passes(tmp_path):
    """No passing set in either direction: a sentence in place of the figure, no empty bars."""
    writer = Writer(tmp_path / "06_function", "desynth", date(2026, 10, 7), {})
    nulls = {k: v for k, v in _sets().items() if k.startswith("NULL")}
    group = m6._featured(
        writer, SPEC.contrasts[0], {SPEC.contrasts[0].name: SAID}, _enrichment(nulls, _gsea(16, 4))
    )
    assert group.items[0] == text.fill("function", "figures.ora.none", alpha=ALPHA)
    assert not group.tabs and isinstance(group.items[1], FigureBlock)


@pytest.mark.parametrize(("passing", "which"), [(16, "all"), (10, "some"), (0, "no")])
def test_the_gsea_highlight_line_names_the_bars_drawn(tmp_path, passing, which):
    """Step-1 drew 15 GSEA bars, all passing, under "the others are shown for context". The
    line now follows the 15 drawn: 16 passing sets, 10, none. Mutation: the line fixed."""
    writer = Writer(tmp_path / "06_function", "desynth", date(2026, 10, 7), {})
    group = m6._featured(
        writer,
        SPEC.contrasts[0],
        {SPEC.contrasts[0].name: SAID},
        _enrichment(_sets(), _gsea(passing, 20)),
    )
    gsea = group.items[1]
    assert isinstance(gsea, FigureBlock) and len(gsea.rows) == m6.TOP_SETS
    assert gsea.text.reading[-1] == text.fill("function", f"figures.gsea.{which}", alpha=ALPHA)


@pytest.mark.parametrize(("tested", "shown"), [(20, 15), (9, 9), (1, 1)])
def test_the_gsea_caption_counts_the_sets_drawn(tested, shown):
    """Package C: the caption said "the 15 sets with the smallest FDR" whatever was drawn. It
    counts the bars: 15 of 20 tested, all 9 of 9, and the one of 1. Mutation: the caption
    filled with TOP_SETS again."""
    data = m6._gsea_bars(SAID, _gsea(tested // 2, tested - tested // 2))
    assert data.frame.height == shown
    want = text.fill("function", "figures.gsea.caption_one" if shown == 1
                     else "figures.gsea.caption", label=SAID.label, permutations=1000, seed=42,
                     shown=shown)  # fmt: skip
    assert data.caption == want
    assert ("the one set tested." if shown == 1 else f"the {shown} sets with") in data.caption


def test_the_background_words_follow_a_union_that_lowers_significance(planted):
    """Against the gene-set union, 6 of Step-1's 7 plotted sets looked less significant, while
    the reading line and the self-check said the union makes overlaps look rarer. Planted here:
    UP_A and UP_B gain 400 and 300 untested members, so both look weaker. Mutations: the
    reading line or the self-check fixed to the more-significant words."""
    found, checks = planted
    (trap,) = found[3].items
    assert isinstance(trap, FigureBlock)
    assert trap.text.lead.startswith(text.fill(
        "function",
        "figures.background.lead_moved",
        union=2999 + 700,
        tested=2999,
        more=0,
        less=2,
        shown=7,
    ))  # fmt: skip
    assert trap.text.reading[1] == text.word("function", "figures.background.less")
    assert checks[0].prompt == text.word("function", "questions.union_less.prompt")


@pytest.mark.parametrize(
    ("planted_sets", "moved", "direction"),
    [(_sets(outside=2000), (2, 0), "more"), (_sets(0, 300, 2000), (1, 1), "mixed")],
)
def test_the_background_words_follow_each_direction(planted_sets, moved, direction):
    """A union that only adds genes outside UP_A and UP_B makes both look more significant;
    one that also gives UP_B 300 untested members moves one each way."""
    results = _results()
    tested = m6._symbols(results)
    wide = background.backgrounds(tested, planted_sets, None)["union"]
    found = background.tested_against(m6._symbols(results, "up"), tested, wide, planted_sets,
                                      "up", "union")  # fmt: skip
    more, less, shown = found.moved("top")
    assert (more, less) == moved and words.union_direction(more, less) == direction
    reading = words.trap_text(len(tested), len(wide), (more, less), shown,
                              found.passing()).reading[1]  # fmt: skip
    assert reading == text.word("function", f"figures.background.{direction}")
    (check, *_) = words.questions(direction)
    assert check.prompt == text.word("function", f"questions.union_{direction}.prompt")
    assert [c.prompt for c in words.questions(None)] == [
        text.word("function", "questions.rank.prompt")
    ]


def test_the_background_caption_counts_calls_and_symbols(planted):
    """Step-1's caption printed "the same 1681 up calls": no separator, and the symbols of
    M5's 1,683 up calls. It names both. Mutation: the symbol count alone, unformatted."""
    found, _ = planted
    (trap,) = found[3].items
    assert "the same 1,200 up calls, 1,199 gene symbols" in trap.caption
    assert "2,999" in trap.caption and "3,699" in trap.caption
