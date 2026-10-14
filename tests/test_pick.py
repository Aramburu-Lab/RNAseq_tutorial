"""The gene picker's rules (ctxPick, templates/_pick.js.j2) and Open beside's messages
(ctxSelect.post and received, GS.6) under node, each against its Python definition.

Each test names the mutation it catches; every one was planted and seen to fail.
"""

import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest
from test_select import ROWS, SRC, TABLE

from ctx.report import genes, text

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src" / "ctx" / "report" / "templates"
HARNESS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const run = {
  apply: (current, picked, mode) => ctxPick.apply(current, picked, mode),
  toggle: (current, n) => ctxPick.toggle(current, n),
  push: (history, before, after, max) => ctxPick.push(history, before, after, max),
  counts: (list, calls) => ctxPick.counts(list, calls),
  word: (words, key, n) => ctxPick.word(words, key, n),
  view: (rows, cfg, state) => ctxPick.view(rows, cfg, state),
  post: (sel) => ctxSelect.post(sel),
  received: (data) => ctxSelect.received(data),
};
console.log(JSON.stringify(input.calls.map(([name, ...args]) => {
  try { return { ok: run[name](...args) }; } catch (e) {
    return e.key ? { refused: e.key } : { error: String(e.message) };
  }
})));
"""


def _node(calls: list[list[Any]]) -> list[dict[str, Any]]:
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = "".join(re.sub(r"\{% (end)?raw %\}", "", (TEMPLATES / f).read_text())
                     for f in ("_select.js.j2", "_pick.js.j2"))  # fmt: skip
    done = subprocess.run([node, "-e", script + HARNESS], input=json.dumps({"calls": calls}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return [r.get("ok", r) for r in json.loads(done.stdout)]


def test_a_pick_replaces_adds_or_removes_in_the_one_order():
    """The list after a pick is each gene once, ascending, as genes.select orders a selection.

    Mutations caught: replace keeping the old list; remove adding; repeats kept; the list in
    pick order, which ctxSelect.select would sort and a link would then differ from."""
    calls = [
        ["apply", [5, 9], [9, 3, 3], "replace"],
        ["apply", [5, 9], [9, 3, 3], "add"],
        ["apply", [5, 9], [9, 3], "remove"],
        ["toggle", [5, 9], 9],
        ["toggle", [5], 2],
    ]
    assert _node(calls) == [[3, 9], [3, 5, 9], [5], [5], [2, 5]]
    numbers = [abs(n) for n in genes.select(SRC, [r[0] for r in ROWS[:4]][::-1], TABLE).genes]
    assert _node([["apply", [], numbers[::-1], "add"]]) == [numbers]


def test_undo_keeps_fifty_steps_and_skips_a_change_that_changes_nothing():
    """Mutations caught: no cap; the cap dropping the newest step; a no-op taking a step."""
    history = [[i] for i in range(genes.UNDO_STEPS)]
    out = _node([["push", history, [99], [100], genes.UNDO_STEPS],
                 ["push", [[1]], [2], [2], genes.UNDO_STEPS]])  # fmt: skip
    assert len(out[0]) == 50 and out[0][-1] == [99] and out[0][0] == [1]  # noqa: PLR2004
    assert out[1] == [[1]]


def test_counts_follow_the_live_calls_and_words_take_their_number():
    """Mutations caught: a gene with no live call counted up; the one-gene form unused."""
    words = text.raw("ui", "picker")
    out = _node([["counts", [1, 2, 3, 4], {"1": "up", "2": "down", "3": "up"}],
                 ["word", words, "count", 1], ["word", words, "count", 0],
                 ["word", words, "count", 7]])  # fmt: skip
    assert out[0] == {"n": 4, "up": 2, "down": 1, "ns": 1}
    assert out[1:] == [words["count_one"], words["count_none"], words["count"]]


def test_the_table_searches_sorts_and_pages():
    """Rows matching the search, the selected only if asked, by a column then by number, a
    missing value last, 25 a page.

    Mutations caught: ties left in row order; a missing padj sorted first; the last page's
    numbers past the total; the search case-sensitive."""
    rows = [{"number": n, "id": f"ENSG{n:011d}", "symbol": f"G{n % 7}", "padj": (n % 4) / 10
             if n % 9 else None, "lfc": 0.0, "call": "ns"} for n in range(60, 0, -1)]  # fmt: skip
    cfg = {"page_rows": genes.PAGE_ROWS}
    out = _node([["view", rows, cfg, {"list": [], "sort": "padj", "page": 2}],
                 ["view", rows, cfg, {"list": [2, 9], "only": True}],
                 ["view", rows, cfg, {"list": [], "query": "g3", "sort": "padj"}]])  # fmt: skip
    last, only, found = out
    assert (last["first"], last["last"], last["total"]) == (51, 60, 60)
    assert [r["padj"] for r in last["rows"]][-6:] == [None] * 6
    assert [r["number"] for r in only["rows"]] == [2, 9]
    padj = [(r["padj"], r["number"]) for r in found["rows"]]
    assert all(r["symbol"] == "G3" for r in found["rows"]) and padj == sorted(
        padj, key=lambda p: (p[0] is None, p[0] or 0, p[1])
    )


@pytest.mark.parametrize(
    "data",
    [
        {"ctx": "genes", "v": 1, "src": SRC, "n": 2, "g": "-111335,120217"},
        {"ctx": "genes", "v": 1, "src": SRC, "n": 3, "g": "-111335,120217"},
        {"ctx": "genes", "v": 1, "src": SRC + "&g=1", "n": 1, "g": "5"},
        {"ctx": "genes", "v": True, "src": SRC, "n": 1, "g": "5"},
        {"ctx": "other", "v": 1},
        "genes",
    ],
)
def test_a_message_is_read_as_strictly_as_a_link(data):
    """Open beside's message, in both languages: a list read as decode reads a link; a count
    that differs, a field smuggling another (decode's form refuses it), or a flag for a
    version refused; anything else not a list.

    Mutations caught: the count not checked; a boolean read as 1."""
    try:
        python = genes.received(data)
        expected: Any = None if python is None else {"src": python.src, "genes": list(python.genes)}
    except genes.SelectionError as e:
        expected = {"refused": e.key}
    assert _node([["received", data]]) == [expected]


def test_post_and_received_are_inverses_in_both_languages():
    """Mutation: post's g joined by another separator than the link's."""
    sel = genes.select(SRC, [r[0] for r in ROWS], TABLE)
    posted = genes.post(sel)
    assert _node([["post", {"src": sel.src, "genes": list(sel.genes)}]]) == [posted]
    assert genes.received(posted) == sel
