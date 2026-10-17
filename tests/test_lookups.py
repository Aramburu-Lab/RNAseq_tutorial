"""The view kind's lookups (ctx.report.lookups, D-073): light options, their own exports, grids.

The grid's key is one definition in Python and its page twin ``ctxGrid`` runs under node on the
same combinations. Each test names the mutations it catches.
"""

from dataclasses import replace
from pathlib import Path

import pytest
from synth import counts, demo_bar
from test_page_script import _node

from ctx.report.builder import FigureText
from ctx.report.figures import prepare
from ctx.report.lookups import grid_of, key, lookup
from ctx.viz.save import figure_block

GRID = Path(__file__).resolve().parents[1] / "src" / "ctx" / "report" / "templates" / "_grid.js.j2"
INDEX = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const named = input.choices.map((c) => [ctxGrid.key(c), ctxGrid.index(input.cfg, c)]);
console.log(JSON.stringify(named));
"""


def _options(n: int):
    words = FigureText("Bars", reading=("Each bar is made up.",))
    built = demo_bar(counts())
    options = [("a", built, words)]
    for i in range(1, n):
        other = replace(demo_bar(counts(seed=i)), caption=f"Made-up count {i}; no test.")
        options.append((f"o{i}", other, words))
    block = figure_block(built, "2026-10-17_de_lookup_demo_v1", words)
    return block, options


def test_the_page_names_each_combination_as_python_does():
    """Mutations caught: the key joined in another order or separator; a missing key taken as
    another option rather than the build's."""
    cfg = dict(grid_of([("A", ["x", "y"]), ("B", ["p", "q", "r"])],
                       [(1, 2), (0, 0), (0, 1), (0, 2), (1, 0), (1, 1)], 6))  # fmt: skip
    choices = [[a, b] for a in range(2) for b in range(3)] + [[5, 5]]
    seen = _node({"cfg": cfg, "choices": choices}, INDEX, GRID)
    assert seen == [[key(c), cfg["grid"].get(key(c), 0)] for c in choices]
    assert seen[-1] == ["5,5", 0] and seen[5] == ["1,2", 0] and seen[0] == ["0,0", 1]
    assert cfg["axes"][0]["start"] == 1 and cfg["axes"][1]["start"] == 2


def test_a_grid_must_name_every_combination_once():
    """Mutations caught: a repeated combination accepted; a combination left out; a position
    beyond its axis."""
    axes = [("A", ["x", "y"]), ("B", ["p", "q"])]
    with pytest.raises(ValueError, match="names 3 of 4"):
        grid_of(axes, [(0, 0), (0, 1), (1, 0), (0, 0)], 4)
    with pytest.raises(ValueError, match="4 combinations"):
        grid_of(axes, [(0, 0), (0, 1), (1, 0)], 3)
    with pytest.raises(ValueError, match="outside"):
        grid_of(axes, [(0, 0), (0, 1), (1, 0), (2, 1)], 4)


def test_options_drop_the_template_and_carry_their_own_export():
    """Mutations caught: each option carrying the light template again; another option
    exporting the build's caption; the build's export not kept on the controls."""
    block, options = _options(3)
    shown = lookup(block, options, "Which")
    controls = shown.controls
    assert all("template" not in o["figure"]["layout"] for o in controls["options"])
    assert "export" not in controls["options"][0]
    assert controls["export"]["caption"] == options[0][1].caption
    for (_, data, _), o in zip(options[1:], controls["options"][1:], strict=True):
        assert o["export"]["caption"] == data.caption != controls["export"]["caption"]
    assert prepare(shown, 0, 1, 3)["defaults"]["axes"] == []
    with pytest.raises(ValueError, match="as built"):
        lookup(block, options[1:], "Which")
    with pytest.raises(ValueError, match="two or more"):
        lookup(block, options[:1], "Which")
