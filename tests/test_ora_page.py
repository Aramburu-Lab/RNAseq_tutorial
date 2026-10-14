"""M6's live over-representation (GS.5): ctxOra (templates/_ora.js.j2) against lite.ora under
node, the panel's data, and the default state equal to the build's ORA table.

D-049's bound: the prototypes' live ORA matched lite.ora within 1.52e-10 in p on Step-1. Each
test names the mutation it catches; every one was planted and seen to fail.
"""

import json
import random
import re
import shutil
import subprocess
from pathlib import Path

import polars as pl
import pytest

from ctx.function import lite, live

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "ctx" / "report" / "templates" / "_ora.js.j2"
BOUND = 1.52e-10
HARNESS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(c => ctxOra.ora(c.hits, c.background, c.sets,
                                                           c.min, c.max))));
"""


def _node(cases: list[dict]) -> list[list[dict]]:
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = re.sub(r"\{% (end)?raw %\}", "", SCRIPT.read_text())
    done = subprocess.run([node, "-e", script + HARNESS], input=json.dumps({"cases": cases}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def _case(seed: int) -> dict:
    """A background of 2,000 to 16,000 names, 40 sets of 5 to 600, hits from none to many,
    some outside the background and some repeated."""
    rng = random.Random(seed)
    n = rng.randint(2_000, 16_000)
    background = [f"G{i}" for i in range(n)]
    pool = background + [f"X{i}" for i in range(300)]
    sets = [[f"S{j}", list(dict.fromkeys(rng.sample(pool, rng.randint(5, 600))))]
            for j in range(40)]  # fmt: skip
    hits = rng.sample(pool, rng.choice([0, 3, 50, 400, 1_500]))
    return {"hits": hits + hits[:5], "background": background, "sets": sets,
            "min": lite.MIN_SIZE, "max": lite.MAX_SIZE}  # fmt: skip


@pytest.mark.parametrize("seeds", [range(0, 12)])
def test_the_page_reruns_lite_ora_within_the_bound(seeds):
    """Sets kept, their size, hits and genes equal; p and FDR within D-049's bound.

    Mutations caught: the size cut before the cut to the background; hits not cut to the
    background; sf(k) in place of sf(k - 1); BH without its running minimum; BH's rank off
    by one; the plain running sum of logs (2.6e-10 on Step-1, over the bound)."""
    cases = [_case(s) for s in seeds]
    worst = 0.0
    for case, page in zip(cases, _node(cases), strict=True):
        python = lite.ora(case["hits"], case["background"], dict(case["sets"]))
        assert [(r.name, r.size, r.hits, list(r.genes)) for r in python] == [
            (r["name"], r["size"], r["hits"], r["genes"]) for r in page]  # fmt: skip
        for a, b in zip(python, page, strict=True):
            for x, y in ((a.p, b["p"]), (a.fdr, b["fdr"]), (a.expected, b["expected"])):
                assert abs(x - y) <= BOUND * max(abs(x), 1e-300), (a.name, x, y)
                worst = max(worst, abs(x - y) / max(abs(x), 1e-300))
    assert worst < BOUND


def test_bitsets_read_back():
    """Mutation: the bit order reversed within a byte."""
    members = [0, 1, 7, 8, 9, 4_000, 17_829]
    encoded = live.bitset(members, 17_830)
    assert live.members(encoded, 17_830) == members


def test_the_default_list_is_the_build_ora():
    """With no list from a link, the page tests the featured contrast's up and down calls
    against its tested symbols: the build's ORA table, row for row.

    Mutations caught: the default calls taken from another contrast; the default background
    the union; a call's symbols left with repeats."""
    rows = [
        (
            f"ENSG{i:011d}.1",
            f"S{i % 900}",
            0.01 if i % 3 else 0.5,
            "up" if i % 7 == 0 else "down" if i % 11 == 0 else "ns",
        )
        for i in range(1, 1200)
    ]
    results = pl.DataFrame(rows, schema=["gene_id", "gene_name", "padj", "call"], orient="row")
    other = results.with_columns(pl.lit("ns").alias("call"))
    sets = {f"SET{j}": [f"S{k}" for k in range(j * 13, j * 13 + 40)] for j in range(30)}
    panel = live.panel({"other": other, "featured": results}, "featured", sets,
                       {"other": "Other", "featured": "Featured", "union": "Union"},
                       {s: s for s in sets})  # fmt: skip
    background = dict((b[0], b) for b in panel["backgrounds"])["featured"]
    tested = [panel["symbols"][i] for i in live.members(background[2], len(panel["symbols"]))]
    for call in ("up", "down"):
        hits = [panel["symbols"][i] for i in panel["calls"][call]]
        assert hits == live.symbols(results, call)
        (page,) = _node(
            [
                {
                    "hits": hits,
                    "background": tested,
                    "min": lite.MIN_SIZE,
                    "max": lite.MAX_SIZE,
                    "sets": [[k, v] for k, v in sets.items()],
                }
            ]
        )
        built = lite.ora(live.symbols(results, call), live.symbols(results), sets)
        assert [(r["name"], r["hits"]) for r in page] == [(r.name, r.hits) for r in built]
        assert [r["p"] for r in page] == pytest.approx([r.p for r in built], rel=BOUND)


def test_a_saved_list_is_tested_as_the_page_tests_a_link(tmp_path):
    """``ctx function --genes-from``: a list's genes signed by their log2 FC in the featured
    contrast, each direction its own family against the tested symbols.

    Mutations caught: the sign taken from the list's order; the union as background."""
    rows = [(f"ENSG{i:011d}.1", f"S{i}", 0.01, (-1.0) ** i * i / 100) for i in range(1, 400)]
    results = pl.DataFrame(rows, schema=["gene_id", "gene_name", "padj", "log2FoldChange"],
                           orient="row").with_columns(pl.lit("ns").alias("call"))  # fmt: skip
    # Each set also holds names no gene tested carries, so the union differs from the tested.
    sets = {f"SET{j}": [f"S{k}" for k in range(j * 7 + 1, j * 7 + 30)] + [f"X{j}_{k}" for k in
                                                                          range(40)]
            for j in range(20)}  # fmt: skip
    path = tmp_path / "list.txt"
    path.write_text("\n".join(f"S{i}" for i in range(1, 120)) + "\n")
    table = live.listed(path, results, sets)
    up = [f"S{i}" for i in range(2, 120, 2)]
    down = [f"S{i}" for i in range(1, 120, 2)]
    tested = [f"S{i}" for i in range(1, 400)]
    for direction, hits in (("up", up), ("down", down)):
        expected = lite.ora(hits, tested, sets)
        got = table.filter(pl.col("direction") == direction)
        assert got["gene_set"].to_list() == [r.name for r in expected]
        assert got["p"].to_list() == [r.p for r in expected]
