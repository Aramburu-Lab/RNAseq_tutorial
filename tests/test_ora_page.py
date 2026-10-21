"""M6's live over-representation (GS.5): ctxOra (templates/_ora.js.j2) against lite.ora under
node, the panel's data, and the default state equal to the build's ORA table.

D-049's bound: the prototypes' live ORA matched lite.ora within 1.52e-10 in p on Step-1. Each
test names the mutation it catches; every one was planted and seen to fail.
"""

import json
import math
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


CHOICE = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const cfg = input.cfg, lf = ctxOra.logFactorials(cfg.symbols.length);
const sets = cfg.sets.map(s => [s[0], s[1].map(i => cfg.symbols[i])]);
console.log(JSON.stringify(input.names.map(name => {
  const c = ctxOra.contrast(cfg, name);
  return ["up", "down"].map(d => ctxOra.ora(c[d], c.background, sets, cfg.min_size,
                                            cfg.max_size, lf));
})));
"""


def _contrast(seed: int) -> pl.DataFrame:
    """A contrast's results: symbols shared by several genes, some untested, calls at random."""
    rng = random.Random(seed)
    rows = [(f"ENSG{i:011d}.1", f"S{i % 700}", rng.choice([0.01, 0.5, None]),
             rng.choice(["up", "down", "ns", "ns"]))
            for i in range(1, 900 + 60 * seed)]  # fmt: skip
    frame = pl.DataFrame(rows, schema=["gene_id", "gene_name", "padj", "call"], orient="row")
    return frame.with_columns(call=pl.when(pl.col("padj").is_null()).then(pl.lit("ns"))
                              .otherwise("call"))  # fmt: skip


def test_a_contrast_choice_tests_its_own_calls_as_by_contrast_does():
    """ORA's contrast select (V2): ``ctxOra.contrast`` reads a contrast's own calls and tested
    symbols from the panel; with ``ctxOra.ora`` the rows are ``live.by_contrast``'s.

    Mutations caught: the featured contrast's calls kept for another; another contrast's calls
    tested against the featured contrast's background; up and down swapped."""
    results = {"featured": _contrast(0), "second": _contrast(1), "third": _contrast(2)}
    sets = {f"SET{j}": [f"S{k}" for k in range(j * 11, j * 11 + 45)] for j in range(40)}
    panel = live.panel(results, "featured", sets,
                       {"featured": "F", "second": "S", "third": "T", "union": "U"},
                       {s: s for s in sets})  # fmt: skip
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = re.sub(r"\{% (end)?raw %\}", "", SCRIPT.read_text())
    done = subprocess.run([node, "-e", script + CHOICE],
                          input=json.dumps({"cfg": panel, "names": list(results)}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    seen = json.loads(done.stdout)
    for name, page in zip(results, seen, strict=True):
        table = live.by_contrast(results[name], sets)
        for direction, rows in zip(("up", "down"), page, strict=True):
            built = table.filter(pl.col("direction") == direction)
            assert [(r["name"], r["hits"]) for r in rows] == list(
                zip(built["gene_set"], built["hits"], strict=True)
            ), (name, direction)
            assert [r["p"] for r in rows] == pytest.approx(built["p"].to_list(), rel=BOUND)
    assert seen[1] != seen[2]


BLOCKS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(c => {
  const families = ["up", "down"].map(d => [d, ctxOra.ora(c[d], c.background, c.sets,
                                                          c.min, c.max)]);
  return ctxOra.blocks(families, c.how, c.alpha, c.top).flatMap(b => b[1].map(
    x => [b[0], x.row.name, x.passes, x.row.p, x.length]));
})));
"""


def _short_lists() -> dict[str, dict]:
    """Forty disjoint sets of 40 in 3,000 genes; lists that pass none, one or many sets.

    One gene from each of several sets gives each the same p, so the tie by name is tested;
    "down_first" puts the down family's smallest p first."""
    sets = [[f"S{j:02d}", [f"G{j * 40 + k}" for k in range(40)]] for j in range(40)]
    background = [f"G{i}" for i in range(3_000)]
    one_each = [f"G{j * 40}" for j in range(31, 39)]  # S31 to S38, a hit each
    lists = {
        "none": ([f"G{j * 40}" for j in range(30, 36)], []),
        "one": ([f"G{k}" for k in range(15)] + one_each, [f"G{400 + k}" for k in range(2)]),
        "many": ([f"G{j * 40 + k}" for j in range(10) for k in range(12)],
                 [f"G{800 + k}" for k in range(10)] + [f"G{840 + k}" for k in range(10)]
                 + one_each),
        "down_first": (one_each, [f"G{1200 + k}" for k in range(15)]),
    }  # fmt: skip
    return {name: {"up": up, "down": down, "background": background, "sets": sets,
                   "min": lite.MIN_SIZE, "max": lite.MAX_SIZE}
            for name, (up, down) in lists.items()}  # fmt: skip


def _python_bars(case: dict, how: str) -> pl.DataFrame:
    """``live.bars`` of the case's two families, as ``run._ora`` builds the frame."""
    rows = [(d, r.name, r.size, r.hits, r.expected, r.p, r.fdr, ";".join(r.genes))
            for d in ("up", "down")
            for r in lite.ora(case[d], case["background"], dict(case["sets"]))]  # fmt: skip
    schema = ["direction", "gene_set", "size", "hits", "expected", "p", "fdr", "genes"]
    return live.bars(pl.DataFrame(rows, schema=schema, orient="row"), how)


@pytest.mark.parametrize("how", live.HOW)
def test_the_page_draws_the_bars_live_bars_draws(how):
    """PK.3: ``ctxOra.blocks`` draws ``live.bars``'s sets in its order with its marks, on lists
    passing none, one and many sets; each p within D-049's bound, each length -log10 FDR
    (Rytis 2026-10-21, D-102) within it too.

    Mutations caught: the length -log10 p in ``ctxOra.drawn`` or in ``live.drawn``; the
    failing sets dropped where few pass (Rytis's empty chart); the fill counted past
    ``live.TOP``; a set holding no gene drawn; the tie by name reversed; the families in
    their given order rather than by smallest p."""
    cases = _short_lists()
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = re.sub(r"\{% (end)?raw %\}", "", SCRIPT.read_text())
    plan = [{**c, "how": how, "alpha": live.ALPHA, "top": live.TOP} for c in cases.values()]
    done = subprocess.run([node, "-e", script + BLOCKS], input=json.dumps({"cases": plan}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    seen = dict(zip(cases, json.loads(done.stdout), strict=True))
    counts = {}
    for name, case in cases.items():
        python = _python_bars(case, how)
        assert [(d, s, bool(x)) for d, s, x, *_ in seen[name]] == list(
            zip(python["direction"], python["gene_set"], python["passes"], strict=True)
        ), name
        assert [p for *_, p, _ in seen[name]] == pytest.approx(python["p"].to_list(), rel=BOUND)
        lengths = [-math.log10(f) for f in python["fdr"]]
        assert python[live.LENGTH].to_list() == pytest.approx(lengths, rel=BOUND)
        assert [x for *_, x in seen[name]] == pytest.approx(lengths, rel=BOUND, abs=1e-12)
        assert python["series"].to_list() == [
            d if x else f"{d}_failing"
            for d, x in zip(python["direction"], python["passes"], strict=True)
        ]
        assert (python["hits"] > 0).all()
        per = python.group_by("direction").agg(pl.col("passes").sum(), pl.len())
        counts[name] = {d: (int(x), n) for d, x, n in per.iter_rows()}
    # The fixtures' premises: no set, one set and ten sets pass on the up calls.
    if how == "passing":
        assert counts["none"] == {"up": (0, 6)}
        assert counts["one"]["up"] == (1, live.TOP)
        assert counts["many"] == {"up": (10, 10), "down": (2, live.TOP)}
    else:
        assert all(n <= live.TOP for c in counts.values() for _, n in c.values())
    # The family holding the smallest p first: down on these two, up on "one".
    assert [seen[n][0][0] for n in ("down_first", "many", "one")] == ["down", "down", "up"]
    assert all(seen[n][0][3] == min(p for *_, p, _ in seen[n]) for n in cases)
    # S30 to S35 hold a gene each: an equal p, drawn by name.
    tied = [s for _, s, *_ in seen["none"]]
    assert tied == sorted(tied) and len(tied) == 6  # noqa: PLR2004 - the six genes


def test_drawn_refuses_an_unknown_choice():
    """Mutation caught: an unknown "Sets drawn" read as "passing" with no error."""
    with pytest.raises(ValueError, match="sets drawn"):
        live.drawn(pl.DataFrame({"gene_set": ["S"], "hits": [1], "p": [0.1], "fdr": [0.1]}),
                   "every")  # fmt: skip


NARROW = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(c => ctxOra.narrow(c, input.phone))));
"""
# Labels as M6 draws them: short ones, Hallmark names in words that wrap to two lines, and two
# of the widest, which wrap to four.
LABELS = ("Up set", "Coagulation, MSigDB Hallmark", "Hypoxia, MSigDB Hallmark",
          "Epithelial mesenchymal transition, MSigDB Hallmark",
          "Reactive oxygen species pathway, MSigDB Hallmark", "Down set")  # fmt: skip
# M6's ORA plot at 1,280 px, measured in Chromium on the hand-off fixture (PK.4).
DESKTOP_PLOT_PX = 567


def _ora_bars(labels: tuple[str, ...]):
    """M6's live ORA bars as ``run._featured`` builds them: horizontal, live, every label
    ``drawable``; the build draws the first three."""
    from ctx.viz.bar import BarData  # noqa: PLC0415

    frame = pl.DataFrame({"set": list(labels[:3]), "series": ["up"] * 3, "x": [3.0, 2.0, 1.0]})
    return BarData(frame=frame, columns={"set": ("none", "Set."), "series": ("none", "Family."),
                                         "x": ("none", "-log10 FDR.")},
                   category="set", value="x", category_label="Gene set", value_label="-log10 FDR",
                   series="series", colours={"up": "#E69F00"}, caption="ORA; BH.",
                   horizontal=True, live={"kind": "ora"}, drawable=labels)  # fmt: skip


def test_the_page_builds_a_phones_layout_as_the_build_does():
    """PK.4: ``ctxOra.narrow`` gives the phone's layout ``bar.narrow_categories`` gives, held,
    on the redraw's sets in any order, and none where every label fits; the build's own layout
    keeps plotly's ticks, so a desktop draws it unchanged, and its phone's layout is the one
    the page builds for the same sets.

    Mutations caught: the threshold twice the widest label in either twin (a desktop's 567 px
    plot then wraps); the redraw's ticks left unnamed (plotly skips every other label); the
    build naming its ticks, which a desktop then draws every one of."""
    from ctx.viz.bar import HELD_PAD_PX, narrow_categories  # noqa: PLC0415

    bars = _ora_bars(LABELS)
    phone = bars.controls()["phone"]
    assert set(phone["lines"]) == set(LABELS) and phone["pad"] == HELD_PAD_PX
    rng = random.Random(4)
    cases = [list(LABELS), ["Up set", "Down set"], list(LABELS[:3])]
    cases += [rng.sample(LABELS, rng.randint(1, len(LABELS))) for _ in range(8)]
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = re.sub(r"\{% (end)?raw %\}", "", SCRIPT.read_text())
    done = subprocess.run([node, "-e", script + NARROW],
                          input=json.dumps({"cases": cases, "phone": phone}),
                          capture_output=True, text=True, check=True)  # fmt: skip
    seen = json.loads(done.stdout)
    for cats, page in zip(cases, seen, strict=True):
        assert page == narrow_categories(cats, phone["height"], phone["lines"], held=True), cats
    assert seen[1] is None
    fig = bars.interactive()
    built = fig.layout.meta["narrow"]
    assert fig.layout.yaxis.tickmode == "auto" and fig.layout.yaxis.tickvals is None
    assert built == seen[2]
    widest = max(phone["lines"][c][1] for c in LABELS[:3])
    assert built["below"] == widest + HELD_PAD_PX < DESKTOP_PLOT_PX
