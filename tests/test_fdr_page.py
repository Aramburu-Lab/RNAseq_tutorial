"""M6's FDR cutoff (P2S.5, D-109): ``ctxFdr`` (templates/_fdr.js.j2) and ``ctxOra`` at each
preset against the Python that builds the figures at that cutoff, under node.

Each figure is built by ``ctx function``'s own contracts at the build's FDR, then redrawn by the
page twin at each preset; the result must be the contract built at that preset, trace for
trace. The fixtures hold FDRs equal to the presets, where < and <= part. Each test names the
mutation it catches; every one was planted and seen to fail.
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

from ctx.de.page import Words
from ctx.de.stats import ALPHA
from ctx.function import cutoff, lite, live, run
from ctx.function.page import gsea_text, nes_text, ora_text
from ctx.viz.scatter import CUTOFF_PRESETS

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ROOT / "src" / "ctx" / "report" / "templates"
SAID = Words("Case against control", "case", "control")
CUTS = [float(c) for c in cutoff.PRESETS]
BOUND = 1.52e-10
MARK = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map(c => input.cuts.map(alpha => {
  const m = ctxFdr.mark(c.figure.data, c.figure.layout, c.spec, alpha);
  return {data: m.data, showlegend: m.showlegend, drawn: m.drawn, passing: m.passing,
          reading: ctxFdr.reading(c.spec, m.drawn, m.passing, alpha),
          lead: c.spec.lead ? c.spec.fdrs.map(f => ctxFdr.fill(c.spec.lead,
                  {passing: ctxFdr.count(f, alpha), tested: f.length, alpha})) : null,
          moved: ctxFdr.fill(c.spec.moved, {alpha, build: c.spec.alpha})};
}))));
"""
# The parts of a bar trace the page draws; plotly's JSON holds no others for these bars.
KEYS = ("type", "name", "orientation", "marker", "hovertemplate", "x", "y", "customdata")


def _node(harness: str, payload: dict, *scripts: str) -> list:
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = "".join(re.sub(r"\{% (end)?raw %\}", "", (TEMPLATES / s).read_text()) for s in scripts)
    done = subprocess.run([node, "-e", script + harness], input=json.dumps(payload),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def _plain(fig) -> dict:
    """A plotly figure as the page holds it: JSON, arrays as lists."""
    from plotly.utils import PlotlyJSONEncoder  # noqa: PLC0415 - heavy

    return json.loads(json.dumps(fig.to_plotly_json(), cls=PlotlyJSONEncoder))


def _gsea(seed: int, n: int = 30) -> pl.DataFrame:
    """A GSEA table with FDRs at each preset, at 0, and either side of each preset."""
    rng = random.Random(seed)
    fdrs = [*CUTS, 0.0, *(c * 0.999 for c in CUTS), *(c * 1.001 for c in CUTS)]
    fdrs += [rng.random() ** 3 for _ in range(n - len(fdrs))]
    rng.shuffle(fdrs)
    return pl.DataFrame({"gene_set": [f"HALLMARK_SET_{i:02d}" for i in range(n)],
                         "size": [rng.randint(15, 200) for _ in range(n)],
                         "es": [rng.uniform(-1, 1) for _ in range(n)],
                         "nes": [rng.uniform(-3, 3) for _ in range(n)],
                         "p": [min(1.0, f + 0.001) for f in fdrs], "fdr": fdrs})  # fmt: skip


def _cases() -> list[dict]:
    """The GSEA bars of three contrasts, as one view, and one set's NES in every contrast."""
    tables = {f"c{i}": _gsea(i) for i in range(3)}
    gsea = run._gsea_bars(SAID, tables["c0"])  # noqa: SLF001 - the contract ctx function draws
    said = {name: Words(f"Contrast {name}", "case", "control") for name in tables}
    # One set whose FDR is each preset in turn, 0 and 1 in a contrast each.
    per = {f"k{i}": pl.DataFrame({"gene_set": ["S"], "nes": [i - 2.5], "fdr": [f]})
           for i, f in enumerate([*CUTS, 0.0, 1.0])}  # fmt: skip
    nes = run._featured_bars("S", per, {n: Words(n, "a", "b") for n in per})  # noqa: SLF001
    return [
        {"figure": _plain(gsea.interactive()), "spec": cutoff.gsea(gsea, list(tables.values())),
         "build": lambda a: run._gsea_bars(SAID, tables["c0"], a),  # noqa: SLF001
         "tables": list(tables.values()), "said": said},
        {"figure": _plain(nes.interactive()), "spec": cutoff.nes(nes),
         "build": lambda a: run._featured_bars("S", per, {n: Words(n, "a", "b") for n in per},
                                               a)},  # noqa: SLF001
    ]  # fmt: skip


def test_the_presets_are_the_volcanos_with_the_build_among_them():
    """Mutation caught: a second list of presets, which the volcano's would no longer move."""
    assert tuple(f"{c:g}" for c in CUTOFF_PRESETS) == cutoff.PRESETS
    assert f"{ALPHA:g}" in cutoff.PRESETS


def test_the_page_marks_the_bars_as_the_build_marks_them_at_each_cutoff():
    """``ctxFdr.mark`` on the build's GSEA and NES bars, at each preset, draws the traces
    ``BarData`` draws at that preset: names, colours, hover, bars and their order, legend.

    Mutations caught: the twin's comparison <= in place of < (an FDR equal to a preset then
    passes); the bars kept in trace order rather than the axis's; an empty series drawn; the
    cutoff left out of the passing name; the legend left on with one trace."""
    cases = _cases()
    seen = _node(MARK, {"cuts": CUTS, "cases": [{k: c[k] for k in ("figure", "spec")}
                                                 for c in cases]}, "_fdr.js.j2")  # fmt: skip
    for case, pages in zip(cases, seen, strict=True):
        for alpha, page in zip(CUTS, pages, strict=True):
            built = _plain(case["build"](alpha).interactive())
            want = [{k: t[k] for k in KEYS} for t in built["data"]]
            assert [{k: t[k] for k in KEYS} for t in page["data"]] == want, alpha
            assert page["showlegend"] == built["layout"]["showlegend"], alpha
    # The fixtures' premises: each preset moves some bar, and an FDR equal to it fails.
    counts = [[p["passing"] for p in pages] for pages in seen]
    assert all(len(set(c)) == len(CUTS) for c in counts), counts


def test_the_pages_words_are_the_builds_at_each_cutoff():
    """The GSEA lead counts each view's tested sets and the reading line names the highlight,
    as ``gsea_text`` words them at that cutoff; NES's line as ``nes_text``'s.

    Mutations caught: the lead counting the bars drawn rather than every set tested; the lead
    left on the build's view; the reading line's all and some swapped."""
    gsea, nes = _cases()
    seen = _node(MARK, {"cuts": CUTS, "cases": [{k: c[k] for k in ("figure", "spec")}
                                                 for c in (gsea, nes)]}, "_fdr.js.j2")  # fmt: skip
    for alpha, page in zip(CUTS, seen[0], strict=True):
        for table, lead in zip(gsea["tables"], page["lead"], strict=True):
            drawn = run._gsea_bars(SAID, table, alpha).frame  # noqa: SLF001
            want = gsea_text(
                SAID,
                int((table["fdr"] < alpha).sum()),
                table.height,
                (int((drawn["fdr"] < alpha).sum()), drawn.height),
                alpha,
            )
            assert lead == want.lead
        built = run._gsea_bars(SAID, gsea["tables"][0], alpha).frame  # noqa: SLF001
        shown = (int((built["fdr"] < alpha).sum()), built.height)
        assert page["reading"] == gsea_text(SAID, 0, 1, shown, alpha).reading[-1]
        assert (
            page["moved"] == f"Sets marked at FDR < {alpha:g}; the figure files mark them "
            f"at FDR < {ALPHA:g}."
        )
    for alpha, page in zip(CUTS, seen[1], strict=True):
        assert page["reading"] == nes_text("S", 6, alpha).reading[-1]
    assert len({p["lead"][1] for p in seen[0]}) == len(CUTS)


ORA = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8")), spec = input.spec;
console.log(JSON.stringify(input.cuts.map(alpha => {
  const families = ["up", "down"].map(d => [d, ctxOra.ora(input[d], input.background,
                                                          input.sets, input.min, input.max)]);
  const blocks = ctxOra.blocks(families, "passing", alpha, input.top);
  const counts = {}, failing = blocks.reduce((n, b) => n + b[1].filter(x => !x.passes).length, 0);
  families.forEach(f => { counts[f[0]] = f[1].filter(r => r.fdr < alpha).length; });
  const values = Object.assign({alpha, failing, top: input.top, background: input.size}, counts);
  return {bars: blocks.flatMap(b => b[1].map(x => [b[0], x.row.name, x.passes])),
          shapes: ctxFdr.line(input.shapes, spec.line, alpha),
          lead: ctxFdr.fill(spec.leads.each, values),
          caption: ctxFdr.fill(spec.caption, values),
          export: ctxFdr.fill(spec.export_caption, values),
          reading: ctxFdr.fill(failing ? spec.reading.failing : spec.reading.passing, values)};
})));
"""


def test_the_live_ora_draws_and_words_each_cutoff_as_the_build_does():
    """The live ORA at each preset: ``ctxOra.blocks``' sets and marks are ``live.bars``' at
    that FDR, the dashed line and its name ``_ora_bars``' line there, the lead, caption and
    reading line ``ora_text``'s and ``_ora_bars``' words at it.

    Mutations caught: the blocks at the build's FDR whatever the cutoff; the line left at
    -log10 0.05; the line named for 0.05; the caption's FDR left at the build's."""
    sets = [[f"S{j:02d}", [f"G{j * 40 + k}" for k in range(40)]] for j in range(40)]
    background = [f"G{i}" for i in range(3_000)]
    # Sets passing at 0.1 but not 0.05, and so on: 1 to 14 hits across the sets, 4 more apart
    # (up FDRs 0.0024, 0.0125, 0.0551, 0.2009 among them).
    up = [f"G{j * 40 + k}" for j in range(14) for k in range(j + 1)]
    up += [f"G{1000 + 40 * i}" for i in range(4)]
    down = [f"G{800 + k}" for k in range(9)] + [f"G{840 + k}" for k in range(4)]
    rows = [(d, r.name, r.size, r.hits, r.expected, r.p, r.fdr, ";".join(r.genes))
            for d, h in (("up", up), ("down", down))
            for r in lite.ora(h, background, dict(sets))]  # fmt: skip
    frame = pl.DataFrame(rows, schema=list(run.ORA_COLUMNS), orient="row")
    built = run._ora_bars(SAID, frame, len(background))  # noqa: SLF001
    reading = len(ora_text(SAID, 1, 0, 0).reading) - 1
    spec = cutoff.ora(run._ora_caption(SAID, len(background), cutoff.HOLE),  # noqa: SLF001
                      len(background), reading)  # fmt: skip
    payload = {"up": up, "down": down, "background": background, "sets": sets,
               "min": lite.MIN_SIZE, "max": lite.MAX_SIZE, "top": live.TOP,
               "size": len(background), "spec": spec, "cuts": CUTS,
               "shapes": _plain(built.interactive())["layout"]["shapes"]}  # fmt: skip
    seen = _node(ORA, payload, "_ora.js.j2", "_fdr.js.j2")
    passing = []
    for alpha, page in zip(CUTS, seen, strict=True):
        bars = live.bars(frame, alpha=alpha)
        assert [tuple(b) for b in page["bars"]] == list(
            zip(bars["direction"], bars["gene_set"], bars["passes"], strict=True)
        ), alpha
        at = run._ora_bars(SAID, frame, len(background), alpha)  # noqa: SLF001
        (line,) = _plain(at.interactive())["layout"]["shapes"]
        (drawn,) = page["shapes"]
        assert drawn["name"] == line["name"] and drawn == {**line, "x0": drawn["x0"],
                                                            "x1": drawn["x1"]}  # fmt: skip
        assert math.isclose(drawn["x0"], -math.log10(alpha), rel_tol=BOUND)
        assert math.isclose(drawn["x1"], line["x1"], rel_tol=BOUND)
        up_n, down_n = (int((frame.filter(pl.col("direction") == d)["fdr"] < alpha).sum())
                        for d in ("up", "down"))  # fmt: skip
        words = ora_text(SAID, len(background), up_n, down_n, int((~bars["passes"]).sum()), alpha)
        assert page["lead"] == words.lead and page["reading"] == words.reading[-1], alpha
        assert page["caption"] == at.caption == page["export"], alpha
        passing.append(int(bars["passes"].sum()))
    # The fixture's premise: each preset passes a different number of sets.
    assert len(set(passing)) == len(CUTS), passing


DRAWN = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cuts.map(alpha => ctxOra.drawn(input.rows, input.how, alpha,
                                                                input.top).map(d => d.passes))));
"""


@pytest.mark.parametrize("how", live.HOW)
def test_drawn_takes_the_cutoff_it_is_given(how):
    """``live.drawn`` and its twin ``ctxOra.drawn`` at each preset mark FDR below it, strictly,
    on sets whose FDR is each preset. Mutations caught: ``<=`` in ``live.drawn``; ``<=`` in
    ``ctxOra.drawn``."""
    rows = [{"name": f"S{i}", "hits": 1, "p": c, "fdr": c} for i, c in enumerate(CUTS)]
    seen = _node(DRAWN, {"rows": rows, "how": how, "top": live.TOP, "cuts": CUTS}, "_ora.js.j2")
    family = pl.DataFrame(rows).rename({"name": "gene_set"})
    for alpha, page in zip(CUTS, seen, strict=True):
        assert page == live.drawn(family, how, alpha)["passes"].to_list(), alpha
    family = pl.DataFrame({"gene_set": [f"S{i}" for i in range(len(CUTS))], "hits": [1] * len(CUTS),
                           "p": CUTS, "fdr": CUTS})  # fmt: skip
    for alpha in CUTS:
        marked = live.drawn(family, how, alpha)
        assert marked.filter(pl.col("fdr") == alpha)["passes"].to_list() == [False]
