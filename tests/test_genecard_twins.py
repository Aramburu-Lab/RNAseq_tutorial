"""M9's F3 and F4 page twins (templates/_m9.js.j2) against their Python definitions under node.

F3 recomputes tau and the top cell type over the lineages kept, on either scale; F4 the exon
panel over the tissues kept: each exon's reference, the cells' steps and levels, the CVs and
the qPCR pair. Each test names the mutation it catches; each was planted and seen to fail.
"""

import gzip
import json
import random
import re
import shutil
import subprocess
from pathlib import Path

import polars as pl

from ctx.genecard import annotation, figures, hpa
from ctx.genecard.summary import tau
from ctx.viz.heatmap import diverging_levels

FIXTURES = Path(__file__).parent / "fixtures" / "genecard"
SCRIPT = Path(__file__).resolve().parents[1] / "src" / "ctx" / "report" / "templates" / "_m9.js.j2"
# Python and the page sum in different orders (numpy pairs its sums): a few ulps apart.
CLOSE = 1e-12
HARNESS = r"""
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const out = {
  immune: input.immune.map(c => ctxM9.immune(c.rows, c.lineages, c.log)),
  exons: input.exons.map(c => {
    const u = ctxM9.exonUse(c.rows, c.tissues);
    u.cells.forEach(x => { if (x.ratio === null) delete x.ratio; });
    /* JSON writes NaN as null: name it, so a NaN CV cannot pass for an undefined one. */
    Object.keys(u.cv).forEach(k => { if (Number.isNaN(u.cv[k])) u.cv[k] = "NaN"; });
    return u;
  }),
  levels: (input.levels || []).map(s => ctxM9.level(s)),
  uneven: (input.uneven || []).map(cv => ctxM9.uneven(cv)),
};
console.log(JSON.stringify(out));
"""


def _node(cases: dict) -> dict:
    node = shutil.which("node")
    assert node, "nodejs is pinned in the ctx-dev environment; run tests with pixi"
    script = re.sub(r"\{% (end)?raw %\}", "", SCRIPT.read_text())
    done = subprocess.run([node, "-e", script + HARNESS], input=json.dumps(cases),
                          capture_output=True, text=True, check=True)  # fmt: skip
    return json.loads(done.stdout)


def _close(a: float | None, b: float | None) -> bool:
    if a is None or b is None:
        return a is b
    return abs(a - b) <= CLOSE * max(1.0, abs(a))


def _immune_case(rng: random.Random) -> dict:
    lineages = ["T cells", "B cells", "NK cells", "Monocytes", "Granulocytes", "DCs"]
    rows = [(f"cell{i}", rng.choice(lineages), rng.choice([0.0, 0.3, rng.uniform(0, 400)]))
            for i in range(rng.randint(0, 19))]  # fmt: skip
    if rows and rng.random() < 0.3:
        rows.append((f"tie{len(rows)}", rows[0][1], max(r[2] for r in rows)))  # a tied top
    return {"rows": rows, "lineages": rng.sample(lineages, rng.randint(0, 6)),
            "log": rng.random() < 0.5}  # fmt: skip


def _exon_case(rng: random.Random) -> dict:
    exons = [f"E{k}" for k in range(1, rng.randint(2, 14))]
    tissues = [f"T{k}" for k in range(rng.randint(1, 9))]
    rows = []
    for e in exons:
        silent = rng.random() < 0.2
        for t in tissues:
            share = 0.0 if silent or rng.random() < 0.25 else rng.choice([0.5, rng.random()])
            rows.append((e, t, share))
    return {"rows": rows, "tissues": rng.sample(tissues, rng.randint(0, len(tissues)))}


def test_f3_tau_and_top_equal_python_on_either_scale():
    # Mutations: the log taken on the linear scale; the maximum before the log; a tied top
    # taken last; lineages not filtered; tau defined below the detection limit.
    rng = random.Random(11)
    cases = [_immune_case(rng) for _ in range(300)]
    page = _node({"immune": cases, "exons": []})["immune"]
    for case, got in zip(cases, page, strict=True):
        want = figures.immune_subset(case["rows"], set(case["lineages"]), case["log"])
        assert (got["cells"], got["top"], got["top_ntpm"]) == (
            want["cells"],
            want["top"],
            want["top_ntpm"],
        )
        assert _close(got["tau"], want["tau"]), (case, got, want)
    assert any(c["log"] is False and figures.immune_subset(c["rows"], set(c["lineages"]), False)
               ["tau"] is not None for c in cases)  # fmt: skip


def test_tau_scale_moves_tau():
    # Mutation: the scale ignored. NPTX2's 51 HPA consensus tissues, the real v25 response:
    # 0.6710 on log2(nTPM + 1) and 0.9575 on nTPM, as recorded in D-056 and summary.tau.
    xml = hpa.parse_xml(gzip.decompress((FIXTURES / "hpa_v25_ENSG00000106236.xml.gz")
                                        .read_bytes()), "v25")  # fmt: skip
    values = [v for _, v in xml.tissues]
    assert len(values) == 51  # noqa: PLR2004
    assert (round(tau(values) or 0, 4), round(tau(values, log=False) or 0, 4)) == (0.671, 0.9575)


def test_f4_exon_use_equals_python():
    # Mutations: the reference over every tissue; a CV over one tissue or of an exon without
    # reads; the outer steps lost; the levels' bounds taken left-closed; the qPCR pair over
    # text-ordered exons or taking the last of tied CVs.
    rng = random.Random(7)
    cases = [_exon_case(rng) for _ in range(300)]
    page = _node({"immune": [], "exons": cases})["exons"]
    names = [s.level for s in diverging_levels([])[1]]
    for case, got in zip(cases, page, strict=True):
        shares = pl.DataFrame(case["rows"], schema=["exon", "tissue", "share"], orient="row",
                              schema_overrides={"share": pl.Float64})  # fmt: skip
        use = annotation.exon_use(shares, case["tissues"])
        levels, _ = diverging_levels(list(use.steps))
        want = list(zip(use.frame.iter_rows(named=True), levels, strict=True))
        assert len(got["cells"]) == len(want)
        for cell, (row, level) in zip(got["cells"], want, strict=True):
            assert (cell["exon"], cell["tissue"]) == (row["exon"], row["tissue"])
            assert _close(cell["reference"], row["reference"])
            assert _close(cell.get("ratio"), row["log2_ratio"])
            assert (names[cell["level"]] if cell["level"] is not None else None) == level
        assert set(got["cv"]) == set(use.cv)
        assert all(_close(got["cv"][e], v) for e, v in use.cv.items())
        pair = annotation.qpcr_pair(use.cv)
        if pair is None:
            assert got["pair"] is None
        else:
            assert (got["pair"]["first"], got["pair"]["second"]) == (pair.first, pair.second)
            assert _close(got["pair"]["cv"], pair.cv)
    assert any(p["pair"] is not None for p in page) and any(p["pair"] is None for p in page)


def test_f4_ties_and_bounds_are_planted():
    # Mutations: the first of two equal CVs replaced by the second; a step exactly on a bound
    # put in the lower level; exons ordered as text, E10 before E2.
    rows = [("E1", "A", 0.2), ("E1", "B", 0.4), ("E2", "A", 0.2), ("E2", "B", 0.4),
            ("E10", "A", 0.2), ("E10", "B", 0.4), ("E3", "A", 0.1), ("E3", "B", 0.4)]  # fmt: skip
    cases = {"immune": [], "exons": [{"rows": rows, "tissues": ["A", "B"]}]}
    got = _node(cases)["exons"][0]
    assert (got["pair"]["first"], got["pair"]["second"]) == ("E1", "E2")
    shares = pl.DataFrame(rows, schema=["exon", "tissue", "share"], orient="row")
    pair = annotation.qpcr_pair(annotation.exon_use(shares, ["A", "B"]).cv)
    assert pair is not None and (pair.first, pair.second) == ("E1", "E2")
    level = {(c["exon"], c["tissue"]): c["level"] for c in got["cells"]}
    # E3: reference 0.25, A at 0.1 is log2(0.4) = -1.32; B at 0.4 is log2(1.6) = 0.678.
    assert level[("E3", "A")] == 1 and level[("E3", "B")] == 4  # noqa: PLR2004
    bounds = [0.5, -0.5, 2.0, -2.0, 1.0]
    names = [s.level for s in diverging_levels([])[1]]
    page = _node({"immune": [], "exons": [], "levels": bounds})["levels"]
    assert [names[i] for i in page] == diverging_levels(bounds)[0]
    assert diverging_levels(bounds)[0][:3] == ["0.5 to 1", "-0.5 to 0.5", "≥ 2"]


def test_the_least_even_exon_equals_python_ties_and_undefined_planted():
    # Mutations: the last of tied CVs taken; an undefined CV taken as 0 or as the largest; the
    # exons walked in text order.
    rng = random.Random(5)

    def cv() -> dict:
        exons = rng.sample(range(1, 14), rng.randint(0, 6))
        return {f"E{k}": rng.choice([None, 0.1, 0.4, round(rng.random(), 3)]) for k in exons}

    cases = [cv() for _ in range(300)]
    cases += [{"E2": 0.4, "E10": 0.4, "E1": None}, {"E1": None, "E2": None}, {}]
    page = _node({"immune": [], "exons": [], "uneven": cases})["uneven"]
    for cv, got in zip(cases, page, strict=True):
        want = figures.least_even(cv)
        assert (got is None) is (want is None), (cv, got)
        if want is not None:
            assert (got["exon"], got["cv"]) == want, (cv, got)
    assert figures.least_even({"E2": 0.4, "E10": 0.4, "E1": None}) == ("E2", 0.4)
