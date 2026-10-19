"""M6's running score drawn in the page (WO7): ``ctxRunning`` under node against ``lite``.

The twin of ``lite._running`` and of the enrichment score ``lite.prerank`` reports, on planted
rankings: ties in the statistic, a set at the top, one at the bottom, one spread. And the
ranking the page carries is the one ``lite.prerank`` walks, ties broken by name. Each test
names the mutations it catches.
"""

import base64
import dataclasses
import math
from pathlib import Path

import numpy as np
from test_page_script import _node

from ctx.function import lite, running

TEMPLATES = Path(__file__).resolve().parents[1] / "src" / "ctx" / "report" / "templates"
TWIN = TEMPLATES / "_running.js.j2"
HARNESS = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const w = ctxRunning.values(input.weights);
console.log(JSON.stringify(input.sets.map((p) => {
  const run = ctxRunning.running(w, p);
  return { run: run, es: ctxRunning.score(run) };
})));
"""


def _ranked(seed: int, n: int) -> tuple[list[str], np.ndarray]:
    rng = np.random.default_rng(seed)
    names = [f"G{i:05d}" for i in range(n)]
    # Rounded so that ties occur, which the ranking breaks by name.
    return names, np.round(rng.normal(0, 2, n), 1)


def test_the_page_draws_the_running_sum_and_es_lite_reports():
    """Mutations caught: the misses spread over every gene, members included; the weights
    summed over the ranking rather than the set; the maximum taken on a tie, where lite takes
    the minimum; a member's weight not divided by the set's total."""
    names, scores = _ranked(3, 400)
    ranked, weights = lite.ranking(names, scores)
    sets = {"top": ranked[:20], "bottom": ranked[-25:], "spread": ranked[::17],
            "mixed": ranked[:8] + ranked[-8:]}  # fmt: skip
    found = {r.name: r for r in lite.prerank(names, scores, sets, permutations=20)}
    encoded = base64.b64encode(np.asarray(weights, dtype="<f8").tobytes()).decode("ascii")
    order = list(found)
    seen = _node({"weights": encoded, "sets": [list(found[s].positions) for s in order]},
                 HARNESS, TWIN)  # fmt: skip
    for name, page in zip(order, seen, strict=True):
        assert np.allclose(page["run"], found[name].running, rtol=0, atol=1e-12), name
        assert math.isclose(page["es"], found[name].es, rel_tol=0, abs_tol=1e-12), name
    assert seen[order.index("top")]["es"] > 0 > seen[order.index("bottom")]["es"]


def test_a_tie_between_the_peak_and_the_trough_takes_the_trough():
    """Mutations caught: the larger magnitude taken with >= (the peak on a tie)."""
    weights = base64.b64encode(np.ones(4, dtype="<f8").tobytes()).decode("ascii")
    # Members at ranks 0 and 3 of 4: up 0.5, down 0.5, down 0.5, up 0.5: peak 0.5, trough -0.5.
    (seen,) = _node({"weights": weights, "sets": [[0, 3]]}, HARNESS, TWIN)
    assert seen["run"] == [0.5, 0.0, -0.5, 0.0]
    assert seen["es"] == -0.5


def test_the_ranking_breaks_ties_by_name_and_weighs_the_absolute_score():
    """Mutations caught: a tie left to input order; the score's sign kept in the weight."""
    ranked, weights = lite.ranking(["b", "a", "c", "d"], np.array([1.0, 1.0, -3.0, 2.0]))
    assert ranked == ["d", "a", "b", "c"]
    assert weights.tolist() == [2.0, 1.0, 1.0, 3.0]


KEEP = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify(input.cases.map((c) => ctxRunning.keep(c.names, c.current, c.build))));
"""


def test_a_contrast_change_keeps_the_set_where_that_contrast_tested_it():
    """The running score's contrast select (V2): the twin of ``running.keep`` under node.

    Mutations caught: the build's set taken before the one shown; the first set taken where
    the build's set was tested; a missing set taken as index -1."""
    names = ["HALLMARK_A", "HALLMARK_B", "HALLMARK_C"]
    cases = [{"names": names, "current": "HALLMARK_C", "build": "HALLMARK_A"},
             {"names": names, "current": "HALLMARK_X", "build": "HALLMARK_B"},
             {"names": names, "current": "HALLMARK_X", "build": "HALLMARK_Y"},
             {"names": ["HALLMARK_B"], "current": "HALLMARK_B", "build": "HALLMARK_B"},
             {"names": [], "current": "HALLMARK_A", "build": "HALLMARK_A"}]  # fmt: skip
    expected = [running.keep(c["names"], c["current"], c["build"]) for c in cases]
    assert expected == [2, 1, 0, 0, 0]
    assert _node({"cases": cases}, KEEP, TWIN) == expected


def test_each_contrast_lists_its_own_sets_the_running_set_first():
    """Mutations caught: another contrast listed in the build's order; the running set left
    first in a contrast that did not test it."""
    hit = lite.Enriched("HALLMARK_B", 3, 0.1, 1.0, 0.5, 0.5, (0, 1, 2), (0.0,))
    tested = {n: dataclasses.replace(hit, name=n)
              for n in ("HALLMARK_C", "HALLMARK_B", "HALLMARK_A")}  # fmt: skip
    assert running.order(tested, "HALLMARK_B") == ["HALLMARK_B", "HALLMARK_A", "HALLMARK_C"]
    del tested["HALLMARK_B"]
    assert running.order(tested, "HALLMARK_B") == ["HALLMARK_A", "HALLMARK_C"]
