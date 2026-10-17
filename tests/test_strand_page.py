"""M3's strand call at the reader's thresholds (WO8): ``ctxStrand`` under node against the
registry's ``Metric.strandedness``, on the edges where floats fail. Each test names the
mutations it catches.
"""

from dataclasses import replace
from fractions import Fraction
from pathlib import Path

from test_page_script import _node

from ctx.qc import registry

TEMPLATES = Path(__file__).resolve().parents[1] / "src" / "ctx" / "report" / "templates"
TWIN = TEMPLATES / "_strand.js.j2"
HARNESS = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify({
  calls: input.cases.map(([s1, s2, st, un]) => ctxStrand.call(s1, s2, st, un)),
  clamped: input.typed.map(([v, r, f]) => ctxStrand.threshold(v, r, f)),
}));
"""
# Each pair of totals beside each pair of thresholds: the shares at and either side of each
# edge, 9 of 20 and 11 of 20 among them, mirrors, an empty sample, ten million reads, and a
# share just under 0.95 that a float division rounds onto it, in 2 ** 53 - 1 reads.
TOTALS = [(8, 2), (2, 8), (79, 21), (21, 79), (9, 11), (11, 9), (51, 49), (55, 45), (45, 55),
          (95, 5), (5, 95), (0, 0), (1, 0), (0, 1), (10_000_001, 9_999_999), (7, 13),
          (8_556_839_292_003_941, 450_359_962_737_050)]  # fmt: skip
THRESHOLDS = [("0.8", "0.1"), ("0.75", "0.2"), ("0.95", "0.05"), ("1", "0"), ("0.5", "1"),
              ("0.55", "0.1"), ("0.9", "1e-07")]  # fmt: skip


def _python(s1: int, s2: int, stranded: str, unstranded: str) -> str:
    base = registry.load(registry.resolve(None)).metrics["gene_strand_bias"]
    metric = replace(base, stranded=float(stranded), unstranded=float(unstranded))
    if s1 + s2 == 0:
        return "undetermined"
    return metric.strandedness(Fraction(s1) / (Fraction(s1) + Fraction(s2)))


def test_the_page_calls_every_edge_as_the_registry_does():
    """Mutations caught: forward taken as strictly above the threshold; the unstranded band
    taken as at most its width; the share compared in floats; an exponent read as digits."""
    cases = [[s1, s2, st, un] for s1, s2 in TOTALS for st, un in THRESHOLDS]
    seen = _node({"cases": cases, "typed": []}, HARNESS, TWIN)
    want = [_python(*c) for c in cases]
    assert seen["calls"] == want
    assert {"forward", "reverse", "unstranded", "undetermined"} <= set(want)


def test_a_typed_threshold_is_held_within_its_range():
    """Mutations caught: a box past its range taken as typed; an empty box taken as 0."""
    typed = [["1.3", [0.5, 1], "0.8"], ["0.2", [0.5, 1], "0.8"], ["", [0.5, 1], "0.8"],
             ["0,85", [0.5, 1], "0.8"], ["-1", [0, 1], "0.1"]]  # fmt: skip
    seen = _node({"cases": [], "typed": typed}, HARNESS, TWIN)
    assert seen["clamped"] == ["1", "0.5", "0.8", "0.85", "0"]
