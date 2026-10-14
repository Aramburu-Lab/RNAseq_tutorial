"""ctx.function.lite against hand-worked values, scipy, and gseapy 1.3.1 as the oracle (D-008)."""

import ast
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom

from ctx.function import lite

ROOT = Path(__file__).resolve().parents[1]
PERMUTATIONS = 1000


def test_es_hand_worked():
    names = ["A", "B", "C", "D", "E", "F"]
    scores = np.array([3.0, 2.0, 1.0, -1.0, -2.0, -3.0])
    (result,) = lite.prerank(names, scores, {"s": ["A", "C"]}, permutations=10, min_size=1)
    # Hits weigh |score| / 4; each of the 4 misses costs 1/4.
    assert result.running == pytest.approx((0.75, 0.5, 0.75, 0.5, 0.25, 0.0))
    assert result.es == pytest.approx(0.75)
    assert result.positions == (0, 2)


@pytest.fixture(scope="module")
def ora_case():
    rng = np.random.default_rng(5)
    background = [f"G{i}" for i in range(2000)]
    hits = list(rng.choice(background, 200, replace=False))
    sets = {f"S{j}": list(rng.choice(background, int(rng.integers(15, 100)), replace=False))
            for j in range(30)}  # fmt: skip
    sets["planted"] = hits[:30] + background[-20:]
    return hits, background, sets


def test_ora_matches_scipy_and_gseapy(ora_case):
    import gseapy  # noqa: PLC0415 - the oracle, tests only

    hits, background, sets = ora_case
    ours = {r.name: r for r in lite.ora(hits, background, sets)}
    for r in ours.values():
        direct = hypergeom.sf(r.hits - 1, len(background), r.size, len(hits))
        assert abs(r.p - direct) <= 1e-12
    oracle = gseapy.enrich(
        gene_list=hits, gene_sets=sets, background=background, outdir=None, cutoff=1.0
    ).results
    compared = 0
    for term, p in oracle[["Term", "P-value"]].itertuples(index=False):
        assert abs(ours[term].p - p) <= 1e-12, term
        compared += 1
    assert compared > 20
    assert ours["planted"].fdr < 0.05


def test_background_is_tested_genes():
    genes = ["a", "b", "c", "d"]
    assert lite.background(genes, [0.1, float("nan"), None, 0.9]) == ["a", "d"]


def test_background_trap():
    """A background larger than the tested genes makes the same hits look rarer."""
    rng = np.random.default_rng(11)
    library = [f"G{i}" for i in range(10000)]
    tested = library[:3000]
    member = library[:60] + library[5000:5040]  # 60 of its 100 genes can be tested
    hits = list(rng.choice(tested[60:], 285, replace=False)) + library[:15]
    sets = {"s": member}
    # gseapy's default with background=None is every gene in the gene-set library (its
    # docstring); here the library's other sets cover all 10,000 genes.
    library_sets = {"s": member, "others": library}
    union = sorted({g for genes in library_sets.values() for g in genes})
    right = lite.ora(hits, tested, sets)[0]
    wrong = lite.ora(hits, union, sets)[0]
    assert (right.size, wrong.size) == (60, 100)
    assert wrong.p < right.p / 10, (wrong.p, right.p)


def test_intersect_then_size_filter():
    universe = {f"G{i}" for i in range(400)}
    big = [f"G{i}" for i in range(600)]  # 600 listed, 400 of them ranked
    assert lite.sized({"big": big}, universe, 15, 500) == {"big": tuple(big[:400])}


def test_collapse_order_free():
    rows = [
        ("g10", "A", 5.0, 100.0), ("g1", "A", 1.0, 100.0), ("g2", "A", 9.0, 50.0),
        ("g3", "B", 2.0, 10.0), ("g4", "", 7.0, 999.0), ("g5", "C", float("nan"), 5.0),
    ]  # fmt: skip
    expected = (["A", "B"], [1.0, 2.0])  # g1 beats g10 on the tie; no symbol and NaN leave
    for seed in range(5):
        order = np.random.default_rng(seed).permutation(len(rows))
        names, scores = lite.collapse(*zip(*[rows[i] for i in order], strict=True))
        assert (names, scores.tolist()) == expected


@pytest.fixture(scope="module")
def ranked():
    rng = np.random.default_rng(7)
    names = [f"G{i}" for i in range(5000)]
    scores = rng.normal(size=5000)
    order = np.argsort(-scores)
    sets = {f"R{j}": list(rng.choice(names, int(rng.integers(15, 200)), replace=False))
            for j in range(40)}  # fmt: skip
    sets["UP"] = [names[i] for i in order[:400:4]]
    sets["DOWN"] = [names[i] for i in order[::-1][:400:4]]
    return names, scores, sets


def test_planted_set_recovered(ranked):
    names, scores, sets = ranked
    results = {r.name: r for r in lite.prerank(names, scores, sets)}
    assert results["UP"].nes > 0 and results["UP"].fdr < 0.05
    assert results["DOWN"].nes < 0 and results["DOWN"].fdr < 0.05
    assert np.median([r.fdr for n, r in results.items() if n.startswith("R")]) >= 0.25


def test_seed_determinism(ranked):
    names, scores, sets = ranked
    first = lite.prerank(names, scores, sets, seed=42)
    assert first == lite.prerank(names, scores, sets, seed=42)
    assert [r.p for r in first] != [r.p for r in lite.prerank(names, scores, sets, seed=43)]


def test_prerank_matches_gseapy(ranked):
    import gseapy  # noqa: PLC0415 - the oracle, tests only

    names, scores, sets = ranked
    ours = {r.name: r for r in lite.prerank(names, scores, sets, permutations=PERMUTATIONS)}
    oracle = gseapy.prerank(
        rnk=pd.Series(scores, index=names), gene_sets=sets, permutation_num=PERMUTATIONS,
        seed=42, min_size=15, max_size=500, outdir=None, threads=1, verbose=False,
    ).res2d  # fmt: skip
    # Rule fixed before the first comparison (P2B.5). Both p-values are Monte-Carlo
    # estimates over about half the permutations (the same-sign ones); their difference has
    # SD sqrt(2) * SE, and a p of 0 has SE 0, hence the floor of one permutation.
    same_sign = PERMUTATIONS / 2
    for term, es, p_oracle, fdr_oracle in oracle[
        ["Term", "ES", "NOM p-val", "FDR q-val"]
    ].itertuples(index=False):
        mine = ours[term]
        assert mine.es == pytest.approx(float(es), abs=1e-9), term
        p = (mine.p + float(p_oracle)) / 2
        tolerance = 3 * math.sqrt(2) * math.sqrt(p * (1 - p) / same_sign) + 1 / (PERMUTATIONS + 1)
        assert abs(mine.p - float(p_oracle)) <= tolerance, (term, mine.p, p_oracle)
        if not 0.02 <= float(fdr_oracle) <= 0.1:
            assert (mine.fdr < 0.05) == (float(fdr_oracle) < 0.05), (term, mine.fdr, fdr_oracle)


def test_nes_matches_gseapy_on_skewed_ranking():
    """NES divides by the same-sign null mean. A long positive tail makes nulls mostly
    positive, so only a set whose ES has the minority sign tells the definitions apart."""
    import gseapy  # noqa: PLC0415 - the oracle, tests only

    rng = np.random.default_rng(9)
    names = [f"G{i}" for i in range(5000)]
    scores = rng.lognormal(0.0, 1.2, 5000) - 1.0
    sets = {f"R{j}": list(rng.choice(names, int(rng.integers(15, 200)), replace=False))
            for j in range(30)}  # fmt: skip
    sets["BOTTOM"] = [names[i] for i in np.argsort(scores)[:15]]
    ours = {r.name: r for r in lite.prerank(names, scores, sets)}
    oracle = gseapy.prerank(
        rnk=pd.Series(scores, index=names), gene_sets=sets, permutation_num=PERMUTATIONS,
        seed=42, min_size=15, max_size=500, outdir=None, threads=1, verbose=False,
    ).res2d  # fmt: skip
    compared = []
    for term, nes in oracle[["Term", "NES"]].itertuples(index=False):
        if math.isfinite(float(nes)):
            assert abs(ours[term].nes - float(nes)) <= 0.10 * abs(float(nes)), term
            compared.append(term)
    assert "BOTTOM" in compared and ours["BOTTOM"].nes < 0


def test_lite_imports_pyodide_safe():
    allowed = {"numpy", "scipy"} | set(sys.stdlib_module_names)
    tree = ast.parse((ROOT / "src" / "ctx" / "function" / "lite.py").read_text())
    imported = {
        a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names
    }
    imported |= {
        (n.module or "").split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)
    }
    assert imported and imported <= allowed, imported - allowed


def test_a_rule_names_a_gene_set_and_m6_records_its_nes_and_fdr():
    """A ``[sanity]`` rule naming ``gsea.<contrast>.<SET>.<nes|fdr>`` gets that set's values;
    an untested set or unknown contrast is NaN with its reason, and other metrics are left out.
    Mutations: ``fdr`` read for ``nes``; an untested set recorded as 0.
    """
    from ctx.function.run import rule_sets

    hit = lite.Enriched("S", 20, 0.5, 1.7, 0.001, 0.02, (0,), (0.5,))
    rules = {"rules": [{"conditions": [
        {"metric": "gsea.c.S.nes", "op": ">", "value": 0},
        {"metric": "x", "op": "<", "other": "gsea.c.S.fdr"},
        {"metric": "gsea.c.T.nes", "op": ">", "value": 0},
        {"metric": "gsea.d.S.fdr", "op": "<", "value": 0.05},
        {"metric": "gene.c.S.lfc", "op": ">", "value": 0},
    ]}]}  # fmt: skip
    found = {m: (v, why) for m, v, why in rule_sets(rules, {"c": {"S": hit}})}
    assert set(found) == {"gsea.c.S.nes", "gsea.c.S.fdr", "gsea.c.T.nes", "gsea.d.S.fdr"}
    assert found["gsea.c.S.nes"] == (1.7, "S in c")
    assert found["gsea.c.S.fdr"][0] == 0.02
    assert math.isnan(found["gsea.c.T.nes"][0])
    assert found["gsea.c.T.nes"][1] == "T in c: not tested"
    assert math.isnan(found["gsea.d.S.fdr"][0])
    assert found["gsea.d.S.fdr"][1] == "no contrast d"
