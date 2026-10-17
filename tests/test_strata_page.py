"""M8's volcano by biotype at the reader's cut-offs (WO10): ``ctxStrata`` under node against
``ctx.ncrna.figures``, on planted rows.

The rows travel packed (float64 for every number the call rule reads), the page unpacks them
and calls them with ctxCalls.call, the twin D-012 holds to ``classify``; the live line is filled
as Python fills it. Each test names the mutations it catches.
"""

import base64
from pathlib import Path

import numpy as np
import polars as pl
from test_page_script import _node

from ctx.ncrna import figures

TEMPLATES = Path(__file__).resolve().parents[1] / "src" / "ctx" / "report" / "templates"
SCRIPTS = (TEMPLATES / "_calls.js.j2", TEMPLATES / "_strata.js.j2")
HARNESS = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
const cfg = input.cfg, rows = ctxStrata.rows(cfg.rows, cfg);
console.log(JSON.stringify({
  rows: rows.length,
  first: rows[0],
  settings: input.settings.map(([basis, alpha, t]) => {
    const split = ctxStrata.split(rows, cfg, t, basis, alpha);
    const out = { summary: ctxStrata.summary(split, cfg, t, basis, alpha),
                  ns: split[cfg.ns].rows.length };
    cfg.levels.forEach((b) => { out[b] = [split[b].up, split[b].down]; });
    return out;
  }),
}));
"""
CHOSEN = ("protein_coding", "lncRNA")


def _frame() -> pl.DataFrame:
    rng = np.random.default_rng(5)
    n = 400
    p = np.concatenate([rng.uniform(0, 1e-4, 60), rng.uniform(0, 0.2, 140), rng.uniform(0, 1, 200)])
    frame = pl.DataFrame({
        "gene_id": [f"G{i}" for i in range(n)], "gene_name": [f"N{i}" for i in range(n)],
        "biotype": [CHOSEN[i % 2] for i in range(n)], "baseMean": np.ones(n),
        # Some fold changes sit exactly on the presets a reader moves to.
        "log2FoldChange": np.where(np.arange(n) % 7 == 0, 1.0, rng.normal(0, 1.5, n)),
        "pvalue": p, "padj": np.minimum(p * 3, 1.0),
    })  # fmt: skip
    genes = frame.select("gene_id", "biotype")
    return figures.by_biotype(frame.with_columns(padj=pl.col("padj")), genes, CHOSEN)


def test_a_thinned_volcano_carries_what_a_cutoff_can_call_and_files_the_rest():
    """Mutations caught: every gene carried, which put TCGA-BRCA's page near its cap; the files
    of a moved view drawing the build's thinned genes, called ones among them, as not called."""
    rng = np.random.default_rng(11)
    n = 6000
    p = rng.uniform(0, 1, n)
    frame = figures.by_biotype(pl.DataFrame({
        "gene_id": [f"G{i}" for i in range(n)], "gene_name": [f"N{i}" for i in range(n)],
        "biotype": [CHOSEN[i % 2] for i in range(n)], "baseMean": np.ones(n),
        "log2FoldChange": rng.normal(0, 1.5, n), "pvalue": p, "padj": np.minimum(p * 3, 1.0),
    }), pl.DataFrame({"gene_id": [f"G{i}" for i in range(n)],
                      "biotype": [CHOSEN[i % 2] for i in range(n)]}), CHOSEN)  # fmt: skip
    data = figures.volcano(frame, CHOSEN, {"tested": "t", "reference": "r"}, "Made up.")
    drawn = pl.Series(data.thinning.page)
    assert not drawn.all()
    callable_ = drawn | (frame["pvalue"] < 0.1)
    assert data.live["rows"]["n"] == int(callable_.sum()) < n
    never = np.frombuffer(base64.b64decode(next(iter(data.live["never"].values()))), "<f4")
    assert never.size == 2 * int((~callable_).sum())
    harness = """
const input = JSON.parse(require("fs").readFileSync(0, "utf8"));
console.log(JSON.stringify([ctxStrata.thinned(input.cfg, true, "built"),
                            ctxStrata.thinned(input.cfg, false, "built")]));
"""
    moved, still = _node({"cfg": data.live}, harness, SCRIPTS)
    assert moved == data.live["never"] and still == "built"


def test_the_page_calls_each_biotype_as_python_does():
    """Mutations caught: the stratum's adjusted p swapped for the pooled one; a gene above the
    loosest cut-off carried as callable; the fold change read as float32; the line's numbers
    grouped differently from Python's."""
    frame = _frame()
    data = figures.volcano(frame, CHOSEN, {"tested": "t", "reference": "r"}, "Made up.")
    cfg = data.live
    settings = [(b, a, t) for b in figures.STRATA_BASES for a in (0.1, 0.05, 0.001)
                for t in (0.0, 1.0, 1.5)]  # fmt: skip
    seen = _node({"cfg": cfg, "settings": settings}, HARNESS, SCRIPTS)
    carried = frame.filter(pl.Series(data.thinning.page) | (pl.col("pvalue") < 0.1))
    assert seen["rows"] == carried.height
    for (basis, alpha, t), page in zip(settings, seen["settings"], strict=True):
        assert page["summary"] == figures.strata_summary(frame, CHOSEN, basis, alpha, t)
        for b in CHOSEN:
            mine = frame.filter(pl.col("biotype") == b).select(
                call=figures.classify(pl.col(basis), pl.col("log2FoldChange"), alpha, t)
            )["call"]
            assert page[b] == [int((mine == "up").sum()), int((mine == "down").sum())]
    assert seen["first"]["log2FoldChange"] == carried["log2FoldChange"][0]
