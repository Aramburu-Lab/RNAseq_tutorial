"""The DE engine on planted truth: recovery, the order of steps, shrinkage, and the two lessons.

Fits share module-scoped fixtures; each takes about 3 s at 1,000 genes (P1.5 spike, D-005).
"""

import ast
import base64
import dataclasses
import math
import shutil
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from synth import design_counts, filtered_counts

from ctx.dataset import Contrast, Split
from ctx.de import deseq
from ctx.de.deseq import ContrastError, fit, pca, select_samples
from ctx.de.filters import expressed, min_sample_size
from ctx.de.page import ma_text
from ctx.de.run import _fit_metrics, _ma, volcano
from ctx.de.stats import ALPHA, LFC_MIN, classify
from ctx.report import text
from ctx.viz import scatter

ROOT = Path(__file__).resolve().parents[1]
BASIC = Contrast("basic", "condition", "control", "case", (6, 6))


def planted(synthetic, column: str = "lfc_condition") -> dict[str, float]:
    return dict(synthetic.truth.select("gene_id", column).rows())


@pytest.fixture(scope="module")
def basic():
    synthetic = design_counts()
    return synthetic, fit(synthetic.counts, synthetic.samples, BASIC, {}, threads=1)


# --- the pre-filter ---------------------------------------------------------------------------


def test_filter_hand_worked():
    # Library sizes 1e6, 1e6, 2e6, 4e6: median 1.5e6, so the cutoff is 10 / 1.5 = 6.67 CPM.
    library = np.array([1e6, 1e6, 2e6, 4e6])
    genes = np.array([
        [7, 7, 14, 28],  # 7 CPM everywhere: kept (4 >= 2 samples, total 56)
        [6, 6, 12, 24],  # 6 CPM: under the cutoff everywhere, dropped
        [0, 0, 14, 28],  # 7 CPM in 2 samples: kept, exactly the smallest group
        [0, 0, 0, 14],   # 3.5 CPM in 1 sample: dropped
    ], dtype=float)  # fmt: skip
    counts = np.vstack([genes, library - genes.sum(axis=0)])
    assert expressed(counts, [2, 2])[:4].tolist() == [True, False, True, False]


def test_min_sample_size_eases_above_large_n():
    assert min_sample_size([4, 13]) == 4
    assert min_sample_size([14, 17]) == pytest.approx(10 + (14 - 10) * 0.7)


# --- planted truth ----------------------------------------------------------------------------


def _significant(results: pl.DataFrame) -> pl.DataFrame:
    """Genes at padj < ALPHA at any fold change: the model's own calls, as R3 and R7 read them.

    The fold-change filter in ``call`` removes most false positives by itself, so a test on
    ``call`` alone passes even with the multiple-testing correction dropped.
    """
    return results.filter(classify(pl.col("padj"), pl.col("log2FoldChange"), lfc_min=0.0) != "ns")


def test_planted_recovered(basic):
    synthetic, result = basic
    truth = planted(synthetic)
    # The stored call is the one rule on the adjusted p, never on the raw one.
    recount = result.results.select(classify(pl.col("padj"), pl.col("log2FoldChange")))
    assert result.results["call"].to_list() == recount.to_series().to_list()
    # BH holds the expected FDP to 0.05 x 0.9 only if null p-values are uniform. The Wald test
    # takes the dispersion estimated from 12 samples as known, so they are not: over seeds
    # 0-199 the padj set's FDP averaged 0.083 and reached 0.149. On the earlier flat data at
    # donor SD 0, putting the true dispersion in place brought the mean from 0.073 to 0.048
    # (D-017). 0.20 is the bound one draw supports; dropping the correction gives at least
    # 0.273. The page's call, with |LFC| >= 1, reached 0.084 against its 0.10. The called set
    # goes last: the sign check below reads it.
    for called, bound in ((_significant(result.results), 0.20),
                          (result.results.filter(pl.col("call") != "ns"), 0.10)):  # fmt: skip
        false = sum(truth[g] == 0 for g in called["gene_id"])
        assert called.height > 50
        assert false / called.height <= bound, (false, called.height)
    signs = [np.sign(truth[g]) == np.sign(lfc) for g, lfc in
             called.filter(pl.col("gene_id").is_in([g for g, v in truth.items() if v]))
             .select("gene_id", "log2FoldChange").rows()]  # fmt: skip
    assert all(signs)


def test_no_shrink_changes_nothing_tested(basic):
    # --no-lfc-shrink skips the shrinkage fit: every other column is the same to the bit,
    # so the test, the calls and the enrichment cannot depend on it; the shrunk ones are empty.
    synthetic, result = basic
    plain = fit(synthetic.counts, synthetic.samples, BASIC, {}, threads=1, shrink=False)
    shrunk = ["log2FoldChange_shrunk", "lfcSE_shrunk"]
    assert plain.results.select(shrunk).null_count().row(0) == (plain.results.height,) * 2
    assert plain.results.drop(shrunk).equals(result.results.drop(shrunk))
    assert result.results.select(shrunk).null_count().row(0) == (0, 0)


def test_ma_draws_the_estimates_the_fit_has(basic):
    synthetic, result = basic
    genes = synthetic.counts.select("gene_id", gene_name=pl.col("gene_id"))
    empty = pl.lit(None, pl.Float64)
    plain = dataclasses.replace(
        result, results=result.results.with_columns(log2FoldChange_shrunk=empty, lfcSE_shrunk=empty)
    )
    assert set(_ma(result, genes).frame["estimate"]) == {"raw", "shrunk"}
    assert set(_ma(plain, genes).frame["estimate"]) == {"raw"}
    assert "--no-lfc-shrink" in _ma(plain, genes).caption
    assert "--no-lfc-shrink" in ma_text(BASIC, shrunk=False).lead


# A page thins at this many units here, so a 1,000-gene fit thins as Step-1's 12,796 genes do.
THIN_AT = 100


def _genes(synthetic) -> pl.DataFrame:
    return synthetic.counts.select("gene_id", gene_name=pl.col("gene_id"))


def test_a_thinned_volcano_page_says_what_it_draws(basic, monkeypatch, tmp_path):
    """Item 2, D7: Step-1's volcano caption called its 12,796 tested genes "all drawn", then
    said the page drew 1,842 of the 5,526 no setting can call, from English written in
    scatter.py. The files' caption claims nothing about drawing; the page's ends with ui.toml's
    sentence, numbers from the thinning. Mutations: ", all drawn" restored in volcano's
    caption; the sentence written in scatter.py again, so an edited ui.toml misses the page."""
    synthetic, result = basic
    monkeypatch.setattr(scatter, "THIN_MAX", THIN_AT)
    data = volcano(result, _genes(synthetic))
    other = data.frame.filter(~pl.col("keep")).height
    stride = math.ceil(other / THIN_AT)
    assert stride > 1, "the case must thin"
    sentence = text.fill(
        "ui", "thinned.rows", drawn=math.ceil(other / stride), units=other,
        rest=text.word("ui", "thinned.callable"), nth=scatter._ordinal(stride),
        order="-log10 raw p", total=data.frame.height,
    )  # fmt: skip
    assert data.caption_text() == f"{data.caption} {sentence}"
    assert "all drawn" not in data.caption_text()
    copy = tmp_path / "content"
    shutil.copytree(Path(text.__file__).parent / "content", copy)
    ui = (copy / "ui.toml").read_text()
    (copy / "ui.toml").write_text(ui.replace("This page draws {drawn:,}", "Shown: {drawn:,}"))
    monkeypatch.setenv(text.CONTENT_ENV, str(copy))
    text.load.cache_clear()
    try:
        assert f"{data.caption} Shown: " in data.caption_text()
    finally:
        text.load.cache_clear()


def test_a_thinned_ma_page_draws_both_estimates_of_a_gene_or_neither(basic, monkeypatch):
    """Item 2, D8: thinned row by row, Step-1's MA page drew one of the two estimates of 1,645 of
    its 9,093 genes not called. A gene's raw and shrunk points are now drawn together or not at
    all, the page keeps its budget of points, and the page's points plus those its own PDF and
    SVG add are the CSV's, trace by trace. Mutations: each row its own unit again (the row
    stride); the budget counted in genes, which grew Step-1's M5 page by 818,839 bytes."""
    synthetic, result = basic
    monkeypatch.setattr(scatter, "THIN_MAX", THIN_AT)
    data = _ma(result, _genes(synthetic))
    other = data.frame.filter(~pl.col("keep")).group_by("gene_id").len()
    assert set(other["len"]) == {2}, "every gene not called has both estimates"
    page = data.on_page().filter(~pl.col("keep")).group_by("gene_id").len()
    assert set(page["len"]) == {2}, "a gene drawn is drawn whole"
    stride = math.ceil(2 * other.height / THIN_AT)
    assert stride > 1 and page.height == math.ceil(other.height / stride)
    assert 2 * page.height <= THIN_AT + 2, "at most one unit over the budget"
    sentence = text.fill(
        "ui", "thinned.units", drawn=page.height, units=other.height,
        rest=text.word("de", "figures.ma.thinned"), nth=scatter._ordinal(stride),
        order="log2 fold change", total=data.frame.height,
    )  # fmt: skip
    assert data.caption_text() == f"{data.caption} {sentence}"
    thinned = data.export_layers()["thinned"]
    for index, level in enumerate(data.colours):
        added = np.frombuffer(base64.b64decode(thinned.get(str(index), "")), "<f4")
        shipped = _pairs(data.on_page(), level) + [tuple(p) for p in added.reshape(-1, 2).tolist()]
        assert sorted(shipped) == sorted(_pairs(data.table(), level)), level


def _pairs(rows: pl.DataFrame, level: str) -> list[tuple[float, float]]:
    """An MA trace's points as float32 (x, y), as the export ships them."""
    xy = rows.filter(pl.col("estimate") == level).select("baseMean", "lfc").to_numpy()
    return [tuple(p) for p in xy.astype("<f4").tolist()]


def test_significant_counts_any_fold_change(basic):
    # R3 and R7 read significance at any fold change; up and down also need |LFC| >= 1. One
    # gene made significant at a small fold change tells the two apart on every draw.
    _, result = basic
    small = pl.col("gene_id") == "g00999"
    r = result.results.with_columns(
        padj=pl.when(small).then(0.01).otherwise(pl.col("padj")),
        log2FoldChange=pl.when(small).then(0.5).otherwise(pl.col("log2FoldChange")),
    ).with_columns(call=classify(pl.col("padj"), pl.col("log2FoldChange")))
    metrics = {m: v for m, v, _ in _fit_metrics({"c": dataclasses.replace(result, results=r)})}
    below = r.filter((pl.col("padj") < ALPHA) & (pl.col("log2FoldChange").abs() < LFC_MIN))
    assert below.height >= 1
    assert metrics["calls.c.significant"] == (
        metrics["calls.c.up"] + metrics["calls.c.down"] + below.height
    )


def test_call_uses_adjusted_p():
    # Moderate effects leave genes with raw p < ALPHA and |LFC| >= LFC_MIN that BH does not
    # pass; on them alone, a call made on the raw p differs from the rule.
    synthetic = design_counts(condition_lfc=1.0)
    contrast = Contrast("c", "condition", "control", "case", (6, 6))
    r = fit(synthetic.counts, synthetic.samples, contrast, {}, threads=1).results
    raw_only = r.filter(
        (pl.col("pvalue") < ALPHA) & (pl.col("log2FoldChange").abs() >= LFC_MIN)
        & ~(pl.col("padj") < ALPHA)
    )  # fmt: skip
    assert raw_only.height > 0, "the case must tell the raw from the adjusted p"
    assert set(raw_only["call"]) == {"ns"}
    recount = r.select(classify(pl.col("padj"), pl.col("log2FoldChange"))).to_series()
    assert r["call"].to_list() == recount.to_list()


def test_reference_sign(basic):
    # "case" sorts before "control": a default alphabetical baseline would flip every sign.
    synthetic, result = basic
    assert sorted([BASIC.tested, BASIC.reference])[0] == BASIC.tested
    truth = planted(synthetic)
    rows = result.results.filter(pl.col("gene_id").is_in([g for g, v in truth.items() if v]))
    agree = [
        np.sign(lfc) == np.sign(truth[g])
        for g, lfc in rows.select("gene_id", "log2FoldChange").rows()
    ]
    assert sum(agree) / len(agree) >= 0.95


def test_raw_kept_after_shrink(basic):
    import pandas as pd  # noqa: PLC0415 - the test rebuilds the fit pydeseq2's way
    from pydeseq2.dds import DeseqDataSet  # noqa: PLC0415
    from pydeseq2.default_inference import DefaultInference  # noqa: PLC0415
    from pydeseq2.ds import DeseqStats  # noqa: PLC0415

    synthetic, result = basic
    used = result.samples
    genes = result.results["gene_id"].to_list()
    names = used["sample"].to_list()
    matrix = synthetic.counts.filter(pl.col("gene_id").is_in(genes)).select(names)
    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=pd.DataFrame(matrix.to_numpy().T, index=used["sample"].to_list(), columns=genes),
        metadata=pd.DataFrame({"group": used["group"].to_list()}, index=used["sample"].to_list()),
        design="~group", inference=inference, quiet=True,
    )  # fmt: skip
    dds.deseq2()
    fresh = DeseqStats(
        dds, contrast=["group", "tested", "reference"], inference=inference, quiet=True
    )
    fresh.summary()
    for column in ("log2FoldChange", "padj"):
        expected = fresh.results_df[column].to_numpy()
        assert np.array_equal(
            result.results[column].fill_null(np.nan).to_numpy(), expected, equal_nan=True
        )
    raw, shrunk = result.results["log2FoldChange"], result.results["log2FoldChange_shrunk"]
    assert (raw != shrunk).sum() > 0


# --- order of the steps -----------------------------------------------------------------------


def test_filter_on_subset():
    synthetic = design_counts(
        (("control", "short", 6), ("case", "short", 6), ("control", "long", 6))
    )
    long = synthetic.samples.filter(pl.col("batch") == "long")["sample"].to_list()
    counts = synthetic.counts.with_columns(
        [pl.when(pl.col("gene_id") == "g00999").then(0 if c not in long else 5000)
         .otherwise(pl.col(c)).alias(c) for c in synthetic.counts.columns[1:]]
    )  # fmt: skip
    contrast = Contrast(
        "sub", "condition", "control", "case", (6, 6), filters={"batch": ("short",)}
    )
    result = fit(counts, synthetic.samples, contrast, {}, threads=1)
    assert "g00999" not in result.results["gene_id"].to_list()
    assert set(result.samples["batch"]) == {"short"}


def _with_replicate(synthetic):
    """The first control donor gets a second sample: fu1 in the long batch; its own is fu2."""
    sheet = synthetic.samples
    first = sheet.row(0, named=True)
    extra = {
        **first,
        "run_accession": "XRR000000",
        "sample": "S000",
        "batch": "long",
        "time_point": "fu1",
    }
    sheet = sheet.with_columns(
        time_point=pl.when(pl.col("sample") == first["sample"]).then(pl.lit("fu2"))
        .otherwise(pl.col("time_point"))
    )  # fmt: skip
    return pl.concat([pl.DataFrame([extra], schema=sheet.schema), sheet]), first


def test_dedupe_after_subset():
    synthetic = design_counts()
    sheet, first = _with_replicate(synthetic)
    contrast = Contrast("d", "condition", "control", "case", (6, 6), filters={"batch": ("short",)})
    used, dropped = select_samples(sheet, contrast, {"time_point": ("", "fu1", "fu2")})
    # Ranked on the whole sheet, fu1 would win, then leave with the long batch: the donor lost.
    assert first["sample"] in used["sample"].to_list()
    assert dropped == ()


def test_a_librarys_lanes_are_one_sample_and_a_second_library_is_dropped():
    """Runs of one library share its rank (D-063): a second lane is part of the sample, not a
    later replicate; a second library of the same donor, cell type and condition is.

    Mutation: ranking runs, as before 2026-10-14 (``rank_replicates`` over the subset): the
    first library's second lane is listed in dropped_runs.
    """
    synthetic = design_counts()
    sheet = synthetic.samples
    first = sheet.row(0, named=True)
    lane = {**first, "run_accession": "XRR900001"}
    other = {**first, "run_accession": "XRR900002", "sample": "S999"}
    sheet = pl.concat([sheet, pl.DataFrame([lane, other], schema=sheet.schema)])
    contrast = Contrast("d", "condition", "control", "case", (6, 6), filters={"batch": ("short",)})
    used, dropped = select_samples(sheet, contrast, {})
    assert dropped == ("XRR900002",)
    assert used["sample"].to_list().count(first["sample"]) == 1
    assert "S999" not in used["sample"].to_list()


def test_pair_after_batch():
    synthetic = design_counts(paired=True)
    last = synthetic.samples.filter(pl.col("condition") == "case")["sample"][-1]
    sheet = synthetic.samples.with_columns(
        batch=pl.when(pl.col("sample") == last).then(pl.lit("long")).otherwise(pl.col("batch"))
    )
    lost = sheet.filter(pl.col("sample") == last)["donor"][0]
    contrast = Contrast("p", "condition", "control", "case", (5, 5), design="~donor + group",
                        filters={"batch": ("short",)}, paired=True)  # fmt: skip
    used, dropped = select_samples(sheet, contrast, {})
    assert lost not in used["donor"].to_list()
    assert len(dropped) == 1
    assert used.group_by("donor").len()["len"].to_list() == [2] * 5


def test_rank_deficient_raises():
    synthetic = design_counts()
    contrast = Contrast("r", "condition", "control", "case", (6, 6), design="~donor + group")
    with pytest.raises(ContrastError, match="rank deficient"):
        fit(synthetic.counts, synthetic.samples, contrast, {}, threads=1)


def test_n_mismatch_raises():
    synthetic = design_counts()
    contrast = Contrast("n", "condition", "control", "case", (5, 6))
    with pytest.raises(ContrastError, match=r"n = \(6, 6\), declared \(5, 6\)"):
        fit(synthetic.counts, synthetic.samples, contrast, {}, threads=1)


def test_split_is_seeded_and_stratified():
    synthetic = design_counts((("control", "short", 13), ("case", "short", 1)))
    split = Split(42, "sex")
    contrast = Contrast("s", "condition", "control", "control", (7, 6),
                        filters={"condition": ("control",)}, split=split)  # fmt: skip
    first, _ = select_samples(synthetic.samples, contrast, {})
    again, _ = select_samples(synthetic.samples, contrast, {})
    assert first.equals(again)
    by_sex = first.group_by("sex", "group").len().sort("sex", "group")["len"].to_list()
    assert all(abs(a - b) <= 1 for a, b in zip(by_sex[::2], by_sex[1::2], strict=True))


# --- the two lessons --------------------------------------------------------------------------


def test_adjusted_removes_batch_genes():
    """The D5 lesson at Step-1 sizes: case all long; control 4 long, 13 short."""
    synthetic = design_counts(
        (("control", "short", 13), ("control", "long", 4), ("case", "long", 14)), batch_lfc=2.0
    )
    batch_genes = [g for g, v in planted(synthetic, "lfc_batch").items() if v]
    naive = fit(
        synthetic.counts,
        synthetic.samples,
        Contrast("naive", "condition", "control", "case", (17, 14)),
        {},
        threads=1,
    )
    adjust = Contrast("adj", "condition", "control", "case", (17, 14), design="~batch + group")
    adjusted = fit(synthetic.counts, synthetic.samples, adjust, {}, threads=1)

    def share(rows: pl.DataFrame) -> float:
        return rows.filter(pl.col("gene_id").is_in(batch_genes)).height / len(batch_genes)

    for name, pick in (("significant", _significant),
                       ("called", lambda r: r.filter(pl.col("call") != "ns"))):  # fmt: skip
        naive_share, adjusted_share = share(pick(naive.results)), share(pick(adjusted.results))
        print(f"batch genes {name}: naive {naive_share:.3f}, adjusted {adjusted_share:.3f}")
        assert naive_share >= 0.80, name
        assert adjusted_share <= 0.05, name


def test_paired_beats_unpaired():
    synthetic = design_counts((("control", "short", 9), ("case", "short", 9)), paired=True,
                              donor_sd=1.0, condition_lfc=1.0)  # fmt: skip
    truth = planted(synthetic)
    paired = fit(synthetic.counts, synthetic.samples,
                 Contrast("p", "condition", "control", "case", (9, 9), design="~donor + group",
                          paired=True), {}, threads=1)  # fmt: skip
    unpaired = fit(
        synthetic.counts,
        synthetic.samples,
        Contrast("u", "condition", "control", "case", (9, 9)),
        {},
        threads=1,
    )

    def true_calls(rows: pl.DataFrame) -> int:
        return sum(truth[g] != 0 for g in rows["gene_id"])

    # R7 compares significance at any fold change; the page's calls add the fold-change filter.
    for name, pick in (("significant", _significant),
                       ("called", lambda r: r.filter(pl.col("call") != "ns"))):  # fmt: skip
        a, b = true_calls(pick(paired.results)), true_calls(pick(unpaired.results))
        print(f"true {name}: paired {a}, unpaired {b}")
        assert a > b, name


# --- missing values, threads, PCA, pandas -----------------------------------------------------


def test_na_reasons():
    # Moderate effects and 600 low-count null genes: filtering them lets BH reject more, which
    # is when pydeseq2's independent filtering removes genes (ds.py _independent_filtering).
    synthetic = design_counts(condition_lfc=1.0)
    result = fit(filtered_counts(synthetic), synthetic.samples, BASIC, {}, threads=1).results
    reason = dict(result.select("gene_id", "na_reason").rows())
    assert reason["g00500"] == "cooks_outlier"
    assert set(reason.values()) == {"tested", "cooks_outlier", "independent_filter"}
    no_p = result.filter(pl.col("pvalue").is_null())
    assert no_p["na_reason"].unique().to_list() == ["cooks_outlier"]
    no_padj = result.filter(pl.col("pvalue").is_not_null() & pl.col("padj").is_null())
    assert no_padj["na_reason"].unique().to_list() == ["independent_filter"]


def test_threads_passed(monkeypatch):
    import pydeseq2.default_inference as inference_module  # noqa: PLC0415

    seen_cpus, seen_jobs = [], []
    original_init = inference_module.DefaultInference.__init__
    original_parallel = inference_module.Parallel

    def spy_init(self, *args, **kwargs):
        seen_cpus.append(kwargs.get("n_cpus"))
        original_init(self, *args, **kwargs)

    def spy_parallel(*args, **kwargs):
        seen_jobs.append(kwargs.get("n_jobs"))
        return original_parallel(*args, **kwargs)

    monkeypatch.setattr(inference_module.DefaultInference, "__init__", spy_init)
    monkeypatch.setattr(inference_module, "Parallel", spy_parallel)
    synthetic = design_counts()
    fit(synthetic.counts, synthetic.samples, BASIC, {}, threads=2)
    assert seen_cpus and set(seen_cpus) == {2}, seen_cpus
    assert seen_jobs and set(seen_jobs) == {2}, seen_jobs


def test_pca_top_genes_by_vst_variance():
    synthetic = design_counts()
    library = synthetic.counts.select(pl.exclude("gene_id")).sum()
    depth = pl.DataFrame(
        {"gene_id": ["depth"], **{c: [int(library[c][0] // 2)] for c in library.columns}}
    )
    counts = pl.concat([synthetic.counts, depth.cast(synthetic.counts.schema)])
    raw_variance = counts.select(pl.exclude("gene_id")).to_numpy().var(axis=1)
    assert counts["gene_id"][int(np.argmax(raw_variance))] == "depth"
    result = pca(counts, synthetic.samples["sample"].to_list(), threads=1)
    assert len(result.genes) == deseq.TOP_GENES
    assert "depth" not in result.genes
    assert (
        math.isclose(sum(result.share), sum(result.share))
        and 0 < result.share[1] <= result.share[0]
    )


def test_pandas_only_at_boundary():
    importers = []
    for path in (ROOT / "src" / "ctx").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            names = [a.name for a in node.names] if isinstance(node, ast.Import) else []
            names += [node.module or ""] if isinstance(node, ast.ImportFrom) else []
            if any(n == "pandas" or n.startswith("pandas.") for n in names):
                importers.append(path.relative_to(ROOT).as_posix())
    assert sorted(set(importers)) == ["src/ctx/de/deseq.py"]
