"""The DE engine on planted truth: recovery, the order of steps, shrinkage, and the two lessons.

Fits share module-scoped fixtures; each takes about 3 s at 1,000 genes (P1.5 spike, D-005).
"""

import ast
import math
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from synth import design_counts

from ctx.de import deseq
from ctx.de.deseq import Contrast, ContrastError, fit, pca, select_samples
from ctx.de.filters import expressed, min_sample_size

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


def test_planted_recovered(basic):
    synthetic, result = basic
    truth = planted(synthetic)
    called = result.results.filter(pl.col("call") != "ns")
    false = sum(truth[g] == 0 for g in called["gene_id"])
    # BH bounds the expected FDR, not one draw's: 0.10 allows the draw's binomial spread.
    assert called.height > 50
    assert false / called.height <= 0.10, (false, called.height)
    signs = [np.sign(truth[g]) == np.sign(lfc) for g, lfc in
             called.filter(pl.col("gene_id").is_in([g for g, v in truth.items() if v]))
             .select("gene_id", "log2FoldChange").rows()]  # fmt: skip
    assert all(signs)


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
    split = deseq.Split(42, "sex")
    contrast = Contrast("s", "condition", "control", "control", (7, 6),
                        filters={"condition": ("control",)}, split=split)  # fmt: skip
    first, _ = select_samples(synthetic.samples, contrast, {})
    again, _ = select_samples(synthetic.samples, contrast, {})
    assert first.equals(again)
    by_sex = first.group_by("sex", "group").len().sort("sex", "group")["len"].to_list()
    assert all(abs(a - b) <= 1 for a, b in zip(by_sex[::2], by_sex[1::2], strict=True))


# --- the two lessons --------------------------------------------------------------------------


def test_adjusted_removes_batch_genes():
    """The D5 lesson at tier-1 sizes: case all long; control 4 long, 13 short."""
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

    def share(result) -> float:
        rows = result.results.filter(pl.col("gene_id").is_in(batch_genes))
        return rows.filter(pl.col("call") != "ns").height / len(batch_genes)

    naive_share, adjusted_share = share(naive), share(adjusted)
    print(f"batch genes called: naive {naive_share:.3f}, adjusted {adjusted_share:.3f}")
    assert naive_share >= 0.80
    assert adjusted_share <= 0.05


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

    def true_calls(result) -> int:
        return sum(truth[g] != 0 for g in result.results.filter(pl.col("call") != "ns")["gene_id"])

    print(f"true calls: paired {true_calls(paired)}, unpaired {true_calls(unpaired)}")
    assert true_calls(paired) > true_calls(unpaired)


# --- missing values, threads, PCA, pandas -----------------------------------------------------


def test_na_reasons():
    # Moderate effects and 600 low-count null genes: filtering them lets BH reject more, which
    # is when pydeseq2's independent filtering removes genes (ds.py _independent_filtering).
    synthetic = design_counts(condition_lfc=1.0)
    rng = np.random.default_rng(1)
    names = synthetic.samples["sample"].to_list()
    low = pl.DataFrame(
        {"gene_id": [f"low{i:03d}" for i in range(600)],
         **{c: rng.negative_binomial(10, 10 / 30, 600).astype(np.int32) for c in names}}
    )  # fmt: skip
    first = names[0]
    counts = pl.concat([synthetic.counts, low]).with_columns(
        pl.when(pl.col("gene_id") == "g00500").then(pl.col(first) * 100 + 5000)
        .otherwise(pl.col(first)).alias(first)
    )  # fmt: skip
    result = fit(counts, synthetic.samples, BASIC, {}, threads=1).results
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
