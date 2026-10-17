"""M5's size factors on disk, the one normalised value, and the one count copy pydeseq2 keeps.

E.1b: M9's box and F1 read the size factors M5 wrote and call ``log2_normalised``, never a
normalisation of their own. S.5: the matrix pydeseq2 fits is built once, in the layout the
earlier three-copy build gave.
"""

import math
from types import SimpleNamespace

import numpy as np
import pandas as pd
import polars as pl
import pytest
from synth import design_counts

from ctx.dataset import Contrast
from ctx.de import deseq
from ctx.de.deseq import fit, samples_by_genes
from ctx.de.normalised import COLUMNS, WHAT, log2_normalised, size_factors


def test_log2_normalised_hand_worked():
    """Two samples whose factors are listed in the other order than the counts' columns.

    Mutations caught: the pseudocount dropped (g1 in s2 gives -inf), a count times its
    factor in place of over it, and factors taken by position rather than joined by sample.
    """
    counts = pl.DataFrame({"gene_id": ["g1", "g2", "g3"], "s1": [10, 3, 50], "s2": [0, 7, 50]})
    factors = pl.DataFrame({"sample": ["s2", "s1"], "size_factor": [0.5, 2.0]})
    got = log2_normalised(counts, factors, ["g2", "g1"])
    assert got.select("gene_id", "sample").rows() == [
        ("g1", "s2"), ("g2", "s2"), ("g1", "s1"), ("g2", "s1"),
    ]  # fmt: skip
    # 0 / 0.5 + 1 = 1; 7 / 0.5 + 1 = 15; 10 / 2 + 1 = 6; 3 / 2 + 1 = 2.5.
    expected = [0.0, math.log2(15), math.log2(6), math.log2(2.5)]
    assert got["value"].to_list() == pytest.approx(expected, abs=1e-12)


def test_the_size_factor_table_joins_each_fits_samples_by_name():
    """The fit lists its samples in one order and its factors in another.

    Mutation caught: the factors laid beside the samples by position (hstack), which gives s1
    the factor of s3.
    """
    one = SimpleNamespace(
        samples=pl.DataFrame({"sample": ["s1", "s2", "s3"], "donor": ["d1", "d2", "d3"],
                              "group": ["reference", "tested", "tested"], "extra": [1, 2, 3]}),
        size_factors=pl.DataFrame({"sample": ["s3", "s1", "s2"], "size_factor": [3.0, 1.0, 2.0]}),
    )  # fmt: skip
    two = SimpleNamespace(
        samples=pl.DataFrame({"sample": ["s2"], "donor": ["d2"], "group": ["reference"]}),
        size_factors=pl.DataFrame({"sample": ["s2"], "size_factor": [0.7]}),
    )
    table = size_factors({"first": one, "second": two})  # type: ignore[dict-item]
    assert table.columns == list(COLUMNS)
    assert table.rows() == [
        ("first", "s1", "d1", "reference", 1.0),
        ("first", "s2", "d2", "tested", 2.0),
        ("first", "s3", "d3", "tested", 3.0),
        ("second", "s2", "d2", "reference", 0.7),
    ]


def test_the_counts_pydeseq2_fits_are_one_int64_copy_in_the_old_layout():
    """Equal to the earlier build, Fortran order, and the very buffer pydeseq2 keeps as X.

    Mutations caught: int32 values, which pydeseq2's ``astype(int)`` copies again (no shared
    buffer); and C order (``to_numpy()`` without ``order="c"``), whose reductions round
    differently from the recorded runs'.
    """
    from pydeseq2.dds import DeseqDataSet  # noqa: PLC0415 - heavy

    rng = np.random.default_rng(3)
    counts = pl.DataFrame(
        {"gene_id": [f"g{i}" for i in range(40)],
         **{f"s{j}": rng.poisson(20, 40).astype(np.int32) for j in range(6)}}
    )  # fmt: skip
    samples = ["s4", "s0", "s2", "s1", "s5", "s3"]
    rows = np.arange(40) % 3 != 0
    frame = samples_by_genes(counts, samples, rows)
    earlier = counts.select(samples).to_numpy()[rows].T
    assert frame.index.to_list() == samples
    assert frame.columns.to_list() == counts["gene_id"].filter(pl.Series(rows)).to_list()
    values = frame.to_numpy()
    assert values.dtype == np.int64 and values.flags.f_contiguous
    assert np.array_equal(values, earlier)
    metadata = pd.DataFrame({"all": ["all"] * len(samples)}, index=samples)
    dds = DeseqDataSet(counts=frame, metadata=metadata, design="~1", quiet=True)
    assert np.shares_memory(dds.X, values)
    assert dds.X.flags.f_contiguous


def test_m5_writes_the_size_factors_its_box_reads(de_run):
    """The written table, read back with its README, gives the box CSV's values.

    One row per sample of every contrast, as the contrasts table counts them; every column in
    the data dictionary. Mutation caught: the box computing its own normalisation (size
    factors refitted over every gene), which moves its values off the written table's.
    """
    root, _ = de_run
    out = root / "outputs" / "05_de"
    (path,) = out.glob(f"*_desynth_{WHAT}_v1.parquet")
    table = pl.read_parquet(path)
    assert table.columns == list(COLUMNS)
    readme = path.with_name(path.name.replace(".parquet", "_README.md")).read_text()
    assert all(
        f"| `{column}` | {unit} | {meaning} |" in readme
        for column, (unit, meaning) in COLUMNS.items()
    )
    contrasts = pl.read_parquet(next(out.glob("*_desynth_contrasts_v1.parquet")))
    expected = dict(contrasts.select("contrast", pl.col("n_reference") + pl.col("n_tested")).rows())
    assert dict(table.group_by("contrast").len().rows()) == expected
    assert table.filter(pl.col("size_factor") <= 0).is_empty()
    counts = pl.read_parquet(next((root / "outputs" / "01_gather").glob("*_counts_v1.parquet")))
    genes = pl.read_parquet(next((root / "outputs" / "01_gather").glob("*_genes_v1.parquet")))
    box = pl.read_csv(next(out.glob("*_desynth_box_a_case_v1.csv")))
    named = dict(genes.select("gene_name", "gene_id").unique("gene_name", keep="first").rows())
    ids = [named[s] for s in box["feature"].unique().to_list()]
    again = log2_normalised(counts, table.filter(pl.col("contrast") == "a_case"), ids)
    value = {(g, s): v for g, s, v in again.rows()}
    drawn = box.select("feature", "sample", "value").rows()
    assert len(drawn) == again.height
    assert all(v == value[(named[f], s)] for f, s, v in drawn)


def test_the_lean_inference_gives_pydeseq2s_numbers(monkeypatch):
    """A paired fit, 9 donors (p = 10), and a blind-VST PCA, each run lean and stock.

    The dispersions run in blocks of 7 genes here, so every block edge is crossed. Every
    column of the fit, the PCA scores, shares and genes equal to the bit. Mutations caught:
    the shrinkage keeping the Hessian's first diagonal entry, not the shrunk coefficient's
    (every shrunk SE moves); and blocks that drop each block's last gene, so the estimates
    no longer line up with the genes.
    """
    from pydeseq2.default_inference import DefaultInference  # noqa: PLC0415 - heavy

    monkeypatch.setattr(deseq, "GENE_BLOCK", 7)
    synthetic = design_counts((("control", "short", 9), ("case", "short", 9)), paired=True,
                              donor_sd=1.0, condition_lfc=1.0)  # fmt: skip
    contrast = Contrast("p", "condition", "control", "case", (9, 9), design="~donor + group",
                        paired=True)  # fmt: skip
    names = synthetic.samples["sample"].to_list()
    ours = fit(synthetic.counts, synthetic.samples, contrast, {}, threads=1)
    ours_pca = deseq.pca(synthetic.counts, names, threads=1)
    monkeypatch.setattr(deseq, "lean_inference", lambda n: DefaultInference(n_cpus=n))
    stock = fit(synthetic.counts, synthetic.samples, contrast, {}, threads=1)
    stock_pca = deseq.pca(synthetic.counts, names, threads=1)
    assert ours.results["lfcSE_shrunk"].null_count() == 0
    assert ours.results.equals(stock.results, null_equal=True)
    assert ours.dispersions.equals(stock.dispersions, null_equal=True)
    assert ours_pca.scores.equals(stock_pca.scores)
    assert (ours_pca.share, ours_pca.genes, ours_pca.trend) == (
        stock_pca.share, stock_pca.genes, stock_pca.trend,
    )  # fmt: skip


def test_the_pca_and_the_fit_run_the_lean_inference(monkeypatch):
    """Both of M5's pydeseq2 entry points take their inference from ``lean_inference``.

    Mutation caught: the PCA's DeseqDataSet built with pydeseq2's stock inference, which
    holds 2.36 GB above entry in its initial dispersions at 1,016 samples (S.5).
    """
    asked = []
    lean = deseq.lean_inference
    monkeypatch.setattr(deseq, "lean_inference", lambda n: asked.append(n) or lean(n))
    synthetic = design_counts()
    deseq.pca(synthetic.counts, synthetic.samples["sample"].to_list(), threads=1)
    assert asked == [1]
    fit(synthetic.counts, synthetic.samples, Contrast("b", "condition", "control", "case", (6, 6)),
        {}, threads=2)  # fmt: skip
    assert asked == [1, 2]
