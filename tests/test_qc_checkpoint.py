"""ctx.qc.checkpoint: count metrics against hand-computed values, judged, and tested."""

import dataclasses
import gzip
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from scipy.stats import fisher_exact
from synth import design_counts, write_gather_outputs

from ctx.dataset import Contrast
from ctx.de.stats import adjust
from ctx.qc import QcError, checkpoint, registry

ROOT = Path(__file__).resolve().parents[1]
REG = registry.load(ROOT / "conf" / "qc_metrics.toml")
# gene_id, gene_name, gene_type, chrom; each class once, and a gene of every other kind.
GENES = [
    ("A", "MT-RNR2", "Mt_rRNA", "chrM"),
    ("B", "MT-CO1", "protein_coding", "chrM"),
    ("C", "RNA5-8SN1", "rRNA", "chr21"),
    ("D", "HBB", "protein_coding", "chr11"),
    ("E", "ACTB", "protein_coding", "chr7"),
]


def gtf_line(gene_id: str, name: str, kind: str, chrom: str, feature: str = "gene") -> str:
    attributes = f'gene_id "{gene_id}"; gene_type "{kind}"; gene_name "{name}";'
    return f"{chrom}\tHAVANA\t{feature}\t1\t100\t.\t+\t.\t{attributes}\n"


def write_gtf(path: Path, genes) -> Path:
    lines = ["##description: test\n"]
    for g in genes:
        lines += [gtf_line(*g), gtf_line(*g, feature="exon")]
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as handle:
        handle.writelines(lines)
    return path


def test_gene_types_read_gzip_and_plain_alike(tmp_path):
    plain = checkpoint.gene_types(write_gtf(tmp_path / "a.gtf", GENES))
    zipped = checkpoint.gene_types(write_gtf(tmp_path / "a.gtf.gz", GENES))
    assert plain.equals(zipped)
    assert plain.select("gene_id", "gene_name", "gene_type", "chrom").rows() == GENES


def test_gene_types_coordinates_come_from_the_gene_line(tmp_path):
    """Start, end and strand of the gene record, typed; the exon line after it is ignored.

    Catches: start and end swapped, the strand read from another column, or the columns left
    as text. The exon line, with its own span and strand, stays out.
    """
    gene = 'gene_id "G"; gene_type "lncRNA"; gene_name "LNC";'
    path = tmp_path / "c.gtf"
    path.write_text(f"chr1\tHAVANA\tgene\t1000\t2000\t.\t-\t.\t{gene}\n"
                    f"chr1\tHAVANA\texon\t1500\t1600\t.\t+\t.\t{gene}\n")  # fmt: skip
    found = checkpoint.gene_types(path)
    assert found.row(0, named=True) == {"gene_id": "G", "gene_name": "LNC", "gene_type": "lncRNA",
                                        "chrom": "chr1", "start": 1000, "end": 2000,
                                        "strand": "-"}  # fmt: skip
    assert found.schema["start"] == pl.Int64


def test_count_metrics_by_hand():
    types = pl.DataFrame(GENES, schema=["gene_id", "gene_name", "gene_type", "chrom"], orient="row")
    counts = pl.DataFrame({"gene_id": ["A", "B", "C", "D", "E"],
                           "s1": [10, 20, 30, 40, 900], "s2": [0, 5, 0, 0, 95]})  # fmt: skip
    got = {r["sample"]: r for r in checkpoint.count_metrics(counts, types).iter_rows(named=True)}
    assert got["s1"] == {"sample": "s1", "library_size": 1000.0, "genes_detected": 5,
                         "top_gene_share": 0.9, "mt_fraction": 0.02, "mt_rrna_fraction": 0.01,
                         "rrna_fraction": 0.03, "globin_fraction": 0.04}  # fmt: skip
    assert got["s2"]["genes_detected"] == 2
    assert got["s2"]["mt_fraction"] == 0.05
    assert got["s2"]["top_gene_share"] == 0.95


def test_an_unknown_gene_stops_the_checkpoint():
    types = pl.DataFrame(GENES[:3], schema=["gene_id", "gene_name", "gene_type", "chrom"],
                         orient="row")  # fmt: skip
    counts = pl.DataFrame({"gene_id": ["A", "D", "E"], "s1": [1, 2, 3]})
    with pytest.raises(QcError, match=r"2 gene IDs .* \['D', 'E'\]"):
        checkpoint.count_metrics(counts, types)


def test_replicate_spearman_finds_the_shuffled_sample():
    rng = np.random.default_rng(0)
    base = rng.lognormal(5, 2, 300)
    columns = {f"s{i}": rng.poisson(base) for i in range(5)}
    columns["s4"] = rng.permutation(columns["s4"])
    counts = pl.DataFrame({"gene_id": [f"g{i}" for i in range(300)], **columns, "lone": base})
    samples = pl.DataFrame({"sample": [*columns, "lone"], "group": ["a"] * 5 + ["b"]})
    rho = dict(checkpoint.replicate_spearman(counts, samples, ["group"]).rows())
    assert min(rho, key=lambda s: rho[s] if rho[s] is not None else 2) == "s4"
    assert rho["lone"] is None
    assert all(rho[f"s{i}"] > 0.9 for i in range(4))


def metric(key: str, **changes) -> registry.Metric:
    return dataclasses.replace(REG.metrics[key], **changes)


def reg_of(*metrics) -> registry.Registry:
    return dataclasses.replace(REG, metrics={m.id: m for m in metrics})


def test_judge_flags_the_bad_side_only():
    values = pl.DataFrame({
        "sample": [f"s{i}" for i in range(8)],
        "genes_detected": [100.0, 101, 99, 100, 102, 98, 10, 500],
    })  # fmt: skip
    samples = pl.DataFrame({"sample": values["sample"], "cell": ["x"] * 8})
    bands, flags = checkpoint.judge(values, samples, reg_of(metric("genes_detected")), ["cell"])
    assert flags["sample"].to_list() == ["s6"]
    assert "below" in flags["reason"][0] and "8 samples of x" in flags["reason"][0]
    both = reg_of(metric("genes_detected", direction="none"))
    assert checkpoint.judge(values, samples, both, ["cell"])[1]["sample"].to_list() == ["s6", "s7"]
    assert bands.row(0, named=True)["median"] == 100.0


def test_judge_uses_the_log_scale_and_keeps_keep_metrics_quiet():
    sizes = [1e6, 2e6, 1.5e6, 1.2e6, 1.8e6, 3e7]
    values = pl.DataFrame({"sample": [f"s{i}" for i in range(6)], "library_size": sizes})
    samples = pl.DataFrame({"sample": values["sample"], "cell": ["x"] * 6})
    log_scale = reg_of(metric("library_size", direction="none"))
    bands, flags = checkpoint.judge(values, samples, log_scale, ["cell"])
    assert flags["sample"].to_list() == ["s5"]
    assert bands["median"][0] == pytest.approx(10 ** np.median(np.log10(sizes)))
    kept = reg_of(metric("library_size", direction="none", action="keep"))
    assert checkpoint.judge(values, samples, kept, ["cell"])[1].is_empty()


def test_small_strata_and_zero_spread_get_no_band():
    values = pl.DataFrame({"sample": ["a", "b", "c", "d", "e", "f"],
                           "globin_fraction": [0.0, 0.9, 0.0, 0.0, 0.0, 0.5]})  # fmt: skip
    samples = pl.DataFrame({"sample": values["sample"], "cell": ["x", "x", "y", "y", "y", "y"]})
    bands, flags = checkpoint.judge(values, samples, reg_of(metric("globin_fraction")), ["cell"])
    notes = dict(bands.select("stratum", "note").rows())
    assert notes["x"].startswith("no band: 2 samples, fewer than the 3")
    assert notes["y"].startswith("no band: the MAD is 0")
    assert flags.is_empty()


def test_categories_and_strandedness_are_judged_without_bands():
    values = pl.DataFrame({
        "sample": ["a", "b", "c", "d"],
        "sex_agrees": ["agrees", "disagrees", "not recorded", "agrees"],
        "salmon-strand_mapping_bias": [0.02, 0.03, 0.5, 0.6],
    })  # fmt: skip
    samples = pl.DataFrame({"sample": values["sample"], "cell": ["x"] * 4})
    reg = reg_of(metric("sex_check"), metric("salmon_strand_bias"))
    _, flags = checkpoint.judge(values, samples, reg, ["cell"])
    got = {(s, m) for s, m, _ in flags.rows()}
    # a and b are reverse; c is unstranded against a reverse stratum; d is undetermined.
    assert got == {("b", "sex_check"), ("c", "salmon_strand_bias"), ("d", "salmon_strand_bias")}


def test_status_takes_each_runs_sample_flags():
    runs = pl.DataFrame({"run_accession": ["r1", "r2", "r3"], "sample": ["a", "a", "b"]})
    # Two reasons for one metric, its absolute band and its MAD band, count as one flag.
    flags = pl.DataFrame({"sample": ["a", "a", "a"],
                          "metric": ["mt_fraction", "library_size", "library_size"],
                          "reason": ["one", "two", "three"]})  # fmt: skip
    table = checkpoint.status(runs, flags, REG)
    both = "one; two; three"
    assert table.rows() == [("r1", "a", "flag", both, 2), ("r2", "a", "flag", both, 2),
                            ("r3", "b", "pass", "", 0)]  # fmt: skip
    excluding = reg_of(metric("mt_fraction", action="exclude"), metric("library_size"))
    assert checkpoint.status(runs, flags, excluding)["qc_status"].to_list()[0] == "exclude"


def test_fisher_odds_ratio_is_haldane_corrected_only_at_a_zero():
    p, odds = checkpoint.fisher((3, 2), (1, 4))
    assert p == fisher_exact([[3, 2], [1, 4]]).pvalue
    assert odds == 6.0
    assert checkpoint.fisher((3, 0), (1, 4))[1] == pytest.approx(3.5 * 4.5 / (0.5 * 1.5))


def contrast(name: str, tested: str, **kw) -> Contrast:
    return Contrast(name=name, variable="condition", reference="control", tested=tested,
                    expected_n=(6, 6), **kw)  # fmt: skip


def test_confounding_tests_each_selection_once_with_bh_over_all():
    synthetic = design_counts((("control", "short", 6), ("case", "short", 6)), n_genes=50)
    sheet = synthetic.samples.with_columns(replicate=pl.lit(1, pl.Int32))
    is_case = sheet["condition"] == "case"
    values = pl.DataFrame({
        "sample": sheet["sample"],
        "library_size": np.where(is_case, 2e6, 1e6) + np.arange(12),
        "genes_detected": np.full(12, 100.0),
        "sex_agrees": ["agrees"] * 11 + ["disagrees"],
    })  # fmt: skip
    flags = pl.DataFrame({"sample": [sheet["sample"][11]], "metric": ["sex_check"],
                          "reason": ["r"]})  # fmt: skip
    contrasts = [contrast("one", "case"), contrast("same", "case"),
                 contrast("pair", "case", paired=True, design="~donor + group")]  # fmt: skip
    reg = reg_of(metric("library_size"), metric("genes_detected"), metric("sex_check"))
    table = checkpoint.confounding(values, flags, sheet, contrasts, {}, reg)
    assert table.filter(pl.col("metric") == "library_size")["contrasts"].to_list() == ["one; same"]
    size = table.filter(pl.col("metric") == "library_size").row(0, named=True)
    assert (size["n_reference"], size["n_tested"], size["effect"]) == (6, 6, 1.0)
    notes = dict(table.select("metric", "note").rows())
    assert notes["genes_detected"] == "no test: every sample has the same value"
    assert notes[""].startswith("not tested: paired")
    sex = table.filter(pl.col("metric") == "sex_check").row(0, named=True)
    assert sex["effect"] == pytest.approx(1.5 * 6.5 / (5.5 * 0.5))
    tested = table.filter(pl.col("p").is_not_null())
    assert tested["p_adj"].to_list() == pytest.approx(list(adjust(tested["p"].to_numpy(), "padj")))


def test_sample_metrics_on_gathered_counts(tmp_path):
    synthetic = design_counts()
    write_gather_outputs(tmp_path, "demo", synthetic, date(2026, 9, 30))
    read = {w: pl.read_parquet(next((tmp_path / "outputs" / "01_gather").glob(f"*_{w}_v1.parquet")))
            for w in ("samplesheet", "counts", "genes")}  # fmt: skip
    sheet = read["samplesheet"].with_columns(
        sex=pl.when(pl.col("sample") == "S001").then(pl.lit("M")).otherwise(pl.col("sex"))
    )
    genes = read["genes"]
    types = genes.with_columns(gene_type=pl.lit("protein_coding"), chrom=pl.lit("chr1"))
    measured = checkpoint.sample_metrics(read["counts"], sheet, genes, types,
                                         ["cell_type", "condition"], threads=1)  # fmt: skip
    frame, sexes = measured.values, measured.sexes
    assert sexes is not None
    assert measured.scores.columns == ["sample", "pc1", "pc2"]
    assert 0 < measured.share[1] <= measured.share[0] < 1
    assert frame.height == 12
    assert frame.filter(pl.col("sex_agrees") == "disagrees")["sample"].to_list() == ["S001"]
    assert frame["pca_distance"].null_count() == 0
    assert frame["replicate_spearman"].min() > 0.5


def test_no_flags_still_joins():
    """Nothing flagged: the empty flags table joins, and every run passes."""
    values = pl.DataFrame({"sample": ["a", "b", "c"], "genes_detected": [100.0, 101.0, 99.0]})
    samples = values.select("sample", cell=pl.lit("x"))
    _, flags = checkpoint.judge(values, samples, reg_of(metric("genes_detected")), ["cell"])
    assert flags.is_empty() and flags.schema["sample"] == pl.String
    runs = pl.DataFrame({"run_accession": ["r1", "r2", "r3"], "sample": ["a", "b", "c"]})
    assert checkpoint.status(runs, flags, REG, flags)["qc_status"].to_list() == ["pass"] * 3
