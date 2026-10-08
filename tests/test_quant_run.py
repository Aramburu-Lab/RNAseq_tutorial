"""ctx quant on the quant fixture: the door, the checks, the planted truth and the whole page.

Each test names the mutation of ``src/ctx/quant`` it was seen to fail on (2026-10-08).
"""

import math
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from quant_root import DATASET, SALMON, quant_root

from ctx import runtime
from ctx.io.naming import lint
from ctx.quant import CountsFrom, QuantError, checks, outdir, tables, tximport
from ctx.quant.run import Options, resolve

SE = ("se_a", "se_b", "se_c", "se_d", "se_e")
FIXTURE = Path(__file__).resolve().parent / "fixtures"
# SEC14L2, ENSG00000249590 and MTFP1 overlap on one strand: reads in the overlap fit all
# three, so only their sum is planted exactly.
TRIO = ("ENSG00000100003.18", "ENSG00000249590.7", "ENSG00000242114.8")


@pytest.fixture
def built(quant_built):
    """The session's one run of ctx quant on the fixture (``conftest.quant_built``)."""
    return quant_built


def newest(root: Path, what: str) -> pl.DataFrame:
    """One table ctx quant wrote."""
    return pl.read_parquet(next((root / "outputs" / "04_quant").glob(f"*_{what}_v1.parquet")))


def test_end_to_end(built):
    """Every figure and table, each named by the rules, the page first.

    Mutation: the isoform figure's group dropped from the page (``q1_groups`` returning the
    gene group alone): no isoforms files.
    """
    root, _, written = built
    assert written[0].name.endswith(f"_quant_{DATASET}_overview_v1.html")
    names = {p.name for p in written}
    for what in ("fates", "correlation", "detection", "units", "estimates", "isoforms"):
        assert any(f"_{what}_v1_white.pdf" in n for n in names), what
    for what in ("fates", "library", "sums", "gene_counts", "fold_changes", "isoforms"):
        assert any(f"_{what}_v1.parquet" in n for n in names), what
    _, problems = lint(root / "outputs" / "04_quant")
    assert problems == []
    page = written[0].read_text()
    assert "Not yet drawn" in page


def test_gene_counts_follow_the_option(built):
    """The gene table holds the mode asked for, equal to R's table of that mode.

    Mutation: the table written from ``port.counts`` whatever ``--counts-from-abundance``
    says: it then equals R's gene_counts, not the length-scaled table.
    """
    root, *_ = built
    ours = newest(root, "gene_counts")
    r = pl.read_csv(SALMON / "salmon.merged.gene_counts_length_scaled.tsv", separator="\t")
    assert ours["gene_id"].to_list() == r["gene_id"].to_list()
    np.testing.assert_allclose(ours.select(SE).to_numpy(), r.select(SE).to_numpy(), rtol=1e-12)


def test_fold_change_thinning_keeps_the_far_genes(tmp_path):
    """The genes F6's lead counts, two-fold or more apart, are always drawn on the page.

    Gene g0 keeps its reads in the tested group while its length falls to a quarter, as an
    isoform switch does: its estimated fold change is 0 and its length-scaled one about 2.
    Mutation: ``keep`` set false for every gene: g0 would be thinned with the rest.
    """
    from datetime import date  # noqa: PLC0415

    from ctx.de.run import Writer  # noqa: PLC0415
    from ctx.quant import blocks  # noqa: PLC0415

    genes, samples = [f"g{i:02d}" for i in range(50)], ("r1", "r2", "r3", "t1", "t2", "t3")
    counts = np.full((50, 6), 100.0)
    lengths = np.full((50, 6), 1000.0)
    lengths[0, 3:] = 250.0
    rate = counts / lengths
    port = tximport.GeneTables(tuple(genes), samples, counts, 1e6 * rate / rate.sum(axis=0),
                               lengths, np.ones(50, dtype=bool))  # fmt: skip
    used = pl.DataFrame({"sample": list(samples), "group": ["reference"] * 3 + ["tested"] * 3})
    annotation = pl.DataFrame({"gene_id": genes, "gene_name": genes,
                               "gene_type": ["protein_coding"] * 50})  # fmt: skip
    writer = Writer(tmp_path / "04_quant", "planted", date(2026, 10, 8), {})
    _, frame = blocks.estimates_figure(writer, port, used, annotation, "planted")
    far = frame.filter(pl.col("difference").abs() >= 1)
    assert far["gene_id"].to_list() == ["g00"]
    assert frame.filter("keep")["gene_id"].to_list() == ["g00"]


def test_planted_truth(built):
    """Salmon's library calls and gene totals are what chr22_se planted.

    SR in four samples and U in se_e. In se_a, se_b and se_d every gene equals its planted
    reads but the overlapping trio, which equals its planted total.
    Mutation: NumReads read as TPM: the run itself refuses, its port off nf-core's tables.
    """
    root, *_ = built
    library = newest(root, "library")
    assert dict(library.select("sample", "library_type").iter_rows()) == {
        "se_a": "SR", "se_b": "SR", "se_c": "SR", "se_d": "SR", "se_e": "U"}  # fmt: skip
    truth = (pl.read_csv(FIXTURE / "chr22_se" / "truth.tsv", separator="\t")
             .group_by("sample", "gene_id").agg(pl.col("n_reads").sum()))  # fmt: skip
    port = tximport.summarize(outdir.quants(SALMON, SE), outdir.tx2gene(SALMON))
    ours = port.frame(port.counts).unpivot(index="gene_id", variable_name="sample",
                                           value_name="estimated")  # fmt: skip
    joined = ours.join(truth, on=["sample", "gene_id"], how="left").fill_null(0)
    for sample in ("se_a", "se_b", "se_d"):
        mine = joined.filter(pl.col("sample") == sample)
        alone = mine.filter(~pl.col("gene_id").is_in(TRIO))
        np.testing.assert_allclose(alone["estimated"], alone["n_reads"], atol=0.01)
        trio = mine.filter(pl.col("gene_id").is_in(TRIO))
        assert trio["estimated"].sum() == pytest.approx(trio["n_reads"].sum(), abs=0.01)
    assert joined.filter((pl.col("sample") == "se_a") & pl.col("gene_id").is_in(TRIO))[
        "n_reads"].sum() == 458  # fmt: skip


def salmon_table(**fates: list[int]) -> pl.DataFrame:
    """A Salmon table of two samples: processed 100 each and the given fates."""
    base = {"sample": ["s1", "s2"], "num_processed": [100, 100], "num_mapped": [80, 70],
            "num_decoy_fragments": [5, 10], "num_fragments_filtered_vm": [3, 4]}  # fmt: skip
    return pl.DataFrame({**base, **fates})


def test_fates():
    """Four disjoint fates that add up to what Salmon processed.

    Mutation: the decoy term dropped from the remainder: no_hit is then 17 and 26, and the
    fates add to more than processed.
    """
    found = tables.fates(salmon_table())
    assert found["no_hit"].to_list() == [12, 16]
    assert found.select(pl.sum_horizontal(tables.FATES)).to_series().to_list() == [100, 100]
    assert checks.fate_problems(found) == []
    too_many = tables.fates(salmon_table(num_mapped=[95, 70]))
    assert "s1" in checks.fate_problems(too_many)[0]


def test_sums():
    """A length-scaled sum off the mapped fragments by more than quant.sf's rounding stops it.

    Mutation: length-scaled counts scaled to each sample's TPM sum (``counts_from`` scaling
    by ``abundance``): the fixture's sums are then a million, not the mapped fragments.
    """
    port = tximport.summarize(outdir.quants(SALMON, SE), outdir.tx2gene(SALMON))
    scaled = dict(zip(SE, port.counts_from(CountsFrom.LENGTH_SCALED_TPM).sum(axis=0),
                      strict=True))  # fmt: skip
    mapped = pl.DataFrame({"sample": list(SE), "mapped": [1000, 1000, 873, 1000, 1000]})
    targets = pl.DataFrame({"sample": list(SE), "targets": [36] * 5})
    assert checks.sum_problems(tables.sums(mapped, scaled), targets) == []
    off = {**scaled, "se_c": scaled["se_c"] + 1}
    assert "se_c" in checks.sum_problems(tables.sums(mapped, off), targets)[0]


def donor_rho(seed: int = 3) -> tuple[np.ndarray, list[str], list[str]]:
    """Six donors of six samples: rho 0.9 within a donor, 0.3 between, no cell-type effect.

    Each donor has five samples of one cell type and one of the other, the majority type
    alternating by donor, so a cell type's pairs are mostly one donor's pairs.
    """
    rng = np.random.default_rng(seed)
    donors = [f"d{i}" for i in range(6) for _ in range(6)]
    labels = [("a" if (i % 2 == 0) == (j < 5) else "b") for i in range(6) for j in range(6)]
    same = np.equal.outer(donors, donors)
    noise = rng.normal(0, 0.01, (36, 36))
    rho = np.where(same, 0.9, 0.3) + (noise + noise.T) / 2
    np.fill_diagonal(rho, 1.0)
    return rho, labels, donors


def test_restricted_permutation():
    """A donor effect that a free permutation calls a cell-type effect; within donors it is not.

    Mutation: labels permuted across all samples, the strata ignored: p falls below 0.05.
    """
    rho, labels, donors = donor_rho()
    within = tables.within_minus_between(rho, labels, donors, permutations=999)
    free = tables.within_minus_between(rho, labels, ["all"] * len(labels), permutations=999)
    assert within is not None and free is not None
    assert free.p < 0.05
    assert within.p > 0.2
    assert within.strata == 6
    assert tables.within_minus_between(rho, ["a"] * 36, donors) is None
    one_each = tables.within_minus_between(rho, labels, [f"s{i}" for i in range(36)])
    assert one_each is not None and math.isnan(one_each.p)


def test_door(tmp_path):
    """Every input problem is listed together, and nothing is written.

    Mutation: ``resolve`` raising at the first problem: the message names one, not three.
    """
    root = tmp_path / "root"
    spec = quant_root(root)
    next((root / "resources").rglob("*.gtf")).unlink()
    empty = tmp_path / "salmon"
    empty.mkdir()
    settings = runtime.Settings("laptop", root, 1, None, 0.0)
    options = Options(spec, tmp_path / "nowhere", empty)
    with pytest.raises(QuantError) as caught:
        resolve(settings, DATASET, options)
    message = str(caught.value)
    assert "the annotation is not at" in message
    assert "nowhere is not a folder" in message
    assert "no quant.sf for" in message
    assert not (root / "outputs" / "04_quant").exists()


def test_a_missing_salmon_column_is_refused(tmp_path):
    """A zero from a column MultiQC's table does not carry is not a result.

    Mutation: the column check in ``salmon_table`` removed: a KeyError, not a QuantError.
    """
    from ctx.qc.multiqc_io import MultiqcRun  # noqa: PLC0415

    data = tmp_path / "multiqc_data"
    data.mkdir()
    (data / "multiqc_salmon.txt").write_text("Sample\tnum_processed\tnum_mapped\nse_a\t10\t8\n")
    (data / "multiqc_general_stats.txt").write_text("Sample\nse_a\n")
    with pytest.raises(QuantError, match="has no"):
        outdir.salmon_table(MultiqcRun(data, "1.35"))


def test_ceiling():
    """The matrix may count no gene the annotation lacks."""
    assert checks.ceiling_problems(["g1", "g2"], ["g1", "g2", "g3"]) == []
    assert "1 counted genes" in checks.ceiling_problems(["g1", "g9"], ["g1"])[0]


def test_port_differing_from_nf_core_is_refused(tmp_path):
    """A merged gene table that the port does not reproduce stops the run.

    Mutation: ``tximport_problems`` comparing only gene order: the edited value passes.
    """
    folder = tmp_path / "salmon"
    folder.mkdir()
    for path in SALMON.glob("salmon.merged.gene_*.tsv"):
        (folder / path.name).write_bytes(path.read_bytes())
    port = tximport.summarize(outdir.quants(SALMON, SE), outdir.tx2gene(SALMON))
    assert checks.tximport_problems(port, folder) == []
    table = folder / "salmon.merged.gene_tpm.tsv"
    frame = pl.read_csv(table, separator="\t")
    frame.with_columns(pl.col("se_a") * 1.01).write_csv(table, separator="\t")
    assert any(p.startswith("gene_tpm") for p in checks.tximport_problems(port, folder))
