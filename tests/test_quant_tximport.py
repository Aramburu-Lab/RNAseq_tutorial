"""ctx's tximport port against R tximport 1.30.0, run in nf-core/rnaseq 3.27.0's own image.

R's tables are in the quant fixture (``scripts/python/make_quant_fixture.py``): five Salmon
runs on the chr22_se reads, and four planted ``quant.sf`` files. The chr22_se runs share one
fragment-length prior, so every sample's effective lengths are equal there and two orders of
operations agree by accident; the planted set's lengths differ by sample, a gene is absent
from one sample, a gene from all and a transcript from the gene table, which tells them apart.

Mutations each test was seen to fail on, planted in ``ctx/quant/tximport.py`` one at a time
(2026-10-08):

- the arithmetic mean for a partly missing gene's length: planted only;
- the mean length taken before the gaps are filled: planted only;
- scaled to each sample's TPM sum, not its count sum: both;
- lengths weighted by counts, not TPM: both;
- ``Length`` read for ``EffectiveLength``: both;
- a transcript without a gene dropped before scaling, not after: planted only.
"""

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from ctx.quant import CountsFrom, QuantError, outdir, tximport

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "quant"
SETS = {
    "chr22_se": (FIXTURE / "rnaseq" / "salmon", ("se_a", "se_b", "se_c", "se_d", "se_e"),
                 "salmon.merged"),
    "planted": (FIXTURE / "planted", ("s1", "s2", "s3", "s4"), "r"),
}  # fmt: skip
# R and Python sum in the same order; what is left is the last bits of a double.
RELATIVE = 1e-12


def tables(name: str) -> tuple[tximport.GeneTables, Path, tuple[str, ...], str]:
    """ctx's gene tables of one set, and where R's are."""
    folder, samples, prefix = SETS[name]
    return tximport.summarize(outdir.quants(folder, samples), outdir.tx2gene(folder)), folder, \
        samples, prefix


@pytest.mark.parametrize("name", list(SETS))
@pytest.mark.parametrize("table", ["gene_counts", "gene_tpm", "gene_lengths",
                                   "gene_counts_length_scaled", "gene_counts_scaled"])
def test_tximport_matches_r(name, table):
    ours, folder, samples, prefix = tables(name)
    values = {
        "gene_counts": ours.counts_from(CountsFrom.NO), "gene_tpm": ours.tpm,
        "gene_lengths": ours.lengths,
        "gene_counts_length_scaled": ours.counts_from(CountsFrom.LENGTH_SCALED_TPM),
        "gene_counts_scaled": ours.counts_from(CountsFrom.SCALED_TPM),
    }[table]  # fmt: skip
    r = pl.read_csv(folder / f"{prefix}.{table}.tsv", separator="\t")
    assert r["gene_id"].to_list() == list(ours.genes)
    expected = r.select(samples).to_numpy()
    np.testing.assert_allclose(values, expected, rtol=RELATIVE, atol=1e-9)


def test_planted_cases_are_there():
    """The planted set holds every case the order test needs, or it proves nothing."""
    ours, *_ = tables("planted")
    tpm = dict(zip(ours.genes, ours.tpm, strict=True))
    assert (tpm["G1"] == 0).sum() == 1  # one gap: the geometric mean fills it
    assert (tpm["G2"] == 0).all()  # no length anywhere: the transcripts' mean fills it
    assert "T9" not in ours.genes and "T9" in ours.rows  # a transcript in no gene
    lengths = outdir.quant(FIXTURE / "planted" / "s1" / "quant.sf")["EffectiveLength"]
    assert not lengths.equals(outdir.quant(FIXTURE / "planted" / "s2" / "quant.sf")
                              ["EffectiveLength"])  # fmt: skip


def test_every_sample_keeps_its_counts():
    """Counts from abundance keep each sample's estimated counts, the ones in a gene included."""
    ours, *_ = tables("chr22_se")
    for mode in CountsFrom:
        np.testing.assert_allclose(ours.counts_from(mode).sum(axis=0), ours.counts.sum(axis=0))


def test_refuses_another_index(tmp_path):
    folder, samples, _ = SETS["planted"]
    other = outdir.quant(folder / "s2" / "quant.sf").reverse()
    quants = [("s1", outdir.quant(folder / "s1" / "quant.sf")), ("s2", other)]
    with pytest.raises(QuantError, match="s2: its quant.sf lists other transcripts"):
        tximport.summarize(quants, outdir.tx2gene(folder))
    with pytest.raises(QuantError, match="no sample"):
        tximport.summarize([], outdir.tx2gene(folder))
