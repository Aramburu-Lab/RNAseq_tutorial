"""ctx.process.tables on hand-made rows: outcomes, comparisons, agreement and the half reference."""

from pathlib import Path

import polars as pl
import pytest

from ctx.process import Aligner, ProcessError, Strandedness, tables

STAR_ROW = {"total_reads": 1000.0, "uniquely_mapped": 921.0, "multimapped": 49.0,
            "multimapped_toomany": 0.0, "unmapped_mismatches": 0.0, "unmapped_tooshort": 23.0,
            "unmapped_other": 7.0}  # fmt: skip


def test_star_outcomes_cover_every_input_read():
    """se_c's STAR numbers: six outcomes adding up to its 1,000 reads; the shared three put
    too many loci, mismatches, too short and other under unmapped.

    Mutation caught: an outcome left out of OUTCOME_COLUMNS (the sum check refuses the row)."""
    reads = tables.outcome_reads(STAR_ROW, Aligner.STAR, paired=False)
    assert reads == {"unique": 921, "multi": 49, "too_many": 0, "mismatches": 0, "too_short": 23,
                     "other": 7, "input": 1000}  # fmt: skip
    assert tables.shared(reads) == {"unique": 921, "multi": 49, "unmapped": 30, "input": 1000}


def test_outcomes_that_miss_reads_are_refused():
    """A row whose outcomes add up to 999 of 1,000 reads stops the page.

    Mutation caught: the sum check dropped (a read with no outcome drawn as if counted)."""
    with pytest.raises(ProcessError, match="add up to 999, not the 1000"):
        tables.outcome_reads({**STAR_ROW, "unmapped_other": 6.0}, Aligner.STAR, paired=False)


def test_a_hisat2_pair_aligned_discordantly_once_is_unique():
    """The paired fixture's HISAT2 row: 2,493 concordant once and 2 discordant once are unique,
    2 multi, 3 neither: 2,500 pairs.

    Mutation caught: discordant pairs left out (the sum check refuses the row)."""
    row = {"paired_total": 2500.0, "paired_aligned_one": 2493.0, "paired_aligned_multi": 2.0,
           "paired_aligned_discord_one": 2.0, "paired_aligned_none": 3.0}  # fmt: skip
    reads = tables.outcome_reads(row, Aligner.HISAT2, paired=True)
    assert reads == {"unique": 2495, "multi": 2, "unmapped": 3, "input": 2500}


def matrix(**columns: list[float]) -> pl.DataFrame:
    """A featureCounts matrix of three genes, one column per BAM."""
    return pl.DataFrame({"gene_id": ["g1", "g2", "g3"], **columns})


def test_each_sample_is_compared_at_its_own_strand():
    """se_a reverse and se_e unstranded: the strand comparison's b is se_a's -s 2 and se_e's
    -s 0, and without STAR the tool comparison is left out.

    Mutation caught: every sample read at one strand (se_e's b taken from -s 2)."""
    mats = {"s0": matrix(**{"a.bam": [5.0, 0.0, 1.0], "e.bam": [7.0, 1.0, 0.0]}),
            "s1": matrix(**{"a.bam": [0.0, 0.0, 0.0], "e.bam": [3.0, 0.0, 0.0]}),
            "s2": matrix(**{"a.bam": [6.0, 0.0, 1.0], "e.bam": [4.0, 1.0, 0.0]}),
            "O": matrix(**{"a.bam": [9.0, 2.0, 1.0], "e.bam": [9.0, 1.0, 0.0]}),
            "O_fraction": matrix(**{"a.bam": [6.5, 0.5, 1.0], "e.bam": [7.0, 1.0, 0.0]}),
            "M": matrix(**{"a.bam": [6.0, 0.0, 3.0], "e.bam": [7.0, 1.0, 2.0]})}  # fmt: skip
    used = {"a": Strandedness.REVERSE, "e": Strandedness.UNSTRANDED}
    long = tables.comparisons(mats, used, None)
    assert set(long["comparison"]) == {"strand", "overlap", "fraction", "multimap"}
    strand = long.filter(pl.col("comparison") == "strand")
    assert strand.filter(pl.col("sample") == "a")["b"].to_list() == [6.0, 0.0, 1.0]
    assert strand.filter(pl.col("sample") == "e")["b"].to_list() == [7.0, 1.0, 0.0]
    multimap = long.filter((pl.col("comparison") == "multimap") & (pl.col("sample") == "e"))
    assert multimap.select("a", "b").rows() == [(7.0, 7.0), (1.0, 1.0), (0.0, 2.0)]


def test_discordant_genes_are_those_whose_sums_differ():
    """Two samples, three genes: g1 agrees in each, g2 differs in each but sums alike, g3 differs.
    Only g3 is discordant, with d on its sums; g2 counts as differing per sample, not in all.

    Mutation caught: discordance judged per sample row, not on the gene sums (g2 listed)."""
    long = pl.DataFrame({"comparison": "multimap", "sample": ["x"] * 3 + ["y"] * 3,
                         "gene_id": ["g1", "g2", "g3"] * 2,
                         "a": [5.0, 2.0, 1.0, 5.0, 0.0, 1.0],
                         "b": [5.0, 0.0, 1.0, 5.0, 2.0, 7.0]})  # fmt: skip
    genes = pl.DataFrame({"gene_id": ["g1", "g2", "g3"], "gene_name": ["A", "B", "C"],
                          "gene_type": ["protein_coding", "lncRNA", "snRNA"]})  # fmt: skip
    found = tables.agreements(long, genes)
    assert found.discordant.select("gene_id", "a", "b").rows() == [("g3", 2.0, 8.0)]
    assert found.discordant["d"].item() == pytest.approx(1.5849625)  # log2(9) - log2(3)
    every = found.concordance.filter(pl.col("unit") == "all").row(0, named=True)
    assert (every["n_genes"], every["n_identical"]) == (3, 2)
    groups = dict(found.sums.select("gene_id", "type_group").rows())
    assert groups == {"g1": "protein_coding", "g2": "lncRNA", "g3": "other"}


def test_the_half_reference_sums_gains_losses_and_removed_reads():
    """Three genes on the whole, two on the half: g1 gains 2, g2 loses 3, g3 is removed with 4.

    Mutations caught: gained and lost swapped; the removed gene's reads counted as lost."""
    called = pl.DataFrame({"gene_id": ["g1", "g2", "g3"], "s": [10, 5, 4]})
    cut = pl.DataFrame({"gene_id": ["g1", "g2"], "s.bam": [12.0, 2.0]})
    row = tables.incomplete(called["s"], cut, called, (99.0, 60.0)).row(0, named=True)
    assert row == {"genes_kept": 2, "genes_removed": 1, "kept_reads_main": 15.0,
                   "kept_reads_cut": 14.0, "gained": 2.0, "lost": 3.0,
                   "removed_reads_main": 4.0, "unique_pct_main": 99.0,
                   "unique_pct_cut": 60.0}  # fmt: skip


def test_a_multimapper_s_alignments_are_no_share_of_reads():
    """se_c under STAR: -M's Assigned is 708 alignments where at most 628 + 49 = 677 reads could
    be assigned, so the M row keeps featureCounts' number and gives no share of input reads;
    the plain rows keep theirs.

    Mutation caught: the M row's share computed like the others (70.8 % of 1,000 reads, though
    at most 67.7 % could be)."""
    counters = tables.COUNTERS
    counted = pl.DataFrame({"sample": ["se_c"] * len(counters), "setting": list(counters),
                            "assigned": [700.0, 49.0, 628.0, 640.0, 628.0, 708.0]})  # fmt: skip
    mats = {
        c: pl.DataFrame({"gene_id": ["g1"], "se_c.bam": [v]})
        for c, v in zip(counters, [700.0, 49.0, 628.0, 700.0, 628.0, 708.0], strict=True)
    }
    table = tables.counting(counted, mats, {"se_c": 1000})
    shares = dict(table.select("counter", "assigned_pct").rows())
    assert shares["M"] is None and shares["s2"] == pytest.approx(62.8)
    assert all(shares[c] is not None for c in counters if c != "M")


def fastp(runs: dict[str, tuple[int, int]]) -> pl.DataFrame:
    """``align.fastp_reads``' table: run, reads before and after fastp."""
    return pl.DataFrame([(r, raw, passed) for r, (raw, passed) in runs.items()],
                        schema=["run_accession", "raw_reads", "passed_reads"],
                        orient="row")  # fmt: skip


def test_the_alignment_table_carries_the_cap_and_star_s_splices():
    """A sample of two runs capped at 500 reads each: fastp_raw sums the capped reads, not the
    passed ones. Under STAR, 310 splices in 921 uniquely mapped reads are 33.7 per 100, 230 of
    them annotated (0.742); HISAT2 has no splice count, so both are null.

    Mutations caught: fastp_raw summed from passed reads; splices per 100 taken over input
    reads; the annotated share inverted."""
    from ctx.process import Aligner, align  # noqa: PLC0415

    files = pl.DataFrame({"run_accession": ["r1", "r2"], "sample": ["se_c", "se_c"],
                          "mate": [1, 1], "path": ["a.fq", "b.fq"]})  # fmt: skip
    row = {**STAR_ROW, "mismatch_rate": 0.3, "num_splices": 310.0,
           "num_annotated_splices": 230.0, "num_noncanonical_splices": 20.0}  # fmt: skip
    flags = {"se_c": {"total_passed": 1049.0, "mapped_passed": 1019.0, "mapped_passed_pct": 97.1,
                      "secondary_passed": 49.0, "supplementary_passed": 0.0}}  # fmt: skip
    star = align.Alignment(Path("."), False, Aligner.STAR, ("se_c",), False)
    table = tables.alignment_table(
        {"se_c": row}, star, files, fastp({"r1": (500, 480), "r2": (500, 470)}), flags
    )
    got = table.row(0, named=True)
    assert (got["fastp_raw"], got["fastp_passed"]) == (1000, 950)
    assert got["splices_per_100_unique"] == pytest.approx(100 * 310 / 921)
    assert got["splices_annotated_share"] == pytest.approx(230 / 310)
    hisat = tables.splicing(row, 921, Aligner.HISAT2)
    assert hisat == {"splices_per_100_unique": None, "splices_annotated_share": None}
