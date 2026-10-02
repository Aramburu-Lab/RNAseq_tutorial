"""M3's tables, built from the readers of ``align`` and ``count`` and nothing else (D-038).

Every number here is a tool's own, read through MultiQC (the logs) or Polars (per-gene
counts), or arithmetic on those. The page, the figures and the files all read these frames.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import polars as pl

from ctx import concordance
from ctx.process import CATEGORIES, STRAND, Aligner, ProcessError, Strandedness, align, count
from ctx.qc import multiqc_io as mqc_io

# MultiQC's columns behind each outcome, summed, and the column of the aligner's input reads.
# A HISAT2 pair aligned discordantly once is placed once, so it counts as unique.
OUTCOME_COLUMNS: dict[tuple[Aligner, bool], tuple[str, dict[str, tuple[str, ...]]]] = {
    (Aligner.STAR, False): ("total_reads", {
        "unique": ("uniquely_mapped",), "multi": ("multimapped",),
        "too_many": ("multimapped_toomany",), "mismatches": ("unmapped_mismatches",),
        "too_short": ("unmapped_tooshort",), "other": ("unmapped_other",)}),
    (Aligner.HISAT2, False): ("unpaired_total", {
        "unique": ("unpaired_aligned_one",), "multi": ("unpaired_aligned_multi",),
        "unmapped": ("unpaired_aligned_none",)}),
    (Aligner.HISAT2, True): ("paired_total", {
        "unique": ("paired_aligned_one", "paired_aligned_discord_one"),
        "multi": ("paired_aligned_multi",), "unmapped": ("paired_aligned_none",)}),
}  # fmt: skip
OUTCOME_COLUMNS[(Aligner.STAR, True)] = OUTCOME_COLUMNS[(Aligner.STAR, False)]
# STAR's alignment details, kept in the alignment table; null under HISAT2.
STAR_DETAILS = {"mismatch_rate": "mismatch_rate", "splices": "num_splices",
                "splices_annotated": "num_annotated_splices",
                "splices_noncanonical": "num_noncanonical_splices"}  # fmt: skip
# What the page derives from STAR's splice counts; null under HISAT2, whose summary has none.
SPLICING = ("splices_per_100_unique", "splices_annotated_share")
# Gene types the scatters colour by; the rest are "other".
GENE_TYPES = ("protein_coding", "lncRNA")
OTHER = "other"
COUNTERS = ("s0", "s1", "s2", "O", "O_fraction", "M")
# Counters whose Assigned counts alignments, not reads: under -M every reported alignment of a
# multimapper is assigned, so a read counts once per locus (se_c under STAR: 708 Assigned where
# at most 628 + 49 = 677 reads could be), and no share of the aligner's input reads follows.
ALIGNMENT_COUNTERS = frozenset({"M"})
# Each comparison: (a, b), d = b - a. "called" is plain featureCounts at the sample's strand.
COMPARISONS = {"tool": ("star", "called"), "strand": ("s0", "called"), "overlap": ("called", "O"),
               "fraction": ("O", "O_fraction"), "multimap": ("called", "M")}  # fmt: skip
FLAGSTAT = ("total_passed", "mapped_passed", "mapped_passed_pct", "secondary_passed",
            "supplementary_passed")  # fmt: skip


def aligner_rows(run: mqc_io.MultiqcRun, alignment: align.Alignment) -> dict[str, dict[str, float]]:
    """Each sample's row of MultiQC's aligner table, as numbers, with every column it needs.

    Raises:
        ProcessError: If a column or a sample is absent: a zero from a column the table does
            not carry is not a result.
    """
    given, kinds = OUTCOME_COLUMNS[(alignment.aligner, alignment.paired)]
    needed = {given, *(c for columns in kinds.values() for c in columns)}
    if alignment.aligner is Aligner.STAR:
        needed |= set(STAR_DETAILS.values())
    table = mqc_io.table(run, str(alignment.aligner))
    if missing := sorted(needed - set(table.columns)):
        raise ProcessError(f"MultiQC's {alignment.aligner} table has no {missing}")
    if absent := sorted(set(alignment.samples) - set(table["Sample"])):
        raise ProcessError(f"no {alignment.aligner} report for {absent}")
    rows = mqc_io.numeric(table, needed).filter(pl.col("Sample").is_in(alignment.samples))
    return {r["Sample"]: r for r in rows.iter_rows(named=True)}


def outcome_reads(row: Mapping[str, float], aligner: Aligner, paired: bool) -> dict[str, int]:
    """A sample's reads per outcome (``CATEGORIES``) and its input reads, under ``"input"``.

    Raises:
        ProcessError: If the outcomes do not add up to the input: every read has one outcome.
    """
    given, kinds = OUTCOME_COLUMNS[(aligner, paired)]
    reads = {k: int(sum(row[c] for c in columns)) for k, columns in kinds.items()}
    if sum(reads.values()) != int(row[given]):
        raise ProcessError(f"the {aligner} outcomes {reads} add up to {sum(reads.values())}, not "
                           f"the {int(row[given])} input reads")  # fmt: skip
    return {**reads, "input": int(row[given])}


def outcomes_long(reads: Mapping[str, Mapping[str, int]], aligner: Aligner) -> pl.DataFrame:
    """One row per bar label and outcome: reads, input reads and percent of input."""
    rows = [(label, k, r[k], r["input"]) for label, r in reads.items() for k in CATEGORIES[aligner]]
    frame = pl.DataFrame(rows, schema={"label": pl.String, "outcome": pl.String,
                                       "reads": pl.Int64, "input_reads": pl.Int64},
                         orient="row")  # fmt: skip
    return frame.with_columns(percent=pl.col("reads") / pl.col("input_reads") * 100)


def _percents(reads: Mapping[str, int], names: Sequence[str]) -> dict[str, float | None]:
    """Each named outcome as a percent of input reads; None where the aligner has none."""
    return {f"{k}_pct": (100 * reads[k] / reads["input"] if k in reads and reads["input"]
                         else None) for k in names}  # fmt: skip


def shared(reads: Mapping[str, int]) -> dict[str, int]:
    """Unique, multi and unmapped: the outcomes both aligners report, unmapped as the rest."""
    unique, multi = reads["unique"], reads["multi"]
    return {"unique": unique, "multi": multi, "unmapped": reads["input"] - unique - multi,
            "input": reads["input"]}  # fmt: skip


def alignment_table(rows: Mapping[str, Mapping[str, float]], alignment: align.Alignment,
                    files: pl.DataFrame, fastp: pl.DataFrame,
                    flags: Mapping[str, Mapping[str, float]]) -> pl.DataFrame:  # fmt: skip
    """One row per sample at the main setting: outcomes, STAR's details, fastp, flagstat.

    The columns are the same under both aligners, STAR's own null under HISAT2, so a later
    judge (M2b) reads one table shape. fastp's reads before filtering are the realised count
    of ``--reads-per-run``, every read without it.
    """
    star_only = CATEGORIES[Aligner.STAR][2:]
    by_sample = (files.unique("run_accession", keep="first").join(fastp, on="run_accession")
                 .group_by("sample").agg(pl.col("raw_reads", "passed_reads").sum()))  # fmt: skip
    raw = dict(by_sample.select("sample", "raw_reads").rows())
    passed = dict(by_sample.select("sample", "passed_reads").rows())
    runs = dict(files.unique("run_accession", keep="first", maintain_order=True)
                .group_by("sample", maintain_order=True)
                .agg(pl.col("run_accession").str.join(",")).rows())  # fmt: skip
    out = []
    for sample in alignment.samples:
        row, flag = rows[sample], flags[sample]
        reads = outcome_reads(row, alignment.aligner, alignment.paired)
        both = shared(reads)
        extra = flag["secondary_passed"] + flag["supplementary_passed"]
        records = flag["total_passed"] - extra
        out.append({
            "sample": sample, "run_accessions": runs[sample], "input_reads": reads["input"],
            **{k: both[k] for k in ("unique", "multi", "unmapped")},
            **_percents(both, ("unique", "multi", "unmapped")),
            **{k: reads.get(k) for k in star_only}, **_percents(reads, star_only),
            **{k: (row[c] if alignment.aligner is Aligner.STAR else None)
               for k, c in STAR_DETAILS.items()},
            **splicing(row, reads["unique"], alignment.aligner),
            "fastp_raw": int(raw[sample]), "fastp_passed": int(passed[sample]),
            "records": int(flag["total_passed"]),
            "mapped_records_pct": flag["mapped_passed_pct"],
            "primary_mapped_pct": 100 * (flag["mapped_passed"] - extra) / records if records
            else None,
        })  # fmt: skip
    # Typed whatever the aligner, so a HISAT2 table's STAR-only nulls keep STAR's types.
    return pl.DataFrame(out, infer_schema_length=None).with_columns(
        pl.col(star_only).cast(pl.Int64),
        pl.col([*(f"{k}_pct" for k in star_only), *STAR_DETAILS, *SPLICING]).cast(pl.Float64),
    )  # fmt: skip


def splicing(row: Mapping[str, float], unique: int, aligner: Aligner) -> dict[str, float | None]:
    """STAR's splices per 100 uniquely mapped reads and the share at annotated junctions.

    HISAT2's summary counts no splice, so both are null under it, as STAR's own columns are.
    """
    if aligner is not Aligner.STAR:
        return dict.fromkeys(SPLICING)
    splices = row[STAR_DETAILS["splices"]]
    return {"splices_per_100_unique": 100 * splices / unique if unique else None,
            "splices_annotated_share": row[STAR_DETAILS["splices_annotated"]] / splices
            if splices else None}  # fmt: skip


def flagstat(run: mqc_io.MultiqcRun) -> dict[str, dict[str, float]]:
    """MultiQC's flagstat rows by sample, the columns the alignment table reads checked first.

    Raises:
        ProcessError: If a column is absent.
    """
    table = mqc_io.table(run, "samtools_flagstat")
    if missing := [c for c in FLAGSTAT if c not in table.columns]:
        raise ProcessError(f"MultiQC's flagstat table has no {missing}")
    return {r["Sample"]: r for r in mqc_io.numeric(table, FLAGSTAT).iter_rows(named=True)}


def strandedness(counted: pl.DataFrame, calls: pl.DataFrame, used: Mapping[str, Strandedness],
                 star: pl.DataFrame | None, salmon: pl.DataFrame | None
                 ) -> pl.DataFrame:  # fmt: skip
    """One row per sample: what featureCounts assigned at each -s, the call, and its cross-checks.

    Args:
        counted: ``count.assigned`` over the plain settings.
        calls: ``count.strand_calls``.
        used: The strand each sample was counted at.
        star: STAR's gene reads per ReadsPerGene column (sample, unstranded, forward, reverse),
            or None under HISAT2.
        salmon: ``context.expected``'s table, or None without ``--compare-dataset``.
    """
    wide = counted.filter(pl.col("setting").is_in(["s0", "s1", "s2"])).pivot(
        "setting", index="sample", values=["assigned", "no_feature", "ambiguous"],
        maintain_order=True)  # fmt: skip
    frame = wide.join(calls.select("sample", "forward_share", "call"), on="sample").with_columns(
        used=pl.col("sample").replace_strict({s: str(v) for s, v in used.items()}))  # fmt: skip
    columns = {f"star_{c}": pl.lit(None, pl.Int64) for c in ("unstranded", "forward", "reverse")}
    if star is not None:
        frame = frame.join(star.rename({c: f"star_{c}" for c in star.columns[1:]}), on="sample")
    else:
        frame = frame.with_columns(**columns)
    if salmon is None:
        return frame.with_columns(compare_bias=pl.lit(None, pl.Float64),
                                  compare_call=pl.lit(None, pl.String),
                                  agrees=pl.lit(None, pl.Boolean))  # fmt: skip
    return frame.join(salmon.select("sample", compare_bias="salmon_strand_bias",
                                    compare_call="salmon_call"), on="sample").with_columns(
        agrees=pl.col("call") == pl.col("compare_call"))  # fmt: skip


def star_columns(alignment: align.Alignment) -> pl.DataFrame:
    """STAR's reads on genes per ReadsPerGene column, one row per sample."""
    rows = []
    for sample in alignment.samples:
        genes = count.star_genes(alignment.folder / f"{sample}.ReadsPerGene.out.tab")
        rows.append({"sample": sample, **{c: int(genes[c].sum()) for c in genes.columns[1:]}})
    return pl.DataFrame(rows)


def matrices(counts: count.Counts) -> dict[str, pl.DataFrame]:
    """Each setting's gene by BAM matrix, read by its one reader."""
    return {s: count.fc_genes(counts.matrix(s), fraction=o.fraction)
            for s, o in counts.settings.items()}  # fmt: skip


def counting(counted: pl.DataFrame, tables: Mapping[str, pl.DataFrame],
             inputs: Mapping[str, int]) -> pl.DataFrame:  # fmt: skip
    """Per sample and counter: Assigned, the matrix's sum and nonzero genes, Assigned % of input.

    Under -M Assigned counts alignments (``ALIGNMENT_COUNTERS``), so its share of input reads
    is null rather than a share of reads that overstates them.
    """
    rows = []
    for setting in COUNTERS:
        said = dict(counted.filter(pl.col("setting") == setting)
                    .select("sample", "assigned").rows())  # fmt: skip
        for sample, n in inputs.items():
            genes = tables[setting][align.bam_name(sample)]
            rows.append({"sample": sample, "counter": setting, "assigned": said[sample],
                         "gene_count_sum": float(genes.sum()),
                         "genes_nonzero": int((genes > 0).sum()),
                         "assigned_pct": 100 * said[sample] / n
                         if n and setting not in ALIGNMENT_COUNTERS else None})  # fmt: skip
    return pl.DataFrame(rows)


def at_called(tables: Mapping[str, pl.DataFrame], used: Mapping[str, Strandedness]
              ) -> pl.DataFrame:  # fmt: skip
    """Plain featureCounts at each sample's own strand: gene_id, then one column per sample."""
    columns = [tables[f"s{STRAND[strand][1]}"][align.bam_name(s)].alias(s)
               for s, strand in used.items()]  # fmt: skip
    return tables["s0"].select("gene_id").with_columns(columns)


def gene_counts(called: pl.DataFrame, genes: pl.DataFrame) -> pl.DataFrame:
    """The called counts with each gene's name and type, as Int32 where every count fits.

    Raises:
        ProcessError: If a count does not fit Int32 (2,147,483,647).
    """
    samples = called.columns[1:]
    if called.select(pl.max_horizontal(samples).max()).item() > 2**31 - 1:
        raise ProcessError("a gene count exceeds Int32; the table would truncate it")
    frame = genes.select("gene_id", "gene_name", "gene_type").join(called, on="gene_id")
    return frame.with_columns(pl.col(samples).cast(pl.Int32))


def comparisons(tables: Mapping[str, pl.DataFrame], used: Mapping[str, Strandedness],
                star: align.Alignment | None) -> pl.DataFrame:  # fmt: skip
    """The long form behind every comparison: comparison, sample, gene_id, a, b.

    Counts split by --fraction stay fractions; the rest are whole reads, as floats.
    """
    frames = []
    for name, (left, right) in COMPARISONS.items():
        if "star" in (left, right) and star is None:
            continue
        for sample, strand in used.items():
            a, b = (_column(tables, side, sample, strand, star) for side in (left, right))
            # Every table holds exactly the annotation's genes (I4), so the join loses none.
            joined = a.join(b, on="gene_id", how="inner", suffix="_b")
            frames.append(joined.select(
                comparison=pl.lit(name), sample=pl.lit(sample), gene_id="gene_id", a="count",
                b="count_b"))  # fmt: skip
    return pl.concat(frames)


def _column(tables: Mapping[str, pl.DataFrame], side: str, sample: str, strand: Strandedness,
            star: align.Alignment | None) -> pl.DataFrame:  # fmt: skip
    """One side of a comparison for one sample: gene_id and its count, as a float."""
    if side == "star":
        assert star is not None  # noqa: S101 - comparisons skips the tool comparison without STAR
        genes = count.star_genes(star.folder / f"{sample}.ReadsPerGene.out.tab")
        return genes.select("gene_id", count=pl.col(str(strand)).cast(pl.Float64))
    setting = f"s{STRAND[strand][1]}" if side == "called" else side
    return tables[setting].select("gene_id", count=pl.col(align.bam_name(sample))
                                  .cast(pl.Float64))  # fmt: skip


@dataclass(frozen=True)
class Agreements:
    """Each comparison's agreement per sample and on gene sums, and its discordant genes."""

    concordance: pl.DataFrame
    discordant: pl.DataFrame
    sums: pl.DataFrame


def agreements(long: pl.DataFrame, genes: pl.DataFrame) -> Agreements:
    """``concordance.per_unit`` on every comparison, and the genes whose summed counts differ.

    Args:
        long: ``comparisons``' table.
        genes: gene_id, gene_name, gene_type of the annotation's genes.

    Returns:
        The concordance table (comparison, unit, the fields of ``Agreement``), the discordant
        genes over ``all`` and every gene's sums per comparison, typed for the scatters.
    """
    rows, sums = [], []
    for (name,), part in long.group_by("comparison", maintain_order=True):
        rows.append(concordance.per_unit(part, "a", "b", unit="sample", gene="gene_id")
                    .with_columns(comparison=pl.lit(name)))  # fmt: skip
        sums.append(part.group_by("gene_id", maintain_order=True)
                    .agg(pl.col("a").sum(), pl.col("b").sum())
                    .with_columns(comparison=pl.lit(name)))  # fmt: skip
    typed = genes.select(
        "gene_id",
        "gene_name",
        "gene_type",
        type_group=pl.when(pl.col("gene_type").is_in(GENE_TYPES))
        .then("gene_type")
        .otherwise(pl.lit(OTHER)),
    )
    summed = pl.concat(sums).join(typed, on="gene_id", how="left")
    discordant = summed.filter(pl.col("a") != pl.col("b")).select(
        "comparison", "gene_id", "gene_name", "gene_type", "a", "b",
        d=(pl.col("b") + 1).log(2) - (pl.col("a") + 1).log(2))  # fmt: skip
    table = pl.concat(rows).rename({"sample": "unit"})
    return Agreements(table.select("comparison", *table.columns[:-1]), discordant, summed)


def incomplete(main: pl.Series, cut: pl.DataFrame, called: pl.DataFrame,
               unique_pct: tuple[float, float]) -> pl.DataFrame:  # fmt: skip
    """The featured sample on the whole and on the half reference, gene by gene, summed.

    Args:
        main: The featured sample's called counts on the whole reference, by ``called``'s rows.
        cut: gene_id and the featured sample's count on the half reference.
        called: gene_id of the whole reference's genes, in ``main``'s order.
        unique_pct: Its uniquely mapped percent on the whole and on the half reference.

    Returns:
        One row: genes kept and removed; reads on kept genes on each reference; reads kept
        genes gained and lost on the half; reads on the removed genes on the whole.
    """
    whole = called.select("gene_id").with_columns(main=main.cast(pl.Float64))
    both = whole.join(cut.rename({cut.columns[1]: "cut"}), on="gene_id", how="left")
    kept = both.filter(pl.col("cut").is_not_null())
    removed = both.filter(pl.col("cut").is_null())
    change = kept["cut"] - kept["main"]
    return pl.DataFrame([{
        "genes_kept": kept.height, "genes_removed": removed.height,
        "kept_reads_main": float(kept["main"].sum()), "kept_reads_cut": float(kept["cut"].sum()),
        "gained": float(change.clip(lower_bound=0).sum()),
        "lost": float((-change).clip(lower_bound=0).sum()),
        "removed_reads_main": float(removed["main"].sum()),
        "unique_pct_main": unique_pct[0], "unique_pct_cut": unique_pct[1],
    }])  # fmt: skip
