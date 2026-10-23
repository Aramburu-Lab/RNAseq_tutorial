"""Read an nf-core/rnaseq 3.27.0 outdir aligned with STAR, every number nf-core's own (D-116).

M3's route from the cluster's alignment (`ctx process --rnaseq-dir`): STAR, samtools, RSeQC,
Qualimap, dupRadar, Preseq and featureCounts' biotypes as nf-core's own MultiQC wrote them
(``multiqc_io.nfcore``), Salmon's two gene matrices and the run's params and trace. Nothing
is aligned, counted or run again (D-113). Each table's samples are held to the sheet's, and
every column a page reads is checked present before it is read: a zero from a column a table
does not carry is not a result. Nothing here writes into the outdir.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ctx.gather import rnaseq, trace
from ctx.process import ProcessError
from ctx.qc import QcError, fates
from ctx.qc import multiqc_io as mqc_io

# STAR's outcomes for a read, the order the mapping figure stacks them in: MultiQC's STAR
# columns of the reads and of STAR's own percent of input reads (Log.final.out).
OUTCOMES = {
    "unique": ("uniquely_mapped", "uniquely_mapped_percent"),
    "multi": ("multimapped", "multimapped_percent"),
    "too_many": ("multimapped_toomany", "multimapped_toomany_percent"),
    "mismatches": ("unmapped_mismatches", "unmapped_mismatches_percent"),
    "too_short": ("unmapped_tooshort", "unmapped_tooshort_percent"),
    "other": ("unmapped_other", "unmapped_other_percent"),
}
STAR_INPUT = "total_reads"
STAR_DETAILS = {"mismatch_rate": "mismatch_rate", "splices": "num_splices",
                "splices_annotated": "num_annotated_splices",
                "splices_noncanonical": "num_noncanonical_splices"}  # fmt: skip
FLAGSTAT = ("total_passed", "mapped_passed", "mapped_passed_pct", "secondary_passed",
            "supplementary_passed")  # fmt: skip
TRIM = ("total_reads", "reads_written")
# RSeQC infer_experiment as MultiQC plots it, % of the reads sampled, and nf-core's strand check.
STRAND_PLOT = "rseqc_infer_experiment_plot"
STRAND_FRACTIONS = {"Sense": "rseqc_sense", "Antisense": "rseqc_antisense",
                    "Undetermined": "rseqc_undetermined"}  # fmt: skip
STRAND_CHECK = "strand_check_summary_table"
STRAND_CALLS = ("provided", "rseqc_inferred", "rseqc_pct", "salmon_inferred", "status")
READ_DISTRIBUTION = "rseqc_read_distribution_plot"
JUNCTIONS = "rseqc_junction_annotation"
# Each junction kind's count and RSeQC's percent, as events and as distinct junctions.
JUNCTION_KINDS = ("known", "partial_novel", "novel")
JUNCTION_MEASURES = ("events", "junctions")
COVERAGE = "qualimap_gene_coverage_profile_Normalised"
DUPRADAR = "dupradar"
PRESEQ = "preseq_complexity_plot_molecules"
# Preseq's sample names in MultiQC end with its output's suffix; one row is MultiQC's own line.
PRESEQ_SUFFIX = ".lc_extrap"
PRESEQ_PERFECT = "A perfect library where each read is unique"
BIOTYPES = "featurecounts_biotype_plot"
# The two Salmon gene matrices M3 compares: Salmon alone (salmon/) and on STAR's alignment.
MATRICES = {"salmon": rnaseq.MATRIX_FILE, "star_salmon": rnaseq.ALIGNED_MATRIX}
# The trace's processes M3's numbers come from, as Nextflow names them after the last colon.
STEPS = ("STAR_ALIGN", "SAMTOOLS_SORT", "SAMTOOLS_INDEX", "SAMTOOLS_STATS",
         "SAMTOOLS_FLAGSTAT", "SAMTOOLS_IDXSTATS", "PICARD_MARKDUPLICATES", "SALMON_QUANT",
         "RSEQC_INFEREXPERIMENT", "RSEQC_READDISTRIBUTION", "RSEQC_JUNCTIONANNOTATION",
         "QUALIMAP_RNASEQ", "DUPRADAR", "PRESEQ_LCEXTRAP", "SUBREAD_FEATURECOUNTS",
         "CUSTOM_MULTIQCCUSTOMBIOTYPE", "MULTIQC")  # fmt: skip


@dataclass(frozen=True)
class Outdir:
    """One aligned outdir, checked: its MultiQC tables, its params and its provenance.

    Attributes:
        path: The outdir.
        mqc: nf-core's own MultiQC tables.
        samples: The sheet's samples, in sheet order.
        params: The first params file's settings.
        provenance: Release, Salmon, the params that agree and each file's SHA-256.
    """

    path: Path
    mqc: mqc_io.MultiqcRun
    samples: tuple[str, ...]
    params: Mapping[str, object]
    provenance: Mapping[str, str]


def open_outdir(path: Path, samples: Sequence[str]) -> Outdir:
    """An outdir aligned with STAR, its release, params and MultiQC tables checked.

    Raises:
        ProcessError: If it is not nf-core/rnaseq 3.27.0, was not aligned with STAR, or holds
            no MultiQC tables of the pinned MultiQC.
    """
    try:
        provenance = rnaseq.provenance(path)
        if not rnaseq.aligned(path):
            raise ProcessError(f"{path} skipped alignment (skip_alignment true): M3 reads an "
                               f"outdir aligned with STAR, aligner {rnaseq.ALIGNER}")  # fmt: skip
        mqc = mqc_io.nfcore(path)
    except (rnaseq.RnaseqError, QcError) as error:
        raise ProcessError(str(error)) from error
    params = rnaseq.first_params(path)
    return Outdir(path, mqc, tuple(samples), params, provenance)


def _held(table: pl.DataFrame, samples: Sequence[str], what: str,
          names: Mapping[str, str] | None = None) -> pl.DataFrame:  # fmt: skip
    """The table's rows keyed by ``sample``, every sample of the sheet in it and no other.

    ``names`` maps a MultiQC name to its sample, where a tool adds a suffix.
    """
    named = table.with_columns(sample=pl.col("Sample").replace_strict(
        names, default=pl.col("Sample")) if names else pl.col("Sample"))  # fmt: skip
    found = set(named["sample"])
    problems = [f"{s}: in the sheet, not in {what}" for s in samples if s not in found]
    problems += [f"{s}: in {what}, not in the sheet" for s in sorted(found - set(samples))]
    if named["sample"].is_duplicated().any():
        problems.append(f"{what}: a sample appears twice")
    if problems:
        raise ProcessError("nf-core's tables and the sheet differ:\n  " + "\n  ".join(problems))
    order = {s: i for i, s in enumerate(samples)}
    return named.drop("Sample").sort(pl.col("sample").replace_strict(order)).select(
        "sample", pl.exclude("sample"))  # fmt: skip


def _columns(table: pl.DataFrame, needed: Sequence[str], what: str) -> pl.DataFrame:
    """The table with ``needed`` cast to numbers, or refused naming every column it lacks."""
    if missing := [c for c in needed if c not in table.columns]:
        raise ProcessError(f"nf-core's {what} has no {missing}")
    try:
        return mqc_io.numeric(table, needed)
    except pl.exceptions.InvalidOperationError as error:
        raise ProcessError(
            f"nf-core's {what} holds a value that is not a number: {error}"
        ) from error


def _table(outdir: Outdir, module: str) -> pl.DataFrame:
    """A module table of nf-core's MultiQC, ``multiqc_<module>.txt``."""
    try:
        return mqc_io.table(outdir.mqc, module)
    except QcError as error:
        raise ProcessError(str(error)) from error


def star(outdir: Outdir) -> pl.DataFrame:
    """STAR's Log.final.out per sample: input reads, each outcome's reads and STAR's percents.

    Raises:
        ProcessError: If a column or a sample is absent, or a sample's outcomes do not add up
            to its input reads: every read has one outcome.
    """
    needed = [STAR_INPUT, *(c for pair in OUTCOMES.values() for c in pair),
              *STAR_DETAILS.values()]  # fmt: skip
    table = _held(_columns(_table(outdir, "star"), needed, "STAR table"), outdir.samples,
                  "STAR's table")  # fmt: skip
    frame = table.select(
        "sample", input_reads=pl.col(STAR_INPUT).cast(pl.Int64),
        **{k: pl.col(reads).cast(pl.Int64) for k, (reads, _) in OUTCOMES.items()},
        **{f"{k}_pct": pl.col(pct) for k, (_, pct) in OUTCOMES.items()},
        **{k: pl.col(c) for k, c in STAR_DETAILS.items()},
    )  # fmt: skip
    total = pl.sum_horizontal(list(OUTCOMES))
    if (unequal := frame.filter(total != pl.col("input_reads"))["sample"]).len():
        raise ProcessError(f"STAR's outcomes do not add up to its input reads for "
                           f"{unequal.to_list()}")  # fmt: skip
    return frame


def flagstat(outdir: Outdir) -> pl.DataFrame:
    """The flagstat of samtools per sample, on nf-core's duplicate-marked BAM, through MultiQC."""
    table = _columns(_table(outdir, "samtools_flagstat"), FLAGSTAT, "flagstat table")
    return _held(table.select("Sample", *FLAGSTAT), outdir.samples, "samtools' flagstat")


def trimmed(outdir: Outdir) -> pl.DataFrame:
    """Trim Galore's reads in and written per sample: read 1 of a pair, as fragments."""
    table = _columns(_table(outdir, "trim_galore"), TRIM, "Trim Galore table")
    ranks = fates.nfcore_ranks(outdir.samples)
    try:
        frame = mqc_io.by_unit(table.select("Sample", *TRIM), ranks, "sample",
                               fates.FIRST_UNUSED_RANK)  # fmt: skip
    except QcError as error:
        raise ProcessError(str(error)) from error
    if missing := sorted(set(outdir.samples) - set(frame["sample"])):
        raise ProcessError(f"no Trim Galore report for {missing}")
    return frame.select("sample", trim_in=pl.col("total_reads").cast(pl.Int64),
                        trim_out=pl.col("reads_written").cast(pl.Int64))  # fmt: skip


def strand(outdir: Outdir) -> pl.DataFrame:
    """RSeQC infer_experiment's three shares (%) and nf-core's strand check, per sample.

    The call, its share and its pass or fail are nf-core's (``multiqc_strand_check_summary_
    table``); no threshold is applied here.
    """
    try:
        shares = mqc_io.plot_table(outdir.mqc, STRAND_PLOT)
    except QcError as error:
        raise ProcessError(str(error)) from error
    if missing := [c for c in STRAND_FRACTIONS if c not in shares.columns]:
        raise ProcessError(f"nf-core's {STRAND_PLOT} has no {missing}")
    shares = _held(shares.select("Sample", *STRAND_FRACTIONS).rename(STRAND_FRACTIONS),
                   outdir.samples, "RSeQC infer_experiment")  # fmt: skip
    check = _table(outdir, STRAND_CHECK)
    if missing := [c for c in STRAND_CALLS if c not in check.columns]:
        raise ProcessError(f"nf-core's strand check table has no {missing}")
    check = _columns(check.select("Sample", *STRAND_CALLS), ["rseqc_pct"], "strand check table")
    check = _held(check, outdir.samples, "nf-core's strand check")
    return shares.join(check, on="sample", maintain_order="left")


def read_distribution(outdir: Outdir) -> pl.DataFrame:
    """RSeQC read_distribution's tags per category, as MultiQC plots them, in long form."""
    try:
        table = mqc_io.plot_table(outdir.mqc, READ_DISTRIBUTION)
    except QcError as error:
        raise ProcessError(str(error)) from error
    held = _held(table, outdir.samples, "RSeQC read_distribution")
    return held.unpivot(index="sample", variable_name="category", value_name="tags")


def junctions(outdir: Outdir) -> pl.DataFrame:
    """RSeQC junction_annotation per sample: each kind's count and percent, events and junctions."""
    needed = [f"{k}_splicing_{m}{p}" for m in JUNCTION_MEASURES for k in JUNCTION_KINDS
              for p in ("", "_pct")]  # fmt: skip
    table = _held(_columns(_table(outdir, JUNCTIONS), needed, "junction annotation table"),
                  outdir.samples, "RSeQC junction_annotation")  # fmt: skip
    rows = [table.select("sample", measure=pl.lit(m), kind=pl.lit(k),
                         count=pl.col(f"{k}_splicing_{m}"),
                         percent=pl.col(f"{k}_splicing_{m}_pct"))
            for m in JUNCTION_MEASURES for k in JUNCTION_KINDS]  # fmt: skip
    return pl.concat(rows).sort(pl.col("sample").replace_strict(
        {s: i for i, s in enumerate(outdir.samples)}), maintain_order=True)  # fmt: skip


def _curve(outdir: Outdir, plot: str, what: str, names: Mapping[str, str] | None = None,
           drop: Sequence[str] = ()) -> pl.DataFrame:  # fmt: skip
    """A line plot's points per sample: sample, x, y; ``drop`` names MultiQC's own rows."""
    try:
        found = mqc_io.points(outdir.mqc, plot)
    except QcError as error:
        raise ProcessError(str(error)) from error
    kept = found.filter(~pl.col("Sample").is_in(list(drop)))
    named = kept.with_columns(sample=pl.col("Sample").replace(names or {}))
    present = set(named["sample"])
    if missing := [s for s in outdir.samples if s not in present]:
        raise ProcessError(f"no {what} curve for {missing}")
    if extra := sorted(present - set(outdir.samples)):
        raise ProcessError(f"{what} curves of samples not in the sheet: {extra}")
    return named.select("sample", "x", "y")


def coverage(outdir: Outdir) -> pl.DataFrame:
    """Qualimap's normalised coverage along transcripts, 5' at 0 to 3' at 100, per sample."""
    return _curve(outdir, COVERAGE, "Qualimap coverage")


def duplication(outdir: Outdir) -> pl.DataFrame:
    """The duplication (%) dupRadar fits against expression (reads per kilobase), per sample."""
    try:
        found = mqc_io.ragged(outdir.mqc, DUPRADAR)
    except QcError as error:
        raise ProcessError(str(error)) from error
    present = set(found["Sample"])
    if missing := [s for s in outdir.samples if s not in present]:
        raise ProcessError(f"no dupRadar curve for {missing}")
    if extra := sorted(present - set(outdir.samples)):
        raise ProcessError(f"dupRadar curves of samples not in the sheet: {extra}")
    order = {s: i for i, s in enumerate(outdir.samples)}
    return found.rename({"Sample": "sample"}).sort(pl.col("sample").replace_strict(order), "x")


def complexity(outdir: Outdir) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Preseq's expected distinct molecules against total molecules, as MultiQC plots them.

    Returns:
        The curves per sample, and MultiQC's line for a library whose every read is unique.
    """
    names = {f"{s}{PRESEQ_SUFFIX}": s for s in outdir.samples}
    curves = _curve(outdir, PRESEQ, "Preseq", names, drop=(PRESEQ_PERFECT,))
    try:
        found = mqc_io.points(outdir.mqc, PRESEQ)
    except QcError as error:
        raise ProcessError(str(error)) from error
    perfect = found.filter(pl.col("Sample") == PRESEQ_PERFECT).select("x", "y")
    return curves, perfect


def biotypes(outdir: Outdir) -> pl.DataFrame:
    """The reads featureCounts counted per biotype, as nf-core's biotype QC gives MultiQC."""
    table = _table(outdir, BIOTYPES)
    table = _columns(table, table.columns[1:], "biotype table")
    held = _held(table, outdir.samples, "featureCounts' biotypes")
    return held.unpivot(index="sample", variable_name="biotype", value_name="reads")


def matrices(outdir: Outdir) -> pl.DataFrame:
    """Salmon alone and Salmon on STAR's alignment, gene by gene: sample, gene_id, gene_name, a, b.

    ``a`` is Salmon alone, ``b`` Salmon on STAR's alignment, each nf-core's length-scaled count
    as written, not rounded: the two are compared, not fitted.

    Raises:
        ProcessError: If a matrix is absent, or the two differ in samples or genes.
    """
    frames = {}
    for name, relative in MATRICES.items():
        path = outdir.path / relative
        if not path.is_file():
            raise ProcessError(f"{path} is missing")
        header = pl.scan_csv(path, separator="\t", n_rows=0).collect_schema().names()
        if header[:2] != ["gene_id", "gene_name"] or set(header[2:]) != set(outdir.samples):
            raise ProcessError(f"{relative}: its header is not gene_id, gene_name and the "
                               "sheet's samples")  # fmt: skip
        schema = {"gene_id": pl.String, "gene_name": pl.String,
                  **dict.fromkeys(header[2:], pl.Float64)}  # fmt: skip
        frames[name] = pl.read_csv(path, separator="\t", schema=schema)
    alone, star_salmon = frames["salmon"], frames["star_salmon"]
    if set(alone["gene_id"]) != set(star_salmon["gene_id"]):
        raise ProcessError(f"{MATRICES['salmon']} and {MATRICES['star_salmon']} hold other genes")
    long = [f.unpivot(index=["gene_id", "gene_name"], variable_name="sample", value_name=side)
            for side, f in (("a", alone), ("b", star_salmon))]  # fmt: skip
    joined = long[0].join(long[1].drop("gene_name"), on=["gene_id", "sample"],
                          maintain_order="left")  # fmt: skip
    return joined.select("sample", "gene_id", "gene_name", "a", "b")


def steps(outdir: Outdir) -> pl.DataFrame:
    """The trace's cost of each process M3's numbers come from: tasks, samples, peak, time."""
    try:
        frame = trace.read(outdir.path)
    except trace.TraceError as error:
        raise ProcessError(str(error)) from error
    costs = trace.costs(frame.filter(pl.col("process").is_in(list(STEPS))))
    order = {p: i for i, p in enumerate(STEPS)}
    return costs.sort(pl.col("process").replace_strict(order))
