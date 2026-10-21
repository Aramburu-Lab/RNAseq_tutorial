"""M3's alignments to gene counts: featureCounts, the strand call, and the checks I1 to I6 (D-038).

Per-gene counts are data, read with Polars: STAR's ReadsPerGene rows (``align.reads_per_gene``)
and featureCounts' matrix (``fc_genes``), each by its one reader. Every tool log goes through
MultiQC, with one MultiQC run per featureCounts setting: in one run, the settings of one BAM
share a sample name and collapse into one row (M3 review 3). The checks tie the two kinds
together, so a parser that drifts from its tool stops the run instead of drawing a figure.
"""

import re
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import partial
from pathlib import Path

import polars as pl

from ctx import tools
from ctx.process import STRAND, Aligner, ProcessError, Strandedness, cache
from ctx.process.align import RPG_COLUMNS, Alignment, bam_name, reads_per_gene
from ctx.process.reference import sha256
from ctx.qc import multiqc_io as mqc_io
from ctx.qc.checkpoint import gene_types
from ctx.qc.registry import Metric

PLAIN = tuple(f"s{strand}" for _, strand in STRAND.values())
# featureCounts prints a split count to two decimals, so a gene printed with a fraction is within
# 0.005 of the count it summed. A gene printed whole summed a whole number while no read on it
# was split over 7 or more genes: shares of 1/1 to 1/6 sum to multiples of 1/60, and a multiple
# that is not whole lies at least 1/60 from one.
FRACTION_PRINTED = 0.005
FC_META = ("Geneid", "Chr", "Start", "End", "Strand", "Length")
FC_STRANDS = re.compile(r'"-s" "([^"]+)"')
# MultiQC 1.35's STAR plots, one per ReadsPerGene column, and their categories by our row names.
STAR_PLOTS = {
    Strandedness.UNSTRANDED: "star_gene_counts_Unstranded",
    Strandedness.FORWARD: "star_gene_counts_Same_Stranded",
    Strandedness.REVERSE: "star_gene_counts_Reverse_Stranded",
}
GENES = "genes"
STAR_CATEGORIES = {"N_unmapped": "Unmapped", "N_multimapping": "Multimapping",
                   "N_noFeature": "No Feature", "N_ambiguous": "Ambiguous Features",
                   GENES: "Overlapping Genes"}  # fmt: skip
FC_COLUMNS = {"Assigned": "assigned", "Unassigned_NoFeatures": "no_feature",
              "Unassigned_Ambiguity": "ambiguous"}  # fmt: skip
SHOWN = 5


@dataclass(frozen=True)
class FcOptions:
    """One featureCounts setting: a -s per BAM, in BAM order, and the counting switches."""

    strands: tuple[int, ...]
    paired: bool = False
    overlap: bool = False
    fraction: bool = False
    multimap: bool = False

    def args(self) -> list[str]:
        """The setting as featureCounts options; a pair is counted once (--countReadPairs).

        Since subread 2.0.2, -p alone only declares pairs and counts each read: on the paired
        fixture, CCDC157's 506 pairs gave 1,012.
        """
        args = ["-s", ",".join(map(str, self.strands))]
        args += ["-p", "--countReadPairs"] if self.paired else []
        args += ["-O"] if self.overlap else []
        args += ["-M"] if self.multimap else []
        return args + (["--fraction"] if self.fraction else [])

    @property
    def sums_to_assigned(self) -> bool:
        """Whether the matrix sums to Assigned: not when -O or -M counts a read more than once.

        -O with --fraction splits each read over its genes, so it sums to Assigned again.
        """
        return not self.multimap and (self.fraction or not self.overlap)


@dataclass(frozen=True)
class Counts:
    """A keyed count folder: one featureCounts matrix and summary per setting."""

    folder: Path
    built: bool
    settings: Mapping[str, FcOptions]

    def matrix(self, setting: str) -> Path:
        """The setting's gene by BAM matrix."""
        return self.folder / f"{setting}.featureCounts.txt"

    def summary(self, setting: str) -> Path:
        """The setting's summary, which MultiQC reads."""
        return self.folder / f"{setting}.featureCounts.txt.summary"


def plain_settings(alignment: Alignment) -> dict[str, FcOptions]:
    """Settings -s 0, 1 and 2 on every BAM: what the strand call is made from."""
    n = len(alignment.samples)
    return {f"s{s}": FcOptions((s,) * n, alignment.paired) for _, s in STRAND.values()}


def variant_settings(alignment: Alignment, strands: Mapping[str, Strandedness]
                     ) -> dict[str, FcOptions]:  # fmt: skip
    """-O, -O --fraction and -M, each BAM at its own called strand."""
    called = tuple(STRAND[strands[s]][1] for s in alignment.samples)
    return {
        "O": FcOptions(called, alignment.paired, overlap=True),
        "O_fraction": FcOptions(called, alignment.paired, overlap=True, fraction=True),
        "M": FcOptions(called, alignment.paired, multimap=True),
    }  # fmt: skip


def feature_counts(alignment: Alignment, samples: Sequence[str], gtf: Path, out: Path,
                   options: FcOptions, threads: int) -> None:  # fmt: skip
    """Count the samples' BAMs at once into ``out``, in the order given.

    featureCounts runs in the alignment folder and is given each BAM by its name there, which
    heads the BAM's column, so a matrix names no folder and a moved project reuses it. MultiQC
    names each sample from that name (``<sample>.bam`` gives ``<sample>``).
    """
    out = out.resolve()
    tools.run("featureCounts", ["-T", str(threads), "-a", str(gtf.resolve()), *options.args(),
                                "-o", str(out), *(bam_name(s) for s in samples)],
              cwd=alignment.folder.resolve(), log=out.with_suffix(".log"))  # fmt: skip


def count(alignment: Alignment, gtf: Path, parent: Path, settings: Mapping[str, FcOptions],
          threads: int) -> Counts:  # fmt: skip
    """Every setting over every BAM of an alignment, in one keyed folder.

    Args:
        alignment: The BAMs, counted in its sample order.
        gtf: The annotation the reference was built with.
        parent: Where keyed folders sit.
        settings: Each setting's name and options (``plain_settings``, ``variant_settings``).
        threads: featureCounts' threads.

    Returns:
        The folder and its settings.
    """
    # The BAMs in order: the matrix columns and each -s follow it, so a reordered sheet
    # builds a new folder (I5).
    inputs = {"tool": "featureCounts", "version": tools.TOOLS["featureCounts"].version,
              "alignment": alignment.folder.name, "gtf_sha256": sha256(gtf),
              "bams": [bam_name(s) for s in alignment.samples],
              "settings": {name: asdict(o) for name, o in settings.items()}}  # fmt: skip

    def build(dest: Path) -> None:
        """Run featureCounts once per setting."""
        for name, options in settings.items():
            feature_counts(alignment, alignment.samples, gtf, dest / f"{name}.featureCounts.txt",
                           options, threads)  # fmt: skip

    folder, built = cache.keyed(parent, "count", inputs, build)
    return Counts(folder, built, dict(settings))


def star_genes(tab: Path) -> pl.DataFrame:
    """STAR's gene rows of ReadsPerGene.out.tab, its four N_ rows left out."""
    return reads_per_gene(tab).filter(~pl.col("gene_id").str.starts_with("N_"))


def fc_genes(matrix: Path, *, fraction: bool) -> pl.DataFrame:
    """The featureCounts matrix: gene_id, then one count column per BAM, headed by its name.

    Args:
        matrix: The ``.featureCounts.txt`` file.
        fraction: Whether counts were split (--fraction), so read as floats; else integers.

    Raises:
        ProcessError: If a count is not a number, or not whole where counts are integers.
    """
    table = pl.read_csv(matrix, separator="\t", comment_prefix="#", infer_schema=False)
    kind = pl.Float64 if fraction else pl.Int64
    bams = [c for c in table.columns if c not in FC_META]
    try:
        return table.select(pl.col("Geneid").alias("gene_id"), pl.col(bams).cast(kind, strict=True))
    except pl.exceptions.InvalidOperationError as error:
        raise ProcessError(f"{matrix} holds a count that is not {kind}: {error}") from error


def summaries(counts: Counts, alignment: Alignment, dest: Path, config: Path,
              jobs: int = 1) -> dict[str, mqc_io.MultiqcRun]:  # fmt: skip
    """One MultiQC run per setting into ``dest/fc_<setting>/``, each with one row per BAM.

    ``jobs`` of them run at a time (``mqc_io.at_once``).

    Raises:
        ProcessError: If a report's rows are not exactly the alignment's samples.
    """

    def one(setting: str) -> mqc_io.MultiqcRun:
        run = mqc_io.run([counts.summary(setting)], dest / f"fc_{setting}", config=config,
                         modules=("featurecounts",))  # fmt: skip
        names = mqc_io.table(run, "featurecounts")["Sample"].to_list()
        if sorted(names) != sorted(alignment.samples):
            raise ProcessError(f"MultiQC's featureCounts rows for {setting} are {names}, not "
                               f"one per BAM of {list(alignment.samples)}")  # fmt: skip
        return run

    return mqc_io.at_once({s: partial(one, s) for s in counts.settings}, jobs)


def assigned(runs: Mapping[str, mqc_io.MultiqcRun]) -> pl.DataFrame:
    """Assigned, no-feature and ambiguous reads per sample and setting, from MultiQC.

    Raises:
        ProcessError: If a report lacks one of those columns.
    """
    frames = []
    for setting, run in runs.items():
        table = mqc_io.table(run, "featurecounts")
        if missing := [c for c in FC_COLUMNS if c not in table.columns]:
            raise ProcessError(f"MultiQC's featureCounts table for {setting} has no {missing}")
        frames.append(mqc_io.numeric(table, FC_COLUMNS).select(
            sample=pl.col("Sample"), setting=pl.lit(setting), **{v: pl.col(k) for k, v in
                                                                 FC_COLUMNS.items()}))  # fmt: skip
    return pl.concat(frames)


def strand_calls(counted: pl.DataFrame, metric: Metric) -> pl.DataFrame:
    """Each sample's forward share, Assigned at -s 1 over -s 1 plus -s 2, and its call.

    Args:
        counted: ``assigned``'s table, holding settings s1 and s2.
        metric: The registry's ``gene_strand_bias``, whose bands make the call (D-039).

    Returns:
        sample, assigned_s1, assigned_s2, forward_share (null when neither assigned a read)
        and call: forward, reverse, unstranded or undetermined. The call is made on the two
        counts exactly (``Metric.strandedness``), so (9, 11) and (11, 9) get mirrored calls.
    """
    wide = counted.filter(pl.col("setting").is_in(["s1", "s2"])).pivot(
        "setting", index="sample", values="assigned", maintain_order=True
    )  # fmt: skip
    shares = wide.select(
        "sample", assigned_s1=pl.col("s1"), assigned_s2=pl.col("s2"),
        forward_share=pl.when(pl.col("s1") + pl.col("s2") > 0)
        .then(pl.col("s1") / (pl.col("s1") + pl.col("s2"))),
    )  # fmt: skip
    counts = shares.select("assigned_s1", "assigned_s2", "forward_share").rows()
    calls = [metric.strandedness(Fraction(s1) / (Fraction(s1) + Fraction(s2)))
             if share is not None else "undetermined" for s1, s2, share in counts]  # fmt: skip
    return shares.with_columns(call=pl.Series(calls, dtype=pl.String))


def called_strands(calls: pl.DataFrame, given: Strandedness) -> dict[str, Strandedness]:
    """The strand each sample is counted at: its call under auto, else the one given.

    Raises:
        ProcessError: If a sample had no read assigned to a gene at -s 1 or -s 2, whatever is
            given: it has no share to call and nothing to count. Under auto, if any call is
            undetermined; the samples and shares listed.
    """
    if (empty := calls.filter(pl.col("forward_share").is_null())).height:
        raise ProcessError(
            f"no read of {empty['sample'].to_list()} was assigned to a gene at -s 1 or -s 2, "
            "so no strand can be called; check that the sample aligned and that the "
            "annotation is the reference's"
        )
    if given is not Strandedness.AUTO:
        return dict.fromkeys(calls["sample"], given)
    if (open_ := calls.filter(pl.col("call") == "undetermined")).height:
        shares = dict(open_.select("sample", "forward_share").rows())
        raise ProcessError(
            f"the strand of {open_.height} sample(s) is undetermined, a forward share between "
            f"the registry's gene_strand_bias bands: {shares}; pass --strandedness "
            "unstranded, forward or reverse"
        )
    return {s: Strandedness(c) for s, c in calls.select("sample", "call").rows()}


def expected_input(fastp: pl.DataFrame, files: pl.DataFrame, column: str, *,
                   paired: bool) -> dict[str, int]:  # fmt: skip
    """The reads each sample's aligner should read: its runs' fastp counts summed.

    Args:
        fastp: ``align.fastp_reads``' table.
        files: ``reads.fastq_files``' table, which runs make which sample.
        column: passed_reads for trimmed reads, raw_reads for the untrimmed setting.
        paired: Whether runs are pairs: fastp counts both mates, an aligner counts the pair once.
    """
    runs = files.unique("run_accession", keep="first", maintain_order=True)
    joined = runs.join(fastp, on="run_accession", how="left")
    totals = joined.group_by("sample", maintain_order=True).agg(pl.col(column).sum())
    return {s: int(n) // (2 if paired else 1) for s, n in totals.rows()}


def input_problems(aligned: pl.DataFrame, expected: Mapping[str, int]) -> list[str]:
    """I1: each sample's aligner read exactly the reads fastp passed (or the capped raw ones).

    The guard D-036 needs: a STAR that reads nothing passes ``--version`` and exits 0.
    """
    said = dict(aligned.select("sample", "input_reads").rows())
    return [f"I1: {s}'s aligner read {said.get(s)} reads of the {n} it was given"
            for s, n in expected.items() if said.get(s) != n]  # fmt: skip


def _star_categories(tab: Path) -> dict[str, int]:
    """A ReadsPerGene file's totals by column and category, keyed ``<column>|<category>``."""
    rows = reads_per_gene(tab)
    genes = rows.filter(~pl.col("gene_id").str.starts_with("N_"))
    out = {f"{c}|{g}": int(n) for g, *ns in rows.filter(pl.col("gene_id").str.starts_with("N_"))
           .rows() for c, n in zip(RPG_COLUMNS[1:], ns, strict=True)}  # fmt: skip
    return out | {f"{c}|{GENES}": int(genes[c].sum()) for c in RPG_COLUMNS[1:]}


def star_parser_problems(alignment: Alignment, run: mqc_io.MultiqcRun,
                         aligned: pl.DataFrame) -> list[str]:  # fmt: skip
    """I2: our ReadsPerGene totals equal MultiQC's, category by category, and sum to the input.

    MultiQC 1.35 leaves out a category that is zero in every sample, so an absent column must
    be zero in every sample of ours too, never assumed. Two parsers agreeing, not biology.
    """
    problems, said = [], dict(aligned.select("sample", "input_reads").rows())
    ours = {s: _star_categories(alignment.folder / f"{s}.ReadsPerGene.out.tab")
            for s in alignment.samples}  # fmt: skip
    for strand, plot in STAR_PLOTS.items():
        table = mqc_io.plot_table(run, plot)
        theirs = {row["Sample"]: row for row in table.iter_rows(named=True)}
        for sample in alignment.samples:
            if sample not in theirs:
                problems.append(f"I2: MultiQC's {plot} has no row for {sample}")
                continue
            mine = {k: ours[sample][f"{strand}|{k}"] for k in STAR_CATEGORIES}
            got = {k: theirs[sample].get(v, 0.0 if v not in table.columns else None)
                   for k, v in STAR_CATEGORIES.items()}  # fmt: skip
            if mine != got:
                problems.append(f"I2: {sample} {strand}: ReadsPerGene {mine}, MultiQC {got}")
            if sum(mine.values()) != said.get(sample):
                problems.append(f"I2: {sample} {strand}: categories sum to "
                                f"{sum(mine.values())}, STAR read {said.get(sample)}")  # fmt: skip
    return problems


def assigned_problems(counts: Counts, counted: pl.DataFrame, alignment: Alignment) -> list[str]:
    """I3: each matrix column sums to MultiQC's Assigned, where the setting counts reads once.

    A split count is printed to two decimals, so the sum may differ by 0.005 for each gene
    printed with a fraction (``FRACTION_PRINTED``); whole counts must match exactly.
    """
    problems = []
    for setting, options in counts.settings.items():
        if not options.sums_to_assigned:
            continue
        matrix = fc_genes(counts.matrix(setting), fraction=options.fraction)
        said = dict(counted.filter(pl.col("setting") == setting)
                    .select("sample", "assigned").rows())  # fmt: skip
        for sample in alignment.samples:
            if (column := bam_name(sample)) not in matrix.columns or sample not in said:
                problems.append(f"I3: {sample} {setting}: no matrix column {column} or no "
                                f"Assigned; columns {matrix.columns[1:]}, Assigned for "
                                f"{sorted(said)}")  # fmt: skip
                continue
            genes = matrix[column]
            slack = FRACTION_PRINTED * int((genes != genes.round(0)).sum())
            if abs(float(genes.sum()) - said[sample]) > slack:
                problems.append(f"I3: {sample} {setting}: matrix sums to {genes.sum()}, "
                                f"MultiQC's Assigned is {said[sample]}")  # fmt: skip
    return problems


def gene_set_problems(gtf: Path, tables: Mapping[str, pl.DataFrame]) -> list[str]:
    """I4: every count table holds exactly the annotation's genes, by gene_id."""
    genes = set(gene_types(gtf)["gene_id"])
    problems = []
    for name, table in tables.items():
        found = set(table["gene_id"])
        if found != genes:
            problems.append(f"I4: {name} has {len(found)} genes, the annotation {len(genes)}; "
                            f"missing {sorted(genes - found)[:SHOWN]}, extra "
                            f"{sorted(found - genes)[:SHOWN]}")  # fmt: skip
    return problems


def order_problems(counts: Counts, alignment: Alignment) -> list[str]:
    """I5: each matrix's columns are the BAMs in sample order, at the -s each was given.

    featureCounts records its command line in the matrix's first line; its -s list must be
    the setting's, value by value, so a strand cannot land on another sample's BAM.
    """
    bams = [bam_name(s) for s in alignment.samples]
    problems = []
    for setting, options in counts.settings.items():
        with counts.matrix(setting).open() as handle:
            command, header = handle.readline(), handle.readline().rstrip("\n").split("\t")
        header = header[len(FC_META) :]
        if header != bams:
            problems.append(f"I5: {setting}'s columns are {header}, not {bams}")
        found = FC_STRANDS.search(command)
        given = found[1] if found else None
        if given != ",".join(map(str, options.strands)):
            problems.append(f"I5: {setting} ran with -s {given}, not {options.strands}")
    return problems


def primary_problems(run: mqc_io.MultiqcRun, aligned: pl.DataFrame,
                     alignment: Alignment) -> list[str]:  # fmt: skip
    """I6: flagstat's primary records equal the aligner's reads, in all and mapped.

    Primary = records - secondary - supplementary, from MultiQC's flagstat table. Without
    unmapped reads in the BAM (STAR's ``Within``, HISAT2 without ``--no-unal``), flagstat's
    total falls short of the input. A pair is two records, and mapped records are the
    aligner's mapped mates (``align.aligner_reads``). STAR counts a pair mapped when one mate
    maps, so under STAR flagstat's singletons, mapped mates of an unmapped mate, add that mate
    back: on the paired fixture with 300 mates made random 20-mers, 4,974 mapped records and
    20 singletons against STAR's 2,497 mapped pairs (measured, linux STAR 2.7.11b).
    """
    halves = alignment.paired and alignment.aligner is Aligner.STAR
    columns = ["total_passed", "secondary_passed", "supplementary_passed", "mapped_passed"]
    columns += ["singletons_passed"] if halves else []
    table = mqc_io.table(run, "samtools_flagstat")
    if missing := [c for c in columns if c not in table.columns]:
        raise ProcessError(f"MultiQC's flagstat table has no {missing}")
    flags = {r["Sample"]: r for r in mqc_io.numeric(table, columns).iter_rows(named=True)}
    mates, problems = 2 if alignment.paired else 1, []
    for sample, reads, mapped in aligned.select("sample", "input_reads", "mapped_mates").rows():
        row = flags.get(sample)
        if row is None:
            problems.append(f"I6: no flagstat for {sample}")
            continue
        extra = row["secondary_passed"] + row["supplementary_passed"]
        back = row["singletons_passed"] if halves else 0
        primary = (row["total_passed"] - extra, row["mapped_passed"] - extra + back)
        if primary != (mates * reads, mapped):
            problems.append(f"I6: {sample}: flagstat's primary records {primary} (all, mapped), "
                            f"the aligner's {(mates * reads, mapped)}")  # fmt: skip
    return problems


def check(*problems: list[str]) -> None:
    """Stop on any failed invariant, every failure listed.

    Raises:
        ProcessError: If any list holds a problem.
    """
    if found := [p for group in problems for p in group]:
        raise ProcessError("checks failed (D-038):\n  " + "\n  ".join(found))
