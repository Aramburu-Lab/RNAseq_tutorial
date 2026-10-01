"""ctx.process.count on the chr22 fixtures: featureCounts, the strand call, the checks I1 to I6.

Tests marked needs_star skip only when the probe fails on darwin-arm64 with CI unset, and fail
everywhere else (R6.2).
"""

import random
import shutil
import tomllib
from pathlib import Path

import polars as pl
import pytest

from ctx import concordance
from ctx.process import Aligner, ProcessError, Strandedness, align, count, reference
from ctx.qc import Fastp, multiqc_io, registry
from ctx.qc.checkpoint import gene_types
from ctx.qc.reads import raw_config

ROOT = Path(__file__).resolve().parents[1]
SLICE = ROOT / "tests" / "fixtures" / "chr22" / "slice.fa"
GTF = ROOT / "tests" / "fixtures" / "chr22" / "slice.gtf"
PAIRED = ROOT / "tests" / "fixtures" / "chr22"
SE = ROOT / "tests" / "fixtures" / "chr22_se"
SAMPLES = ("se_a", "se_b", "se_c", "se_d", "se_e")
REGISTRY = ROOT / "conf" / "qc_metrics.toml"
NAMES = dict(gene_types(GTF).select("gene_name", "gene_id").rows())


def frame(rows: list[tuple[str, str, list[Path]]]) -> pl.DataFrame:
    """A ``reads.fastq_files`` table: (run, sample, files read 1 first) per run."""
    return pl.DataFrame(
        [(run, sample, mate, str(p)) for run, sample, paths in rows
         for mate, p in enumerate(paths, 1)],
        schema={"run_accession": pl.String, "sample": pl.String, "mate": pl.Int32,
                "path": pl.String}, orient="row",
    )  # fmt: skip


class Run:
    """The five single-end samples through trim, HISAT2, both count folders and MultiQC."""

    def __init__(self, parent: Path) -> None:
        self.ref = reference.prepare(SLICE, GTF, parent)
        index = reference.hisat2_index(self.ref, parent, 2)[0]
        self.files = frame([(s, s, [SE / f"{s}.fq"]) for s in SAMPLES])
        self.trimmed = align.trim(self.files, parent, fastp=Fastp(), cap=None, threads=2)
        reads = align.by_sample(self.files, {r: self.trimmed.reads(r) for r in SAMPLES})
        self.alignment = align.align(reads, parent, align.AlignOptions(Aligner.HISAT2, index, 2))
        self.config = raw_config(parent / "config.json")
        self.main = align.summarise(self.alignment, parent / "main", self.config, self.trimmed)
        self.plain = count.count(self.alignment, self.ref.gtf, parent,
                                 count.plain_settings(self.alignment), 2)  # fmt: skip
        self.plain_runs = count.summaries(self.plain, self.alignment, parent, self.config)
        counted = count.assigned(self.plain_runs)
        metric = registry.load(REGISTRY).metrics["gene_strand_bias"]
        self.calls = count.strand_calls(counted, metric)
        strands = count.called_strands(self.calls, Strandedness.AUTO)
        self.variants = count.count(self.alignment, self.ref.gtf, parent,
                                    count.variant_settings(self.alignment, strands), 2)  # fmt: skip
        self.variant_runs = count.summaries(self.variants, self.alignment, parent, self.config)
        self.counted = pl.concat([counted, count.assigned(self.variant_runs)])


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    return Run(tmp_path_factory.mktemp("count"))


def test_the_readers_on_hand_written_files(tmp_path):
    """ReadsPerGene's N_ rows are not genes, and its columns 3 and 4 are forward and reverse;
    the featureCounts matrix gives one column per BAM, integers unless split.

    Mutations caught: N_ rows read as genes; STAR's columns 3 and 4 swapped in STRAND."""
    tab = tmp_path / "s.ReadsPerGene.out.tab"
    tab.write_text("N_unmapped\t3\t3\t3\nN_multimapping\t1\t1\t1\nN_noFeature\t4\t950\t2\n"
                   "N_ambiguous\t437\t0\t388\nG1\t500\t40\t600\nG2\t55\t6\t7\n")  # fmt: skip
    genes = count.star_genes(tab)
    assert genes.columns == ["gene_id", "unstranded", "forward", "reverse"]
    assert genes.rows() == [("G1", 500, 40, 600), ("G2", 55, 6, 7)]
    matrix = tmp_path / "s0.featureCounts.txt"
    matrix.write_text('# Program:featureCounts v2.1.1; Command:"featureCounts" "-s" "2,0"\n'
                      "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/x/a.bam\t/x/b.bam\n"
                      "G1\tc\t1\t9\t+\t9\t5\t0\nG2\tc;c\t3;7\t5;9\t-;-\t6\t1\t2\n")  # fmt: skip
    table = count.fc_genes(matrix, fraction=False)
    assert table.columns == ["gene_id", "/x/a.bam", "/x/b.bam"]
    assert table.rows() == [("G1", 5, 0), ("G2", 1, 2)]
    assert table.schema["/x/a.bam"] == pl.Int64
    assert count.fc_genes(matrix, fraction=True).schema["/x/a.bam"] == pl.Float64


def test_a_split_count_in_an_integer_matrix_is_refused(tmp_path):
    """Mutation caught: the cast not strict (0.50 read as a null or a 0)."""
    matrix = tmp_path / "s0.featureCounts.txt"
    matrix.write_text("Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\nG1\tc\t1\t9\t+\t9\t0.50\n")
    with pytest.raises(ProcessError, match="not Int64"):
        count.fc_genes(matrix, fraction=False)
    assert count.fc_genes(matrix, fraction=True)["a.bam"].to_list() == [0.5]


def test_the_strand_is_called_from_the_forward_share(run):
    """se_a to se_d are reverse-stranded and se_e unstranded (samples.tsv), from featureCounts'
    Assigned at -s 1 over -s 1 plus -s 2: se_a 49 / (49 + 612), se_e 368 / (368 + 329).

    Mutation caught: the share inverted, -s 2 over the sum (se_a turns forward)."""
    calls = dict(run.calls.select("sample", "call").rows())
    assert calls == {"se_a": "reverse", "se_b": "reverse", "se_c": "reverse", "se_d": "reverse",
                     "se_e": "unstranded"}  # fmt: skip
    se_a = run.calls.filter(pl.col("sample") == "se_a").row(0, named=True)
    assert (se_a["assigned_s1"], se_a["assigned_s2"]) == (49.0, 612.0)
    assert se_a["forward_share"] == pytest.approx(49 / 661)
    se_e = run.calls.filter(pl.col("sample") == "se_e")["forward_share"].item()
    assert se_e == pytest.approx(368 / 697)


def test_an_undetermined_strand_stops_the_run(run, tmp_path):
    """The registry's bands decide, never a threshold in code: stranded at 0.95 leaves se_a
    (reverse share 1 - 0.0741 = 0.926) undetermined, unstranded at 0.01 leaves se_e (share
    0.528) undetermined, and the run stops naming both; a given strand overrides.

    Mutations caught: a hard-coded 0.8 in place of the registry's stranded band (se_a passes);
    a hard-coded 0.1 in place of its unstranded band (se_e passes)."""
    text = REGISTRY.read_text()
    narrowed = tmp_path / "qc_metrics.toml"
    narrowed.write_text(text.replace("stranded = 0.8\nunstranded = 0.1",
                                     "stranded = 0.95\nunstranded = 0.01"))  # fmt: skip
    metric = registry.load(narrowed).metrics["gene_strand_bias"]
    bands = tomllib.loads(narrowed.read_text())["metrics"]["gene_strand_bias"]
    assert (bands["stranded"], bands["unstranded"]) == (0.95, 0.01)
    calls = count.strand_calls(count.assigned(run.plain_runs), metric)
    with pytest.raises(ProcessError, match=r"undetermined.*'se_a': 0\.07.*'se_e': 0\.52"):
        count.called_strands(calls, Strandedness.AUTO)
    assert set(count.called_strands(calls, Strandedness.REVERSE).values()) == {
        Strandedness.REVERSE}  # fmt: skip


def synthetic(counts: dict[str, tuple[float, float]]) -> pl.DataFrame:
    """``assigned``'s table for samples with these Assigned reads at -s 1 and -s 2."""
    rows = [
        (s, f"s{n}", v) for s, pair in counts.items() for n, v in zip((1, 2), pair, strict=True)
    ]
    return pl.DataFrame(rows, schema=["sample", "setting", "assigned"], orient="row")


def test_a_share_and_its_mirror_get_mirrored_calls():
    """9 and 11 reads, and 11 and 9, at -s 1 and -s 2: both undetermined, on the band's edge
    (|2 x 0.45 - 1| = 0.1, not below it), as are 45,000 and 55,000 both ways.

    Mutation caught: the registry deciding in floats (9 / 20 unstranded, 11 / 20
    undetermined: 2 x 0.45 - 1 is -0.09999999999999998 in floats). Passing the float share
    here instead of the two counts is equivalent, since the registry reads a float as it
    prints."""
    metric = registry.load(REGISTRY).metrics["gene_strand_bias"]
    pairs = {"a": (9, 11), "b": (11, 9), "c": (45000, 55000), "d": (55000, 45000), "e": (12, 8)}
    calls = dict(count.strand_calls(synthetic(pairs), metric).select("sample", "call").rows())
    assert calls == {"a": "undetermined", "b": "undetermined", "c": "undetermined",
                     "d": "undetermined", "e": "undetermined"}  # fmt: skip
    assert count.strand_calls(synthetic({"f": (11, 10), "g": (10, 11)}), metric)[
        "call"].to_list() == ["unstranded", "unstranded"]  # fmt: skip


def test_a_sample_with_no_read_on_a_gene_stops_the_run_whatever_is_given():
    """No read assigned at -s 1 or -s 2: a refusal of its own, which --strandedness cannot
    override, never "a share between the bands".

    Mutation caught: the empty sample handled as undetermined (a given strand lets it pass)."""
    metric = registry.load(REGISTRY).metrics["gene_strand_bias"]
    calls = count.strand_calls(synthetic({"ok": (49, 612), "empty": (0, 0)}), metric)
    for given in (Strandedness.AUTO, Strandedness.REVERSE):
        with pytest.raises(ProcessError, match=r"no read of \['empty'\] was assigned to a gene"):
            count.called_strands(calls, given)


def test_variants_count_at_each_samples_called_strand(run):
    """se_e is counted at -s 0, the reverse-stranded samples at -s 2.

    Mutation caught: one strand for every BAM (se_e counted at -s 2)."""
    assert {o.strands for o in run.variants.settings.values()} == {(2, 2, 2, 2, 0)}


def test_one_multiqc_run_per_setting(run):
    """Each setting's report has one row per BAM, with that setting's own numbers.

    Mutation caught: one MultiQC run over every setting (one BAM's six summaries collapse
    into one row, so -s 1 and -s 2 report the same Assigned)."""
    folders = {r.data.parent.name for r in {**run.plain_runs, **run.variant_runs}.values()}
    assert folders == {f"fc_{s}" for s in ("s0", "s1", "s2", "O", "O_fraction", "M")}
    for setting, mqc in {**run.plain_runs, **run.variant_runs}.items():
        rows = multiqc_io.table(mqc, "featurecounts")["Sample"].to_list()
        assert sorted(rows) == list(SAMPLES), setting
    se_a = dict(run.counted.filter(pl.col("sample") == "se_a").select("setting", "assigned").rows())
    assert (se_a["s0"], se_a["s1"], se_a["s2"]) == (563.0, 49.0, 612.0)


def test_fraction_sums_to_assigned_and_overlap_exceeds_it(run):
    """-O counts a read once per gene it overlaps; --fraction splits it, so the sum returns to
    Assigned: se_a's -O matrix sums to 1,390 of 1,000 assigned, -O --fraction to 1,000.01,
    within featureCounts' two printed decimals.

    Mutations caught: --fraction dropped (the -O --fraction sum is 1,390, and I3 fails); no
    slack for printed fractions (I3 fails on 1,000.01); -O alone held to Assigned (I3 fails)."""
    bam = align.bam_name("se_a")
    overlap = count.fc_genes(run.variants.matrix("O"), fraction=False)[bam].sum()
    split = count.fc_genes(run.variants.matrix("O_fraction"), fraction=True)[bam].sum()
    assigned = dict(run.counted.filter(pl.col("sample") == "se_a")
                    .select("setting", "assigned").rows())  # fmt: skip
    assert (overlap, assigned["O"]) == (1390, 1000.0)
    assert split == pytest.approx(1000.01) and assigned["O_fraction"] == 1000.0
    assert count.assigned_problems(run.variants, run.counted, run.alignment) == []


def test_every_invariant_holds_on_the_fixture(run):
    """I1, I3, I4, I5 and I6 on the five HISAT2 alignments (I2 is STAR's own, below).

    Mutation caught: featureCounts given the -s list in reverse (I5 reads the -s it ran)."""
    aligned = align.aligner_reads(run.main, run.alignment)
    fastp = align.fastp_reads(run.main, SAMPLES)
    tables = {f"{c.folder.name}/{s}": count.fc_genes(c.matrix(s), fraction=o.fraction)
              for c in (run.plain, run.variants) for s, o in c.settings.items()}  # fmt: skip
    count.check(
        count.input_problems(aligned, count.expected_input(fastp, run.files, "passed_reads",
                                                           paired=False)),
        count.assigned_problems(run.plain, run.counted, run.alignment),
        count.assigned_problems(run.variants, run.counted, run.alignment),
        count.gene_set_problems(run.ref.gtf, tables),
        count.order_problems(run.plain, run.alignment),
        count.order_problems(run.variants, run.alignment),
        count.primary_problems(run.main, aligned, run.alignment),
    )  # fmt: skip


def test_a_shuffled_bam_order_is_caught(run, tmp_path):
    """I5: BAMs handed to featureCounts in the reverse of the samples' -s list, the called
    strands (2, 2, 2, 2, 0), so se_a's BAM is counted at se_e's -s 0.

    Mutation caught: I5 without its header check (the -s list matches; the swap is unseen)."""
    strands = {**dict.fromkeys(SAMPLES, Strandedness.REVERSE), "se_e": Strandedness.UNSTRANDED}
    options = count.variant_settings(run.alignment, strands)["M"]
    assert options.strands == (2, 2, 2, 2, 0)
    count.feature_counts(run.alignment, list(reversed(SAMPLES)), run.ref.gtf,
                         tmp_path / "M.featureCounts.txt", options, 1)  # fmt: skip
    counts = count.Counts(tmp_path, True, {"M": options})
    problems = count.order_problems(counts, run.alignment)
    assert len(problems) == 1 and problems[0].startswith("I5: M's columns are")


def test_a_foreign_gene_and_a_wrong_input_are_caught(run, tmp_path):
    """I4 against an annotation whose MTFP1 has another id, so both hold 10 genes; I1 against
    one read more for se_b and one fewer for se_c, so the totals over samples agree.

    Mutations caught: I4 comparing gene counts only; I1 comparing totals over samples."""
    renamed = tmp_path / "renamed.gtf"
    renamed.write_text(run.ref.gtf.read_text().replace(NAMES["MTFP1"], "ENSG00000242114.99"))
    table = count.fc_genes(run.plain.matrix("s0"), fraction=False)
    problems = count.gene_set_problems(renamed, {"s0": table})
    assert len(problems) == 1 and "has 10 genes, the annotation 10" in problems[0]
    assert NAMES["MTFP1"] in problems[0] and "ENSG00000242114.99" in problems[0]
    aligned = align.aligner_reads(run.main, run.alignment)
    off = {s: 1000 + (s == "se_b") - (s == "se_c") for s in SAMPLES}
    assert count.input_problems(aligned, off) == [
        "I1: se_b's aligner read 1000 reads of the 1001 it was given",
        "I1: se_c's aligner read 1000 reads of the 999 it was given",
    ]  # fmt: skip


def fraction_counts(folder: Path, column: list[str], assigned: float):
    """A one-sample -O --fraction count folder holding ``column`` as its matrix, and MultiQC's
    Assigned for it: (counts, counted, alignment) for ``assigned_problems``."""
    rows = "".join(f"G{n}\tc\t1\t9\t+\t9\t{v}\n" for n, v in enumerate(column))
    (folder / "O_fraction.featureCounts.txt").write_text(
        "# featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tx.bam\n" + rows
    )  # fmt: skip
    options = count.FcOptions((2,), overlap=True, fraction=True)
    counts = count.Counts(folder, True, {"O_fraction": options})
    counted = pl.DataFrame({"sample": ["x"], "setting": ["O_fraction"], "assigned": [assigned]})
    return counts, counted, align.Alignment(folder, True, Aligner.HISAT2, ("x",), False)


@pytest.mark.parametrize(
    ("column", "assigned", "holds"),
    [(["0.33", "0.33", "0.33", "1.00"], 2.0, True),
     (["0.33", "0.33", "0.33", "1.00"], 2.0 + 3 * 0.005 + 0.002, False),
     (["992.00", *["0"] * 1746], 1000.0, False),
     (["992.00", *["0"] * 1746], 992.0, True)],
)  # fmt: skip
def test_i3_slack_covers_printed_fractions_only(tmp_path, column, assigned, holds):
    """0.005 per gene printed with a fraction: three thirds printed 0.33 pass against 2 reads,
    and fail beyond 3 x 0.005; a 1,747-gene matrix with one whole gene 8 reads short fails.

    Mutations caught: a whole read of slack per gene; slack of 0.005 over every gene of the
    matrix (8.7 reads at 1,747 genes, the shortfall passes)."""
    problems = count.assigned_problems(*fraction_counts(tmp_path, column, assigned))
    assert (problems == []) is holds, problems


def test_a_missing_matrix_column_is_a_problem_not_a_crash(tmp_path):
    """Mutation caught: the column read without its check (a Polars ColumnNotFoundError)."""
    counts, counted, _ = fraction_counts(tmp_path, ["1.00"], 1.0)
    other = align.Alignment(tmp_path, True, Aligner.HISAT2, ("y",), False)
    problems = count.assigned_problems(counts, counted, other)
    assert len(problems) == 1 and "no matrix column y.bam" in problems[0]


def test_a_moved_root_and_a_reordered_sheet_rerun_clean(tmp_path, run):
    """Two samples built under A; A moved to B and rerun: every folder reused and I3 and I5
    hold. The same root with the samples reversed: a new count folder, and both hold again.

    Mutations caught: BAM paths in the matrix (after the move I5 names A's paths and I3
    finds no column); the key without the BAM order (the reversed sheet reuses the old
    matrix and I5 fails)."""
    index = reference.hisat2_index(run.ref, tmp_path / "index", 2)[0]
    settings = {"s2": count.FcOptions((2, 2))}

    def once(root: Path, samples: tuple[str, ...]) -> tuple[bool, bool, list[str]]:
        alignment = align.align({s: [(SE / f"{s}.fq",)] for s in samples}, root / "tools",
                                align.AlignOptions(Aligner.HISAT2, index, 2))  # fmt: skip
        counts = count.count(alignment, run.ref.gtf, root / "tools", settings, 2)
        mqc = count.summaries(counts, alignment, root / "multiqc" / "_".join(samples), run.config)
        problems = count.order_problems(counts, alignment)
        problems += count.assigned_problems(counts, count.assigned(mqc), alignment)
        return alignment.built, counts.built, problems

    assert once(tmp_path / "A", ("se_a", "se_e")) == (True, True, [])
    shutil.move(tmp_path / "A", tmp_path / "B")
    assert once(tmp_path / "B", ("se_a", "se_e")) == (False, False, [])
    assert once(tmp_path / "B", ("se_e", "se_a")) == (False, True, [])


def test_a_column_the_rows_do_not_carry_stops_the_reader(run, tmp_path):
    """MultiQC tables copied with one column dropped: HISAT2's unpaired_aligned_one, flagstat's
    secondary_passed, featureCounts' Unassigned_NoFeatures. Each reader stops, naming it.

    Mutations caught: any of the three guards replaced by a zero-filled column."""
    cases = [(run.main, "multiqc_hisat2.txt", "unpaired_aligned_one",
              lambda mqc: align.aligner_reads(mqc, run.alignment)),
             (run.main, "multiqc_samtools_flagstat.txt", "secondary_passed",
              lambda mqc: count.primary_problems(mqc, align.aligner_reads(run.main, run.alignment),
                                                 run.alignment)),
             (run.plain_runs["s0"], "multiqc_featurecounts.txt", "Unassigned_NoFeatures",
              lambda mqc: count.assigned({"s0": mqc}))]  # fmt: skip
    for n, (mqc, name, column, read) in enumerate(cases):
        copied = tmp_path / str(n)
        shutil.copytree(mqc.data, copied)
        table = pl.read_csv(copied / name, separator="\t", infer_schema=False)
        table.drop(column).write_csv(copied / name, separator="\t")
        with pytest.raises(ProcessError, match=column):
            read(multiqc_io.MultiqcRun(copied, mqc.version))


def test_i2_finds_a_category_multiqc_left_out_that_is_not_zero(tmp_path):
    """Hand-written ReadsPerGene and MultiQC plots: Multimapping, zero in every column of ours,
    is absent from all three plots (MultiQC's own rule) and passes; Ambiguous Features, 388 in
    ours, is absent from the Reverse_Stranded plot and is one I2 problem.

    Mutations caught: an absent category assumed equal to ours (no problem); an absent
    category read as missing, never as zero (a problem per plot)."""
    (tmp_path / "x.ReadsPerGene.out.tab").write_text(
        "N_unmapped\t3\t3\t3\nN_multimapping\t0\t0\t0\nN_noFeature\t4\t950\t2\n"
        "N_ambiguous\t437\t0\t388\nG1\t500\t40\t600\nG2\t56\t7\t7\n"
    )  # fmt: skip
    data = tmp_path / "multiqc_data"
    data.mkdir()
    plots = {Strandedness.UNSTRANDED: (556, 4, 437, 3), Strandedness.FORWARD:
             (47, 950, 0, 3), Strandedness.REVERSE: (607, 2, None, 3)}  # fmt: skip
    for strand, values in plots.items():
        header = ["Sample", "Overlapping Genes", "No Feature", "Ambiguous Features", "Unmapped"]
        kept = [(h, v) for h, v in zip(header[1:], values, strict=True) if v is not None]
        (data / f"{count.STAR_PLOTS[strand]}.txt").write_text(
            "\t".join(["Sample", *(h for h, _ in kept)]) + "\n"
            + "\t".join(["x", *(f"{v}.0" for _, v in kept)]) + "\n"
        )  # fmt: skip
    alignment = align.Alignment(tmp_path, True, Aligner.STAR, ("x",), False)
    aligned = pl.DataFrame({"sample": ["x"], "input_reads": [1000]})
    problems = count.star_parser_problems(alignment, multiqc_io.MultiqcRun(data, "1.35"), aligned)
    assert len(problems) == 1 and "x reverse" in problems[0] and "'N_ambiguous': 388" in problems[0]


@pytest.fixture(scope="module")
def star(tmp_path_factory):
    """STAR's index of the slice, and a folder for its alignments."""
    parent = tmp_path_factory.mktemp("star")
    ref = reference.prepare(SLICE, GTF, parent)
    return ref, reference.star_index(ref, parent, 2)[0], parent


def by_gene(frame: pl.DataFrame, column: str) -> dict[str, float]:
    """A count column by gene name."""
    ids = {v: k for k, v in NAMES.items()}
    return {ids[g]: n for g, n in frame.select("gene_id", column).rows()}


@pytest.mark.needs_star
def test_star_and_featurecounts_agree_on_single_end_reads(star):
    """se_a to se_d: STAR's column 4 equals featureCounts -s 2 on every gene (rho 1, every
    gene identical), and STAR's ReadsPerGene agrees with MultiQC's totals (I2).

    Mutations caught: STAR's column 3 read as reverse (rho falls, genes differ); MultiQC's
    Same_Stranded plot read for column 4 (I2); a category MultiQC leaves out, zero in every
    sample, read as missing instead of zero (I2); flagstat's secondary records not subtracted
    (I6)."""
    ref, index, parent = star
    samples = SAMPLES[:4]
    alignment = align.align({s: [(SE / f"{s}.fq",)] for s in samples}, parent,
                            align.AlignOptions(Aligner.STAR, index, 2))  # fmt: skip
    mqc = align.summarise(alignment, parent / "se_main", raw_config(parent / "se.json"))
    aligned = align.aligner_reads(mqc, alignment)
    assert count.star_parser_problems(alignment, mqc, aligned) == []
    assert count.primary_problems(mqc, aligned, alignment) == []
    settings = {"s2": count.FcOptions((2,) * len(samples))}
    counts = count.count(alignment, ref.gtf, parent, settings, 2)
    matrix = count.fc_genes(counts.matrix("s2"), fraction=False)
    for sample in samples:
        genes = count.star_genes(alignment.folder / f"{sample}.ReadsPerGene.out.tab")
        both = genes.join(matrix.select("gene_id", fc=align.bam_name(sample)), on="gene_id")
        agreed = concordance.agreement(both["reverse"].to_numpy(), both["fc"].to_numpy())
        assert agreed.spearman == 1.0, sample
        assert agreed.n_identical == agreed.n_genes > 0, sample


@pytest.mark.needs_star
def test_star_and_featurecounts_differ_on_pairs(star):
    """The paired fixture: STAR counts a pair only when both mates fall on one gene, while
    featureCounts --countReadPairs also takes a pair with one ambiguous mate (measured).

    | gene | STAR column 4 | featureCounts -s 2 |
    | MTFP1 | 70 | 264 |
    | SEC14L2 | 199 | 262 |
    | ENSG00000249590 | 0 | 3 |

    CCDC157's 506 pairs count 506 with --countReadPairs.

    Mutation caught: --countReadPairs dropped from -p (CCDC157 counts 1,012 reads)."""
    ref, index, parent = star
    alignment = align.align({"pe": [(PAIRED / "reads_R1.fq", PAIRED / "reads_R2.fq")]}, parent,
                            align.AlignOptions(Aligner.STAR, index, 2))  # fmt: skip
    assert alignment.paired
    counts = count.count(alignment, ref.gtf, parent, {"s2": count.FcOptions((2,), paired=True)}, 2)
    fc = by_gene(count.fc_genes(counts.matrix("s2"), fraction=False), align.bam_name("pe"))
    genes = count.star_genes(alignment.folder / "pe.ReadsPerGene.out.tab")
    ours = by_gene(genes, "reverse")
    shown = ("MTFP1", "SEC14L2", "ENSG00000249590")
    assert [ours[g] for g in shown] == [70, 199, 0]
    assert [fc[g] for g in shown] == [264, 262, 3]
    assert fc["CCDC157"] == ours["CCDC157"] == 506
    agreed = concordance.agreement(list(ours.values()), [fc[g] for g in ours])
    assert agreed.spearman is not None and agreed.spearman < 1


@pytest.mark.needs_star
def test_star_pairs_with_one_mate_mapped_pass_i6(star, tmp_path):
    """The paired fixture with 300 mates 2 made random 20-mers: STAR maps 20 pairs by mate 1
    alone and counts them mapped, so flagstat's mapped records fall short of twice STAR's
    mapped pairs by its 20 singletons, which I6 adds back; records stay twice the pairs.

    Mutation caught: I6 for STAR pairs without the singletons (4,974 against 4,994)."""
    _, index, _ = star
    lines = (PAIRED / "reads_R2.fq").read_text().splitlines()
    for n in range(300):
        lines[4 * n + 1] = "".join(random.Random(n).choices("ACGT", k=20))
        lines[4 * n + 3] = "I" * 20
    mate2 = tmp_path / "reads_R2.fq"
    mate2.write_text("\n".join(lines) + "\n")
    alignment = align.align({"pe": [(PAIRED / "reads_R1.fq", mate2)]}, tmp_path,
                            align.AlignOptions(Aligner.STAR, index, 2))  # fmt: skip
    mqc = align.summarise(alignment, tmp_path / "mqc", raw_config(tmp_path / "c.json"))
    flags = multiqc_io.table(mqc, "samtools_flagstat").row(0, named=True)
    assert int(flags["singletons_passed"]) == 20
    aligned = align.aligner_reads(mqc, alignment)
    assert aligned.rows() == [("pe", 2500, 2497, 4994)]
    assert count.primary_problems(mqc, aligned, alignment) == []
