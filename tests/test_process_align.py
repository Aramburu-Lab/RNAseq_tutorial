"""ctx.process.align on the chr22 fixtures: the cap, fastp, HISAT2, samtools, the STAR probe.

Tests marked needs_star skip only when the probe fails on darwin-arm64 with CI unset, and fail
everywhere else (R6.2, ``conftest.star_outcome``).
"""

import gzip
import os
import shutil
import stat
import tempfile
from pathlib import Path

import polars as pl
import pytest
from conftest import star_outcome

from ctx import tools
from ctx.process import Aligner, ProcessError, align, count, reference
from ctx.qc import Fastp, multiqc_io
from ctx.qc.reads import raw_config

ROOT = Path(__file__).resolve().parents[1]
CONFTEST = Path(__file__).resolve().with_name("conftest.py")
SLICE = ROOT / "tests" / "fixtures" / "chr22" / "slice.fa"
GTF = ROOT / "tests" / "fixtures" / "chr22" / "slice.gtf"
PAIRED = ROOT / "tests" / "fixtures" / "chr22"
SE = ROOT / "tests" / "fixtures" / "chr22_se"
SAMPLES = ("se_a", "se_b", "se_c", "se_d", "se_e")


def frame(rows: list[tuple[str, str, list[Path]]]) -> pl.DataFrame:
    """A ``reads.fastq_files`` table: (run, sample, files read 1 first) per run."""
    return pl.DataFrame(
        [(run, sample, mate, str(p)) for run, sample, paths in rows
         for mate, p in enumerate(paths, 1)],
        schema={"run_accession": pl.String, "sample": pl.String, "mate": pl.Int32,
                "path": pl.String}, orient="row",
    )  # fmt: skip


def names(path: Path) -> list[str]:
    """The read names of a FASTQ, plain or gzipped, in file order."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [line.split()[0] for n, line in enumerate(handle) if n % 4 == 0]


def report(alignment: align.Alignment, dest: Path, trimmed: align.Trimmed | None = None):
    """The alignment's MultiQC run with raw counts."""
    dest.mkdir(parents=True, exist_ok=True)
    return align.summarise(alignment, dest / "multiqc", raw_config(dest / "config.json"), trimmed)


@pytest.fixture(scope="module")
def hisat2(tmp_path_factory):
    """The slice's reference and HISAT2 index."""
    parent = tmp_path_factory.mktemp("index")
    ref = reference.prepare(SLICE, GTF, parent)
    return ref, reference.hisat2_index(ref, parent, 2)[0]


@pytest.fixture(scope="module")
def five(tmp_path_factory, hisat2):
    """The five single-end samples trimmed and aligned with HISAT2, and their MultiQC run."""
    parent = tmp_path_factory.mktemp("five")
    files = frame([(s, s, [SE / f"{s}.fq"]) for s in SAMPLES])
    trimmed = align.trim(files, parent, fastp=Fastp(), cap=None, threads=2)
    reads = align.by_sample(files, {r: trimmed.reads(r) for r in align.runs(files)})
    alignment = align.align(reads, parent, align.AlignOptions(Aligner.HISAT2, hisat2[1], 2))
    return files, trimmed, reads, alignment, report(alignment, parent, trimmed)


def fake_rows(*rows: tuple[str, int]) -> str:
    """A fake STAR's shell line writing its three files, ReadsPerGene from (row, count) pairs,
    the same count in each of the three columns."""
    lines = "".join(f"{name}\\t{n}\\t{n}\\t{n}\\n" for name, n in rows)
    return (
        f'printf "{lines}" > "${{prefix}}ReadsPerGene.out.tab"; : > "${{prefix}}Log.final.out"; '
        ': > "${prefix}Aligned.sortedByCoord.out.bam"'
    )


N_ROWS = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
FAKES = {
    # Exits 0 and writes STAR's files with every count 0: what the osx-arm64 build does
    # without GeneCounts (D-036).
    "zero counts": fake_rows(*((n, 0) for n in N_ROWS), ("probe_plus", 0), ("probe_minus", 0)),
    "nothing": ":",
    # Every column sums to 10, but one read is unmapped or none is on a gene.
    "nine on genes": fake_rows(("N_unmapped", 1), *((n, 0) for n in N_ROWS[1:]),
                               ("probe_plus", 5), ("probe_minus", 4)),
    "none on genes": fake_rows(*((n, 10 * (n == "N_noFeature")) for n in N_ROWS),
                               ("probe_plus", 0), ("probe_minus", 0)),
}  # fmt: skip


def fake_star(folder: Path, version: str, writes: str) -> Path:
    """A STAR in ``folder`` that reports ``version``, exits 0, and runs ``writes`` when given
    reads."""
    fake = folder / "STAR"
    fake.write_text(
        "#!/bin/sh\n"
        f'[ "$1" = "--version" ] && {{ echo {version}; exit 0; }}\n'
        'prefix=""; reads=""\n'
        'while [ $# -gt 0 ]; do\n'
        '  [ "$1" = "--outFileNamePrefix" ] && prefix="$2"\n'
        '  [ "$1" = "--readFilesIn" ] && reads="$2"\n'
        "  shift\n"
        "done\n"
        f'[ -n "$reads" ] && {{ {writes}; }}\n'
        "exit 0\n"
    )  # fmt: skip
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    return fake


def refused(monkeypatch, folder: Path) -> str:
    """The probe's refusal with ``folder`` as CTX_TOOLS and the temporary folder in it."""
    monkeypatch.setenv(tools.TOOLS_ENV, str(folder))
    monkeypatch.setattr(tempfile, "tempdir", str(folder))
    tools.require.cache_clear()
    align._probe.cache_clear()
    try:
        with pytest.raises(tools.ToolError) as error:
            align.star_probe()
    finally:
        tools.require.cache_clear()
        align._probe.cache_clear()
    message = str(error.value)
    assert "D-036" in message and "--aligner hisat2" in message
    assert "set CTX_TOOLS to a folder holding a STAR 2.7.11b that maps" in message
    return message


@pytest.mark.parametrize("writes", list(FAKES))
def test_a_star_that_maps_nothing_is_refused(tmp_path, monkeypatch, writes):
    """A STAR on CTX_TOOLS that reports 2.7.11b and exits 0 is refused unless all 10 reads are
    on genes, and the message names D-036, both routes and the probe's kept files.

    Mutations caught: a probe that trusts the exit code (the zero-count STAR passes); 9 of 10
    reads accepted (nine on genes passes); the on-genes rule dropped (none on genes passes);
    the probe's folder deleted while its log is named (the named folder is gone)."""
    fake = fake_star(tmp_path, "2.7.11b", FAKES[writes])
    message = refused(monkeypatch, tmp_path)
    assert str(fake) in message
    kept = Path(message.split("the probe's files are kept in ")[1].split(". D-036")[0])
    assert kept.parent == tmp_path and (kept / "probe.fastq.gz").is_file()


def test_a_star_of_another_version_or_none_is_refused_with_the_routes(tmp_path, monkeypatch):
    """STAR 2.7.10a on CTX_TOOLS, then no STAR anywhere: both refusals name D-036 and both routes.

    Mutation caught: ``tools.require`` called bare in ``star_probe`` (its message names neither)."""
    fake_star(tmp_path, "2.7.10a", ":")
    assert "reports 2.7.10a" in refused(monkeypatch, tmp_path)
    fake = tmp_path / "STAR"
    fake.rename(tmp_path / "STAR.away")
    monkeypatch.setattr(tools, "SIBLING", tmp_path / "no_reads_env")
    monkeypatch.setenv("PATH", os.pathsep.join(("/usr/bin", "/bin")))
    assert "STAR not found" in refused(monkeypatch, tmp_path)


@pytest.mark.parametrize(("system", "machine", "named"),
                         [("Darwin", "arm64", True), ("Linux", "x86_64", False),
                          ("Darwin", "x86_64", False)])  # fmt: skip
def test_the_osx_arm64_faults_are_named_on_darwin_arm64_only(monkeypatch, system, machine, named):
    """Mutation caught: the osx-arm64 build's faults given as the cause of a Linux failure."""
    monkeypatch.setattr("platform.system", lambda: system)
    monkeypatch.setattr("platform.machine", lambda: machine)
    assert ("osx-arm64" in align.refusal()) is named
    assert align.refusal().startswith(align.D036)


@pytest.mark.parametrize(
    ("failure", "system", "machine", "ci", "expected"),
    [
        (None, "Linux", "x86_64", "true", None),
        ("no reads", "Darwin", "arm64", None, "skip"),
        ("no reads", "Darwin", "arm64", "true", "fail"),
        ("no reads", "Darwin", "arm64", "", "fail"),
        ("no reads", "Darwin", "x86_64", None, "fail"),
        ("no reads", "Linux", "arm64", None, "fail"),
        ("no reads", "Linux", "x86_64", "true", "fail"),
    ],
)  # fmt: skip
def test_the_needs_star_rule(monkeypatch, failure, system, machine, ci, expected):
    """R6.2: a failing probe skips only on darwin-arm64 with CI unset; a passing one runs.

    Mutations caught: `or` for `and` (Linux skips); the CI check dropped (CI on a Mac skips);
    the machine check dropped (Rosetta skips); CI judged by its value (an empty CI skips)."""
    monkeypatch.setattr("platform.system", lambda: system)
    monkeypatch.setattr("platform.machine", lambda: machine)
    if ci is None:
        monkeypatch.delenv("CI", raising=False)
    else:
        monkeypatch.setenv("CI", ci)
    outcome = star_outcome(failure)
    assert (outcome[0] if outcome else None) == expected
    if outcome:
        assert "R6.2" in outcome[1] and failure in outcome[1]


@pytest.mark.parametrize(
    ("system", "machine", "ci", "outcome"),
    [("Darwin", "arm64", None, "skipped"), ("Darwin", "arm64", "true", "errors"),
     ("Darwin", "arm64", "", "errors"), ("Linux", "x86_64", None, "errors"),
     ("Darwin", "x86_64", None, "errors"), ("Linux", "arm64", None, "errors")],
)  # fmt: skip
def test_the_hook_skips_a_star_test_only_where_r6_2_allows(
    pytester, monkeypatch, system, machine, ci, outcome
):
    """conftest's pytest_runtest_setup on a needs_star test with STAR unreachable: no CTX_TOOLS,
    no reads environment, and a PATH without Homebrew's STAR.

    Mutation caught: the hook skipping on "fail" too (CI with a broken STAR skips)."""
    monkeypatch.delenv(tools.TOOLS_ENV, raising=False)
    monkeypatch.setattr(tools, "SIBLING", pytester.path / "no_reads_env")
    monkeypatch.setenv("PATH", os.pathsep.join(("/usr/bin", "/bin")))
    monkeypatch.setattr("platform.system", lambda: system)
    monkeypatch.setattr("platform.machine", lambda: machine)
    if ci is None:
        monkeypatch.delenv("CI", raising=False)
    else:
        monkeypatch.setenv("CI", ci)
    tools.require.cache_clear()
    align._probe.cache_clear()
    pytester.makeconftest(CONFTEST.read_text())
    pytester.makepyfile("import pytest\n\n\n@pytest.mark.needs_star\ndef test_star():\n    pass\n")
    result = pytester.runpytest(
        "-p", "no:cacheprovider", "-rsE", "-W", "ignore::pytest.PytestUnknownMarkWarning"
    )
    result.assert_outcomes(**{outcome: 1})
    result.stdout.fnmatch_lines(["*R6.2*STAR not found*D-036*"])


def test_the_cap_draws_a_seeded_random_subset(tmp_path):
    """500 of se_a's 1,000 reads: the same 500 every time, not the first 500, and fastp's
    report (through MultiQC) gives the realised count.

    Mutation caught: `seqkit head` without `seqkit sample` (the first 500)."""
    files = frame([("se_a", "se_a", [SE / "se_a.fq"])])
    first = align.trim(files, tmp_path / "one", fastp=Fastp(), cap=500, threads=1)
    again = align.trim(files, tmp_path / "two", fastp=Fastp(), cap=500, threads=1)
    kept, source = names(first.raw("se_a")[0]), names(SE / "se_a.fq")
    assert first.raw("se_a")[0].parent == first.folder / align.CAPPED
    assert len(kept) == 500 and len(set(kept)) == 500
    assert kept == names(again.raw("se_a")[0])
    assert set(kept) <= set(source) and kept != source[:500]
    assert [n for n in source if n in set(kept)] == kept
    mqc = multiqc_io.run([first.reports], tmp_path / "mqc", modules=("fastp",),
                         config=raw_config(tmp_path / "config.json"))  # fmt: skip
    assert align.fastp_reads(mqc, ["se_a"]).rows() == [("se_a", 500, 500)]


def test_a_run_within_the_cap_is_read_whole(tmp_path):
    """Mutation caught: a cap at or above the run's size still sampled (the run is copied into
    capped/ and fastp reads the copy, not the run's own file)."""
    files = frame([("se_b", "se_b", [SE / "se_b.fq"])])
    trimmed = align.trim(files, tmp_path, fastp=Fastp(), cap=1000, threads=1)
    assert trimmed.raw("se_b") == (SE / "se_b.fq",)
    assert len(names(trimmed.reads("se_b")[0])) == 1000


def test_the_margin_reaches_the_cap_where_cap_over_reads_falls_short(tmp_path):
    """Seed 42 at p = 100 / 1,000 keeps 90 of se_a's reads (measured with seqkit 2.13.0); the
    5 SD margin samples at 0.15, so the cap keeps exactly 100.

    Mutation caught: p = cap / reads, without the margin (90 reads)."""
    assert align.proportion(100, 1000) == pytest.approx(0.15)
    files = frame([("se_a", "se_a", [SE / "se_a.fq"])])
    trimmed = align.trim(files, tmp_path, fastp=Fastp(), cap=100, threads=1)
    assert len(names(trimmed.raw("se_a")[0])) == 100


def test_a_cap_just_below_the_run_is_its_first_reads(tmp_path):
    """Cap 990 of 1,000: p reaches 1, every read is sampled and `seqkit head` keeps the first
    990, as ``cap_run`` says; no random subset.

    Mutation caught: p not capped at 1 (seqkit refuses a proportion above 1)."""
    assert align.proportion(990, 1000) == 1.0
    files = frame([("se_a", "se_a", [SE / "se_a.fq"])])
    trimmed = align.trim(files, tmp_path, fastp=Fastp(), cap=990, threads=1)
    assert names(trimmed.raw("se_a")[0]) == names(SE / "se_a.fq")[:990]


def test_relative_paths_reach_every_tool(tmp_path, monkeypatch, hisat2):
    """From a working folder: cap_run, align, bam_stats and count given relative paths, each
    tool running in its own folder, write where they say.

    Mutations caught: dest not resolved in cap_run (seqkit writes under logs/, exit 0); the
    reads not resolved in _joined (hisat2 exits 1); the BAM not resolved in bam_stats
    (samtools exits 1); the annotation not resolved for featureCounts, which runs in the
    alignment folder (exit 1)."""
    monkeypatch.chdir(tmp_path)
    shutil.copy(SE / "se_a.fq", tmp_path / "se_a.fq")
    shutil.copy(GTF, tmp_path / "slice.gtf")
    capped = align.cap_run([Path("se_a.fq")], 500, Path("capdest"), "x")
    assert capped == (tmp_path / "capdest" / "x.fastq.gz",) and len(names(capped[0])) == 500
    assert align.count_reads(Path("se_a.fq")) == 1000
    alignment = align.align({"se_a": [(Path("se_a.fq"),)]}, Path("parent"),
                            align.AlignOptions(Aligner.HISAT2, hisat2[1], 1))  # fmt: skip
    (tmp_path / "bams").mkdir()
    shutil.copy(alignment.bam("se_a"), tmp_path / "bams" / "t.bam")
    align.bam_stats(Path("bams/t.bam"))
    assert "\n1000 + 0 primary\n" in (tmp_path / "bams" / "t.flagstat").read_text()
    counts = count.count(alignment, Path("slice.gtf"), Path("parent"),
                         {"s0": count.FcOptions((0,))}, 1)  # fmt: skip
    # Untrimmed se_a: 562 assigned at -s 0, against 563 after fastp.
    assert count.fc_genes(counts.matrix("s0"), fraction=False)["se_a.bam"].sum() == 562


def test_a_fastq_of_partial_records_is_refused(tmp_path):
    """Mutations caught: a line count not divisible by 4 floored to reads; a last line without
    its newline left uncounted (the gzipped two reads read as 7 lines)."""
    broken = tmp_path / "broken.fq"
    broken.write_text("".join((SE / "se_a.fq").read_text().splitlines(keepends=True)[:6]))
    with pytest.raises(ProcessError, match="holds 6 lines, not 4 per read"):
        align.count_reads(broken)
    packed = tmp_path / "two.fq.gz"
    with gzip.open(packed, "wt") as handle:
        handle.write("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII")
    assert align.count_reads(packed) == 2


def test_both_mates_keep_the_same_reads(tmp_path):
    """The cap draws one set of pairs: read 1 and read 2 keep the same 1,000 of 2,500.

    Mutation caught: a seed per mate (the mates fall out of step)."""
    files = frame([("pe", "pe", [PAIRED / "reads_R1.fq", PAIRED / "reads_R2.fq"])])
    trimmed = align.trim(files, tmp_path, fastp=Fastp(), cap=1000, threads=1)
    one, two = (names(p) for p in trimmed.raw("pe"))
    assert len(one) == len(two) == 1000
    assert [n.removesuffix("/1") for n in one] == [n.removesuffix("/2") for n in two]
    assert [p.name for p in trimmed.reads("pe")] == ["pe_1.fastq.gz", "pe_2.fastq.gz"]


def test_the_aligner_reads_what_fastp_passed(tmp_path, hisat2):
    """A run of se_d's 101 nt reads and se_a's 76 nt ones: --cut_tail trims se_d's to 81 nt,
    so -l 80 passes those 1,000 of 2,000, and HISAT2 reads exactly them (I1); the untrimmed
    setting reads all 2,000, fastp's reads before filtering.

    Mutation caught: aligning the raw reads (`Trimmed.reads` returning the gathered files)."""
    mixed = tmp_path / "mixed.fastq"
    mixed.write_text((SE / "se_d.fq").read_text() + (SE / "se_a.fq").read_text())
    files = frame([("mix", "mix", [mixed])])
    trimmed = align.trim(files, tmp_path, fastp=Fastp(l=80, cut_tail=True), cap=None, threads=1)
    reads = align.by_sample(files, {"mix": trimmed.reads("mix")})
    alignment = align.align(reads, tmp_path, align.AlignOptions(Aligner.HISAT2, hisat2[1], 1))
    mqc = report(alignment, tmp_path, trimmed)
    fastp = align.fastp_reads(mqc, ["mix"])
    assert fastp.rows() == [("mix", 2000, 1000)]
    aligned = align.aligner_reads(mqc, alignment)
    assert aligned["input_reads"].to_list() == [1000]
    passed = count.expected_input(fastp, files, "passed_reads", paired=False)
    assert count.input_problems(aligned, passed) == []
    raw = count.expected_input(fastp, files, "raw_reads", paired=False)
    assert count.input_problems(aligned, raw) == [
        "I1: mix's aligner read 1000 reads of the 2000 it was given"
    ]
    untrimmed = align.align({"mix": [trimmed.raw("mix")]}, tmp_path / "untrimmed",
                            align.AlignOptions(Aligner.HISAT2, hisat2[1], 1))  # fmt: skip
    as_given = align.aligner_reads(report(untrimmed, tmp_path / "untrimmed"), untrimmed)
    assert count.input_problems(as_given, raw) == []


def test_a_sample_of_two_runs_is_one_alignment(tmp_path, hisat2):
    """Runs se_a and se_b as one sample: one HISAT2 row, reading both runs' 2,000 reads.

    Mutation caught: only a sample's first run joined (1,000 reads)."""
    files = frame([("run1", "x", [SE / "se_a.fq"]), ("run2", "x", [SE / "se_b.fq"])])
    trimmed = align.trim(files, tmp_path, fastp=Fastp(), cap=None, threads=1)
    reads = align.by_sample(files, {r: trimmed.reads(r) for r in ("run1", "run2")})
    assert list(reads) == ["x"] and len(reads["x"]) == 2
    alignment = align.align(reads, tmp_path, align.AlignOptions(Aligner.HISAT2, hisat2[1], 1))
    mqc = report(alignment, tmp_path, trimmed)
    aligned = align.aligner_reads(mqc, alignment)
    assert aligned.select("sample", "input_reads").rows() == [("x", 2000)]
    expected = count.expected_input(align.fastp_reads(mqc, ["run1", "run2"]), files,
                                    "passed_reads", paired=False)  # fmt: skip
    assert expected == {"x": 2000}
    assert count.input_problems(aligned, expected) == []


def test_unspliced_alignment_leaves_spliced_reads_unaligned(tmp_path, five, hisat2):
    """se_a with --no-spliced-alignment: 131 of 1,000 reads unaligned, against none spliced
    (measured 2026-10-01, fastp-trimmed se_a on the plain slice index).

    Mutation caught: the spliced switch ignored (0 unaligned)."""
    _, _, reads, _, main = five
    options = align.AlignOptions(Aligner.HISAT2, hisat2[1], 1, spliced=False)
    unspliced = align.align({"se_a": reads["se_a"]}, tmp_path, options)
    column = "unpaired_aligned_none"
    table = multiqc_io.numeric(multiqc_io.table(report(unspliced, tmp_path), "hisat2"), [column])
    before = multiqc_io.numeric(multiqc_io.table(main, "hisat2"), [column])
    assert table[column].to_list() == [131.0]
    assert before.filter(pl.col("Sample") == "se_a")[column].to_list() == [0.0]


def test_flagstat_primary_records_equal_hisat2_reads(five):
    """I6 under HISAT2: se_c's 864 of 1,000 aligned, and every read in the BAM.

    Mutation caught: --no-unal added to HISAT2 (flagstat's total falls to the aligned reads)."""
    _, _, _, alignment, mqc = five
    aligned = align.aligner_reads(mqc, alignment)
    assert aligned.filter(pl.col("sample") == "se_c").rows() == [("se_c", 1000, 864, 864)]
    assert count.primary_problems(mqc, aligned, alignment) == []


def test_runs_of_both_layouts_are_refused(tmp_path, hisat2):
    """Mutation caught: no layout check (HISAT2 would take one sample's pair as two runs)."""
    reads = {"a": [(SE / "se_a.fq",)],
             "b": [(PAIRED / "reads_R1.fq", PAIRED / "reads_R2.fq")]}  # fmt: skip
    with pytest.raises(ProcessError, match="all single-end or all paired-end"):
        align.align(reads, tmp_path, align.AlignOptions(Aligner.HISAT2, hisat2[1], 1))


@pytest.mark.parametrize(
    ("aligner", "setting", "named"),
    [(Aligner.STAR, {"spliced": False}, "spliced=False"),
     (Aligner.HISAT2, {"multimap": 1}, "multimap=1")],
)  # fmt: skip
def test_a_setting_the_aligner_cannot_apply_is_refused(tmp_path, aligner, setting, named):
    """STAR has no unspliced switch, and HISAT2 here no multimapper limit: refused before any
    tool runs, never keyed as a second folder holding the same alignment.

    Mutations caught: either refusal dropped (a second folder, labelled falsely, is built)."""
    options = align.AlignOptions(aligner, tmp_path / "no_index", 1, **setting)
    with pytest.raises(ProcessError, match=named):
        align.align({"se_a": [(SE / "se_a.fq",)]}, tmp_path, options)
    assert list(tmp_path.iterdir()) == []


@pytest.fixture(scope="module")
def paired(tmp_path_factory, hisat2):
    """The paired fixture trimmed and aligned with HISAT2, and its MultiQC run."""
    parent = tmp_path_factory.mktemp("paired")
    files = frame([("pe", "pe", [PAIRED / "reads_R1.fq", PAIRED / "reads_R2.fq"])])
    trimmed = align.trim(files, parent, fastp=Fastp(), cap=None, threads=2)
    reads = align.by_sample(files, {"pe": trimmed.reads("pe")})
    alignment = align.align(reads, parent, align.AlignOptions(Aligner.HISAT2, hisat2[1], 2))
    return files, reads, alignment, report(alignment, parent, trimmed)


def test_pairs_pass_i1_and_i6_under_hisat2(paired):
    """The paired fixture: HISAT2 reads fastp's passed mates halved (I1); flagstat's primary
    records are twice its pairs, and its mapped records twice the pairs aligned concordantly
    or discordantly plus the mates aligned alone (I6).

    Mutations caught: fastp's mates not halved for pairs (I1); the mates aligned alone left
    out of HISAT2's mapped mates (I6 finds 3 records more than the aligner's)."""
    files, _, alignment, mqc = paired
    aligned = align.aligner_reads(mqc, alignment)
    pairs, mates = aligned.select("input_reads", "mapped_mates").row(0)
    fastp = align.fastp_reads(mqc, ["pe"])
    assert fastp["passed_reads"].item() == 2 * pairs
    assert count.input_problems(aligned, count.expected_input(fastp, files, "passed_reads",
                                                              paired=True)) == []  # fmt: skip
    assert count.input_problems(aligned, count.expected_input(fastp, files, "passed_reads",
                                                              paired=False)) != []  # fmt: skip
    table = multiqc_io.numeric(multiqc_io.table(mqc, "hisat2"), align.HISAT2_ALONE)
    assert table.select(align.HISAT2_ALONE).sum_horizontal().item() > 0
    assert mates < 2 * pairs
    assert count.primary_problems(mqc, aligned, alignment) == []


def test_paired_i6_sees_a_bam_without_its_unaligned_mates(tmp_path, monkeypatch, paired, hisat2):
    """--no-unal added to HISAT2's options: a new keyed folder (the options are keyed), and I6
    finds flagstat's records short of twice the pairs.

    Mutations caught: HISAT2's options left out of the key (the old folder is reused); I6
    skipped for pairs (no problem)."""
    _, reads, before, _ = paired
    monkeypatch.setattr(align, "HISAT2_OPTIONS", (*align.HISAT2_OPTIONS, "--no-unal"))
    alignment = align.align(reads, tmp_path, align.AlignOptions(Aligner.HISAT2, hisat2[1], 2))
    assert alignment.built and alignment.folder.name != before.folder.name
    mqc = report(alignment, tmp_path)
    problems = count.primary_problems(mqc, align.aligner_reads(mqc, alignment), alignment)
    assert len(problems) == 1 and problems[0].startswith("I6: pe:")


def test_the_derived_limit_moves_multimappers_to_too_many_loci():
    """se_c's STAR numbers (49 multimapped, 921 unique) had the limit been 1; and 49 of 160
    reads, a rounding tie, printed as STAR prints it: 30.63, where 100 * 49 / 160 gives 30.62.

    Mutations caught: unique reads moved too; the earlier too-many-loci reads dropped; the
    percent computed as round(100 * v / total, 2)."""
    row = {"total_reads": 1000.0, "uniquely_mapped": 921.0, "multimapped": 49.0,
           "multimapped_toomany": 3.0, "mapped": 970.0}  # fmt: skip
    derived = align.derived_limit_one(row)
    fields = ("uniquely_mapped", "multimapped", "multimapped_toomany", "mapped")
    assert [derived[f] for f in fields] == [921.0, 0.0, 52.0, 921.0]
    assert [derived[f"{f}_percent"] for f in fields[1:]] == [0.0, 5.2, 92.1]
    tie = {**row, "total_reads": 160.0, "uniquely_mapped": 111.0, "multimapped_toomany": 0.0}
    assert align.derived_limit_one(tie)["multimapped_toomany_percent"] == 30.63


def star_limits(parent: Path, index: Path, sample: str, fastq: Path):
    """One sample through STAR at the default limit and at 1, and the MultiQC run of each."""
    out = {}
    for limit in (None, 1):
        options = align.AlignOptions(Aligner.STAR, index, 2, multimap=limit)
        alignment = align.align({sample: [(fastq,)]}, parent, options)
        out[limit] = alignment, report(alignment, parent / f"{sample}_limit_{limit}")
    return out


@pytest.fixture(scope="module")
def star_index(tmp_path_factory):
    """A folder for STAR's alignments, and STAR's index of the slice."""
    parent = tmp_path_factory.mktemp("star")
    return parent, reference.star_index(reference.prepare(SLICE, GTF, parent), parent, 2)[0]


@pytest.fixture(scope="module")
def star_se_c(star_index):
    """se_c through STAR at the default limit and at 1."""
    return star_limits(*star_index, "se_c", SE / "se_c.fq")


@pytest.fixture(scope="module")
def star_tie(star_index, star_se_c):
    """160 of se_c's reads, its 49 multimappers and its first 111 unique reads by NH tag, at
    the default limit and at 1: 49 / 160 is a rounding tie (30.625 %)."""
    parent, index = star_index
    records = parent / "se_c.records.sam"
    tools.run("samtools", ["view", "-F", "260", str(star_se_c[None][0].bam("se_c"))],
              cwd=parent, log=parent / "view.log", stdout=records)  # fmt: skip
    hits = {r.split("\t")[0]: int(r.split("\tNH:i:")[1].split("\t")[0])
            for r in records.read_text().splitlines()}  # fmt: skip
    unique = [n for n in (h.removeprefix("@") for h in names(SE / "se_c.fq")) if hits.get(n) == 1]
    keep = {n for n, h in hits.items() if h > 1} | set(unique[:111])
    lines = (SE / "se_c.fq").read_text().splitlines(keepends=True)
    tie = parent / "tie.fq"
    tie.write_text("".join("".join(lines[i : i + 4]) for i in range(0, len(lines), 4)
                           if lines[i].split()[0].removeprefix("@") in keep))  # fmt: skip
    return star_limits(parent, index, "tie", tie)


@pytest.mark.needs_star
def test_the_probe_maps_its_reads():
    """Mutation caught: a probe that refuses a STAR that maps, such as the one-exon probe
    genome that aborts linux STAR's genomeGenerate (D-036); an error with CI set, a skip on
    darwin-arm64 with CI unset, as R6.2 allows."""
    align.star_probe()


STAR_FIELDS = ("total_reads", "uniquely_mapped", "multimapped", "multimapped_toomany", "mapped",
               "multimapped_percent", "multimapped_toomany_percent", "mapped_percent")  # fmt: skip


@pytest.mark.needs_star
@pytest.mark.parametrize(("runs", "sample"), [("star_se_c", "se_c"), ("star_tie", "tie")])
def test_the_derived_limit_equals_a_real_run_at_one(request, runs, sample):
    """se_c, and its 160-read tie, at --outFilterMultimapNmax 1: the 49 multimappers become
    too many loci, unique reads stay, every field and percent equals STAR's own (the tie's
    30.63 %), and no gene count moves (GeneCounts never counts a multimapper).

    Mutations caught: --outFilterMultimapNmax not passed (the rerun keeps 49 multimapped); the
    percent as round(100 * v / total, 2) (the tie's 30.62 against STAR's 30.63)."""
    star_runs, rows = request.getfixturevalue(runs), {}
    for limit, (_, mqc) in star_runs.items():
        table = multiqc_io.numeric(multiqc_io.table(mqc, "star"), STAR_FIELDS)
        rows[limit] = table.row(0, named=True)
    derived = align.derived_limit_one(rows[None])
    assert {f: derived[f] for f in STAR_FIELDS} == {f: rows[1][f] for f in STAR_FIELDS}
    assert (
        rows[None]["multimapped"],
        rows[1]["multimapped_toomany"] - rows[None]["multimapped_toomany"],
    ) == (49.0, 49.0)
    genes = [count.star_genes(a.folder / f"{sample}.ReadsPerGene.out.tab")
             for a, _ in star_runs.values()]  # fmt: skip
    assert genes[0].equals(genes[1])


@pytest.mark.needs_star
def test_flagstat_primary_records_equal_star_reads(star_se_c):
    """I6 under STAR: se_c's primary records are its 1,000 reads, 970 mapped, 97.00 %.

    Mutation caught: --outSAMunmapped Within dropped (flagstat's total falls to 970)."""
    alignment, mqc = star_se_c[None]
    aligned = align.aligner_reads(mqc, alignment)
    assert aligned.rows() == [("se_c", 1000, 970, 970)]
    assert count.primary_problems(mqc, aligned, alignment) == []
    flags = multiqc_io.table(mqc, "samtools_flagstat").row(0, named=True)
    primary = int(flags["mapped_passed"]) - int(flags["secondary_passed"])
    assert round(100 * primary / 1000, 2) == 97.00
