"""ctx quant's laptop path: the transcriptome, Salmon's probe, M3's reads, the planted truth.

Salmon runs here on the chr22_se reads, planted as M3's trim folder, against the chr22
fixture's slice with its genome as decoys. Each test names the mutation of
``src/ctx/quant`` it was seen to fail on (2026-10-09).
"""

import gzip
import json
import os
import shutil
import stat
from dataclasses import replace
from pathlib import Path

import numpy as np
import polars as pl
import pytest
import typer
from process_root import SLICE
from quant_root import DATASET, RNASEQ, quant_root

from ctx import runtime, tools
from ctx.cli import app, parameter_cards
from ctx.process.cache import MANIFEST
from ctx.process.reference import records
from ctx.quant import CountsFrom, QuantError, SalmonOptions, outdir, salmon, tximport
from ctx.quant.laptop import counts as laptop_counts
from ctx.quant.run import Options, resolve, run, run_salmon

FIXTURE = Path(__file__).resolve().parent / "fixtures"
TRIM = "trim_000000000000"
# SEC14L2, ENSG00000249590 and MTFP1 overlap on one strand; only their sum is planted.
TRIO = ("ENSG00000100003.18", "ENSG00000249590.7", "ENSG00000242114.8")


def plant_trim(root: Path, *, paired: str = "") -> None:
    """M3's newest folders table and its trim folder: each run's reads under its accession."""
    sheet_path = next((root / "outputs" / "01_gather").glob(f"*_{DATASET}_samplesheet_v1.parquet"))
    sheet = pl.read_parquet(sheet_path)
    if paired:
        mate = pl.when(pl.col("sample") == paired).then(pl.lit("fastq/x_2.fastq.gz"))
        sheet = sheet.with_columns(fastq_2=mate.otherwise(pl.col("fastq_2")))
        sheet.write_parquet(sheet_path)
    process = root / "outputs" / "03_process"
    folder = process / "tools" / DATASET / TRIM
    (folder / "reads").mkdir(parents=True)
    files = {}
    for accession, sample in sheet.select("run_accession", "sample").rows():
        dest = folder / "reads" / f"{accession}.fastq.gz"
        shutil.copy(root / "data" / "raw" / DATASET / "fastq" / f"{sample}.fastq.gz", dest)
        files[f"reads/{dest.name}"] = dest.stat().st_size
    (folder / MANIFEST).write_text(json.dumps({"key": "x", "inputs": {}, "files": files}))
    pl.DataFrame({"step": ["trim"], "folder": [TRIM], "built": [True]}).write_parquet(
        process / f"2026-10-09_process_{DATASET}_folders_v1.parquet")  # fmt: skip


def settings(root: Path) -> runtime.Settings:
    """A laptop profile on one thread."""
    return runtime.Settings("laptop", root, 1, None, 0.0)


@pytest.fixture(scope="module")
def laptop(tmp_path_factory):
    """ctx quant's laptop path run twice on one root: (root, first files, second run's)."""
    root = tmp_path_factory.mktemp("laptop")
    spec = quant_root(root)
    plant_trim(root)
    options = Options(spec, None, genome=SLICE)
    group = typer.main.get_command(app)
    cards = parameter_cards(group, group.commands["quant"])  # type: ignore[attr-defined]
    first = run(settings(root), cards, "ctx quant", DATASET, options)
    second = resolve(settings(root), DATASET, options)
    return root, first, second, options


def newest(root: Path, what: str) -> pl.DataFrame:
    """One table ctx quant wrote."""
    return pl.read_parquet(next((root / "outputs" / "04_quant").glob(f"*_{what}_v1.parquet")))


def test_transcriptome_is_the_fixtures():
    """Spliced from the slice and its GTF, every transcript equals the fixture builder's.

    The builder spliced ``transcripts.fa`` with its own code (make_chr22_fixture.py), 11 of
    36 transcripts on the minus strand. Mutations: a minus-strand transcript not reverse-
    complemented; exons joined in GTF order, which GENCODE lists 5' to 3', so a minus-strand
    transcript's exons come reversed; exon starts read 0-based.
    """
    ours = dict(salmon.spliced(SLICE, FIXTURE / "chr22" / "slice.gtf"))
    theirs = dict(records(FIXTURE / "chr22" / "transcripts.fa"))
    assert len(ours) == 36
    assert ours == theirs


def test_laptop_page_and_reuse(laptop):
    """The laptop path writes M4's page and figures; a rerun reuses Salmon's folders.

    Mutations: the transcript TPM table not written (``transcript_tpm`` skipped): no isoform
    figure; Salmon's quantification keyed without its arguments (``quantify``'s inputs): the
    run with ``-l U`` reuses the ``-l A`` folder.
    """
    root, written, inputs, options = laptop
    names = {p.name for p in written}
    assert written[0].name.endswith(f"_quant_{DATASET}_overview_v1.html")
    for what in ("fates", "correlation", "detection", "units", "isoforms"):
        assert any(f"_{what}_v1_white.pdf" in n for n in names), what
    out = root / "outputs" / "04_quant"
    quant = sorted((out / "tools" / DATASET).glob("salmon_quant_*"))
    assert len(quant) == 1
    assert len(sorted((out / "tools" / "index").glob("salmon_index_*"))) == 1
    assert "start codons read ATG" in written[0].read_text()
    _, again = run_salmon(inputs, options, out, 1)
    assert (again.folder, again.built, again.index_built) == (quant[0], False, False)
    _, other = run_salmon(inputs, replace(options, salmon=SalmonOptions(libtype="U")), out, 1)
    assert other.folder != quant[0]
    assert other.built
    assert not other.index_built


def test_laptop_planted_truth(laptop):
    """Salmon here calls the planted library types and recovers the planted gene totals.

    SR in four samples and U in se_e. In se_a, se_b and se_d the genes outside the trio sum
    to their planted reads within 1 %, and all but one estimate theirs within 1 % or 1 read:
    Salmon at -p 1 is not reproducible run to run, and in 1 of 6 runs one gene planted at 4
    reads came out 8.663 (2026-10-09).
    Mutation: Salmon run with ``-l U`` whatever ``--libtype`` says (``quant_args``): every
    sample's library type reads U.
    """
    root, *_ = laptop
    library = newest(root, "library")
    assert dict(library.select("sample", "library_type").iter_rows()) == {
        "se_a": "SR", "se_b": "SR", "se_c": "SR", "se_d": "SR", "se_e": "U"}  # fmt: skip
    truth = (pl.read_csv(FIXTURE / "chr22_se" / "truth.tsv", separator="\t")
             .group_by("sample", "gene_id").agg(pl.col("n_reads").sum()))  # fmt: skip
    assert (newest(root, "sums")["relative_difference"].abs() < 1e-5).all()  # 0.001 rounding
    quant = next((root / "outputs" / "04_quant" / "tools" / DATASET).glob("salmon_quant_*"))
    samples = ("se_a", "se_b", "se_c", "se_d", "se_e")
    port = tximport.summarize(outdir.quants(quant, samples), outdir.tx2gene(quant))
    long = port.frame(port.counts).unpivot(index="gene_id", variable_name="sample",
                                           value_name="estimated")  # fmt: skip
    joined = long.join(truth, on=["sample", "gene_id"], how="left").fill_null(0)
    for sample in ("se_a", "se_b", "se_d"):
        mine = joined.filter(pl.col("sample") == sample)
        alone = mine.filter(~pl.col("gene_id").is_in(TRIO) & (pl.col("n_reads") > 0))
        assert alone.height > 0
        assert alone["estimated"].sum() == pytest.approx(alone["n_reads"].sum(), rel=0.01)
        off = (alone["estimated"] - alone["n_reads"]).abs() > np.maximum(1, 0.01 * alone["n_reads"])
        assert off.sum() <= 1


def test_door(tmp_path):
    """Every laptop-path problem is listed together, and nothing is written.

    Mutation: ``resolve`` checking the reads only when the genome is a file: the paired-end
    run and the manifest's wrong size go unlisted.
    """
    root = tmp_path / "root"
    spec = quant_root(root)
    plant_trim(root, paired="se_b")
    reads = next((root / "outputs" / "03_process" / "tools" / DATASET / TRIM / "reads").glob("*"))
    reads.write_bytes(b"short")
    with pytest.raises(QuantError) as caught:
        resolve(settings(root), DATASET, Options(spec, None, genome=tmp_path / "no.fa"))
    message = str(caught.value)
    assert "no.fa is not a file" in message
    assert "is paired-end" in message
    assert "is not the file M3's manifest lists" in message
    with pytest.raises(QuantError, match="one of the two"):
        resolve(settings(root), DATASET, Options(spec, RNASEQ, genome=SLICE))
    assert not (root / "outputs" / "04_quant").exists()


def test_door_without_m3(tmp_path):
    """No M3 run to take trimmed reads from is a problem, named with its fix.

    Mutation: ``trim_folder`` returning no problem when the folders table is missing: the
    run starts and fails later inside Salmon.
    """
    root = tmp_path / "root"
    spec = quant_root(root)
    with pytest.raises(QuantError, match="run ctx process first"):
        resolve(settings(root), DATASET, Options(spec, None, genome=SLICE))


def test_probe_refuses_a_salmon_that_maps_nothing(tmp_path, monkeypatch):
    """A Salmon that reports 2.7.0 and exits 0 having assigned no read is refused.

    Mutation: ``_probe`` checking only that Salmon exited 0: the fake passes.
    """
    fake = tmp_path / "bin" / "salmon"
    fake.parent.mkdir()
    fake.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "--version" ]; then echo "salmon 2.7.0"; exit 0; fi\n'
        'if [ "$1" = "quant" ]; then mkdir -p out; printf "Name\\tLength\\tEffectiveLength'
        '\\tTPM\\tNumReads\\nprobe0\\t1000\\t800\\t0\\t0\\n" > out/quant.sf; fi\n'
        "exit 0\n"
    )
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv(tools.TOOLS_ENV, str(fake.parent))
    tools.require.cache_clear()
    salmon._probe.cache_clear()
    try:
        with pytest.raises(tools.ToolError, match="assigned 0 of the probe's 10 reads"):
            salmon.salmon_probe()
    finally:
        tools.require.cache_clear()
        salmon._probe.cache_clear()
    assert os.environ[tools.TOOLS_ENV] == str(fake.parent)


def test_probe_reads_are_exact_copies(tmp_path):
    """The probe's reads are cut from its transcripts, so a working Salmon must take all 10.

    Mutation: reads drawn from a third random sequence: no read is a substring.
    """
    salmon._probe_inputs(tmp_path)
    transcripts = "".join(s for _, s in records(tmp_path / "transcripts.fa"))
    with gzip.open(tmp_path / "probe.fastq.gz", "rt") as handle:
        reads = handle.read().splitlines()[1::4]
    assert len(reads) == salmon.PROBE_READS
    assert all(r in transcripts for r in reads)


def test_counts_are_length_scaled():
    """The matrix M4's correlation figure reads is the port's length-scaled counts, rounded.

    On the planted set, whose effective lengths differ by sample, estimated and length-scaled
    counts differ. Mutation: ``laptop.counts`` built from countsFromAbundance "no".
    """
    planted = FIXTURE / "quant" / "planted"
    samples = ("s1", "s2", "s3", "s4")
    port = tximport.summarize(outdir.quants(planted, samples), outdir.tx2gene(planted))
    made = laptop_counts(port, {"G1": "one"}, {})
    expected = np.rint(port.counts_from(CountsFrom.LENGTH_SCALED_TPM))
    np.testing.assert_array_equal(made.matrix.select(samples).to_numpy(), expected)
    assert not np.array_equal(expected, np.rint(port.counts))
    assert made.genes.filter(pl.col("gene_id") == "G1")["gene_name"].item() == "one"
