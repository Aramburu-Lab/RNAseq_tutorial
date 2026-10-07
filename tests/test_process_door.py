"""`ctx process` at the door: every input problem listed together, and nothing aligned.

Each refusal here stops before an alignment, so these run in seconds; the end-to-end run is
``conftest.process_root`` (tests/test_process_run.py). The reference under ``tools/index/``
is the one thing a refusal may leave, since its contigs and start codons are among the checks.
"""

import gzip
import shutil
from pathlib import Path

import polars as pl
import pytest
from conftest import ctx_run
from process_root import (
    COMPARE,
    DATASET,
    ROOT,
    SLICE,
    annotation,
    dataset_file,
    plant,
    plant_compare,
)
from test_process_align import FAKES, fake_star

from ctx import runtime, tools
from ctx.dataset import load
from ctx.process import Aligner, ProcessError, align, run

PAIRED = ROOT / "tests" / "fixtures" / "chr22"
REGISTRY = ROOT / "conf" / "qc_metrics.toml"


def out(root: Path) -> Path:
    """The module's output folder."""
    return root / "outputs" / "03_process"


def fresh(tmp_path: Path) -> tuple[Path, Path]:
    """A planted root and its dataset file, for a command that must stop at the door."""
    root = tmp_path / "root"
    root.mkdir()
    plant(root)
    return root, dataset_file(root, tmp_path)


def nothing_aligned(root: Path) -> bool:
    """No page, no file beside it, and no tool folder of the dataset: only the reference."""
    folder = out(root)
    written = [p for p in folder.iterdir() if p.is_file()] if folder.is_dir() else []
    return not written and not (folder / "tools" / DATASET).exists()


def options(tmp_path: Path, **given: object) -> run.Options:
    """Options for ``run.resolve`` called in this process, HISAT2 unless told otherwise."""
    spec = given.pop("dataset_file", tmp_path / f"{DATASET}.toml")
    return run.Options(spec, given.pop("genome", SLICE), given.pop("aligner", Aligner.HISAT2),
                       **given)  # fmt: skip


def problems(root: Path, opts: run.Options) -> str:
    """The refusal ``run.resolve`` raises, as one text."""
    settings = runtime.Settings("laptop", root, 1, None, 0.0)
    with pytest.raises(ProcessError) as caught:
        run.resolve(settings, DATASET, opts)
    return str(caught.value)


def test_a_star_that_cannot_map_stops_the_command(tmp_path):
    """A fake STAR 2.7.11b on CTX_TOOLS that exits 0 and reads nothing: exit 1, stderr names
    D-036 and --aligner hisat2, and nothing is aligned (M3.4's deferred CLI test).

    Mutation caught: the probe dropped from ``run.resolve`` (STAR exits 0 without its files,
    and the refusal names neither D-036 nor HISAT2)."""
    root, spec = fresh(tmp_path)
    tools_dir = tmp_path / "tools"
    tools_dir.mkdir()
    fake_star(tools_dir, "2.7.11b", FAKES["zero counts"])
    done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                   str(spec), "--genome", str(SLICE), cwd=root,
                   env={"CTX_TOOLS": str(tools_dir), "TMPDIR": str(tmp_path)})  # fmt: skip
    assert done.returncode == 1
    assert "D-036" in done.stderr and "--aligner hisat2" in done.stderr
    assert nothing_aligned(root)


def test_a_failing_probe_a_foreign_genome_and_an_absent_comparison_are_listed_together(tmp_path):
    """A fake STAR that reads nothing, a genome whose one contig is chr2, and a
    --compare-dataset with no tables: exit 1 naming all three, nothing aligned.

    Mutations caught: a failing probe raising at once; a --compare-dataset problem raising at
    once; the reference checked only when every other input is sound."""
    root, spec = fresh(tmp_path)
    tools_dir = tmp_path / "tools"
    tools_dir.mkdir()
    fake_star(tools_dir, "2.7.11b", FAKES["zero counts"])
    genome = tmp_path / "chr2.fa"
    genome.write_text(">chr2\n" + "".join(SLICE.read_text().splitlines()[1:]) + "\n")
    done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                   str(spec), "--genome", str(genome), "--compare-dataset", "absent", cwd=root,
                   env={"CTX_TOOLS": str(tools_dir), "TMPDIR": str(tmp_path)})  # fmt: skip
    assert done.returncode == 1
    for said in ("D-036", "share no contig", "--compare-dataset absent"):
        assert said in done.stderr, said
    assert nothing_aligned(root)


def test_every_input_problem_is_listed_together(tmp_path):
    """Five problems in one command: another slug, an unknown --settings token, --multimap-max
    under HISAT2, an unknown --featured-sample, a --compare-dataset with no tables. Exit 1
    naming all five, and nothing aligned.

    Mutation caught: ``resolve`` stopping at the first problem it finds."""
    root, _ = fresh(tmp_path)
    other = dataset_file(root, tmp_path, slug="another")
    done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                   str(other), "--genome", str(SLICE), "--aligner", "hisat2", "--settings",
                   "cut,shuffled", "--multimap-max", "5", "--featured-sample", "nobody",
                   "--compare-dataset", "absent", cwd=root)  # fmt: skip
    assert done.returncode == 1
    for said in ("slug 'another'", "unknown ['shuffled']", "--multimap-max",
                 "--featured-sample nobody", "--compare-dataset absent"):  # fmt: skip
        assert said in done.stderr, said
    assert nothing_aligned(root)


def test_a_gene_the_annotation_lacks_is_refused_before_anything_is_aligned(tmp_path):
    """The Salmon twin's counts name a gene the annotation does not hold: exit 1 naming it,
    no tool folder of the dataset, no file beside the page (the issue: found only after every
    alignment, the MultiQC config already written).

    Mutation caught: ``context.unknown_genes`` left out of ``run._compare`` (the run trims,
    aligns and counts before ``context.expected`` refuses the gene)."""
    root, spec = fresh(tmp_path)
    plant_compare(root)
    path = next((root / "outputs" / "01_gather").glob(f"*_{COMPARE}_counts_v1.parquet"))
    counts = pl.read_parquet(path)
    extra = pl.DataFrame({"gene_id": ["ENSG99999999999.1"],
                          **{s: [5] for s in counts.columns[1:]}})  # fmt: skip
    pl.concat([counts, extra.cast(counts.schema)]).write_parquet(path)
    done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                   str(spec), "--genome", str(SLICE), "--aligner", "hisat2", "--settings", "",
                   "--compare-dataset", COMPARE, cwd=root)  # fmt: skip
    assert done.returncode == 1
    assert "counted genes are not in the annotation" in done.stderr
    assert "ENSG99999999999.1" in done.stderr
    assert nothing_aligned(root)


def test_a_mixed_sheet_lists_the_problems_that_read_its_samples(tmp_path):
    """se_b given a mate 2, so one run is paired; the annotation edited after the dataset file
    recorded its SHA-256; an unknown --featured-sample; an absent --compare-dataset. Exit 1
    naming all four, nothing aligned.

    Mutations caught: the layout check left to the aligner (one problem listed); the
    annotation's SHA-256 unchecked; a mixed sheet returned as no sheet, so the featured sample
    and the comparison go unchecked."""
    root, spec = fresh(tmp_path)
    sheet_path = next((root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet"))
    mate = root / "data" / "raw" / DATASET / "fastq" / "se_b_2.fastq.gz"
    with (PAIRED / "reads_R2.fq").open("rb") as src, gzip.open(mate, "wb") as dst:
        shutil.copyfileobj(src, dst)
    sheet = pl.read_parquet(sheet_path).with_columns(
        fastq_2=pl.when(pl.col("sample") == "se_b").then(pl.lit("fastq/se_b_2.fastq.gz"))
        .otherwise(pl.lit("")),
        status_2=pl.when(pl.col("sample") == "se_b").then(pl.lit("verified")).otherwise(pl.lit("")),
    )  # fmt: skip
    sheet.write_parquet(sheet_path)
    gtf = annotation(root)
    gtf.write_text(gtf.read_text() + "# edited\n")
    done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                   str(spec), "--genome", str(SLICE), "--aligner", "hisat2",
                   "--featured-sample", "nobody", "--compare-dataset", "absent",
                   cwd=root)  # fmt: skip
    assert done.returncode == 1
    for said in ("mixes single-end and paired-end", "SHA-256", "--featured-sample nobody",
                 "--compare-dataset absent"):  # fmt: skip
        assert said in done.stderr, said
    assert nothing_aligned(root)


def test_a_foreign_genome_is_refused(tmp_path):
    """A genome whose one contig is chr2: exit 1, both sides' contig names listed. The slice
    shifted by one base: exit 1, its start codons below the registry's usable band.

    Mutation caught: the start-codon band unchecked (the shifted genome runs)."""
    root, spec = fresh(tmp_path)
    sequence = "".join(SLICE.read_text().splitlines()[1:])
    for name, text in (("chr2.fa", f">chr2\n{sequence}\n"),
                       ("shifted.fa", f">chr22_slice\nA{sequence[:-1]}\n")):  # fmt: skip
        genome = tmp_path / name
        genome.write_text(text)
        done = ctx_run("--root", str(root), "process", "--dataset", DATASET, "--dataset-file",
                       str(spec), "--genome", str(genome), "--aligner", "hisat2",
                       cwd=root)  # fmt: skip
        assert done.returncode == 1, name
        assert nothing_aligned(root)
        if name == "chr2.fa":
            assert "share no contig" in done.stderr and "chr22_slice" in done.stderr
        else:
            assert "start codons read ATG" in done.stderr and "not one assembly" in done.stderr


@pytest.mark.parametrize("given", [0, 51])
def test_a_multimapper_limit_outside_star_s_range_is_a_problem(given):
    """--multimap-max 0 or 51 under STAR: outside 1 to 50, named with the range.

    Mutation caught: the range check dropped (STAR refuses 51 only after the index is built,
    and 0 discards every read)."""
    found = run.option_problems(run.Options(Path("x.toml"), SLICE, Aligner.STAR, multimap=given))
    assert found == [f"--multimap-max {given} is outside 1 to 50"]


def test_a_genome_that_is_not_a_file_is_a_problem(tmp_path):
    """--genome naming no file: a problem naming it, beside the others.

    Mutation caught: the genome check dropped (the reference step raises alone, later)."""
    missing = tmp_path / "nothere.fa"
    found = run.option_problems(run.Options(Path("x.toml"), missing, Aligner.HISAT2,
                                            settings=("cut", "bogus")))  # fmt: skip
    assert f"--genome {missing} is not a file" in found and len(found) == 2


def test_a_missing_registry_and_dataset_file_join_the_other_problems(tmp_path):
    """No registry at --registry, no --dataset-file, no genome, an unknown token: one refusal
    naming all four, where a missing registry or dataset file used to raise alone.

    Mutations caught: ``reg.load`` called outside ``_registry`` (RegistryError alone); the
    dataset file read before its existence is checked (FileNotFoundError alone)."""
    root, _ = fresh(tmp_path)
    said = problems(root, options(tmp_path, dataset_file=tmp_path / "none.toml",
                                  genome=tmp_path / "nothere.fa", settings=("cut", "bogus"),
                                  registry=tmp_path / "noreg.toml"))  # fmt: skip
    for part in ("No QC registry", "none.toml is not a file", "nothere.fa is not a file",
                 "unknown ['bogus']"):  # fmt: skip
        assert part.lower() in said.lower(), part


def test_a_subset_whose_parent_has_no_dataset_file_is_a_problem(tmp_path):
    """A dataset file declaring subset_of a dataset whose file is not beside it: a problem
    naming the file, before anything is aligned, since the page names the parent by its label
    from that file. With the file there, its label is the parent's name.

    Mutation caught: the parent left unread at the door (the run stopped after every
    alignment, when the page named the parent)."""
    root = tmp_path / "root"
    root.mkdir()
    plant(root)
    spec = dataset_file(root, tmp_path, subset_of="whole")
    said = problems(root, options(tmp_path, dataset_file=spec))
    assert "subset_of 'whole'" in said and str(tmp_path / "whole.toml") in said
    dataset_file(root, tmp_path, "whole", label="The whole cohort")
    parent = load(spec).parent()
    assert parent is not None and parent.name == "The whole cohort"


def test_a_registry_without_a_metric_the_command_reads_is_a_problem(tmp_path):
    """The registry less its reference_start_codons block: a problem naming the metric, not a
    KeyError traceback after the reference is prepared.

    Mutation caught: the metric check dropped from ``_registry`` (``_reference`` raises
    KeyError)."""
    root, spec = fresh(tmp_path)
    text = REGISTRY.read_text()
    start = text.index("[metrics.reference_start_codons]")
    rest = text.find("\n[metrics.", start + 1)
    cut = tmp_path / "registry.toml"
    cut.write_text(text[:start] + (text[rest + 1 :] if rest != -1 else ""))
    said = problems(root, options(tmp_path, dataset_file=spec, registry=cut))
    assert "has no ['reference_start_codons']" in said


def test_star_keeps_its_own_multimapper_limit_when_none_is_given(tmp_path, monkeypatch):
    """With no --multimap-max, star_align passes no --outFilterMultimapNmax, so STAR keeps its
    own 10, as the card says; with 5 it passes 5.

    Mutation caught: a default limit written into the call (--outFilterMultimapNmax 20 when
    unset)."""
    seen: list[list[str]] = []

    def record(name: str, args: list[str], **_: object) -> None:
        seen.append(args)
        raise tools.ToolError("recorded")

    monkeypatch.setattr(tools, "run", record)
    for limit in (None, 5):
        with pytest.raises(tools.ToolError, match="recorded"):
            align.star_align("s", [(SLICE,)], tmp_path, tmp_path / "out", limit, 1)
    assert "--outFilterMultimapNmax" not in seen[0]
    assert seen[1][seen[1].index("--outFilterMultimapNmax") + 1] == "5"
