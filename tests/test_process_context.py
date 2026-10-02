"""ctx.process.context: another dataset's Salmon run, context beside a partial reference."""

from pathlib import Path

import polars as pl
import pytest
from process_root import COMPARE, DAY, GTF, OFF_CONTIG, SAMPLES, annotation, plant_compare

from ctx.process import ProcessError, context
from ctx.qc import registry
from ctx.qc.checkpoint import gene_types

REGISTRY = registry.load(Path(__file__).resolve().parents[1] / "conf" / "qc_metrics.toml")
CONTIGS = ("chr22_slice",)


def planted(tmp_path: Path) -> tuple[Path, pl.DataFrame]:
    """A root holding the compare dataset's tables, and the annotation's genes with their contig."""
    plant_compare(tmp_path)
    return tmp_path, gene_types(annotation(tmp_path))


def test_the_expected_share_is_salmon_on_the_contigs_over_all(tmp_path):
    """900 of each sample's 1,000 counts lie on the slice's genes, 100 on a gene elsewhere: 0.9,
    and 0.72 of every read at Salmon's 80 % mapped; se_e alone is called unstranded.

    Mutations caught: the share taken over the contigs' counts alone (1.0); the off-contig gene
    counted as on the contigs (1.0); the read share not scaled by Salmon's mapped percent."""
    root, genes = planted(tmp_path)
    tables, problems = context.read(root, COMPARE, list(SAMPLES), REGISTRY)
    assert problems == [] and tables is not None
    table = context.expected(tables, genes, CONTIGS, REGISTRY)
    assert table["expected_share"].to_list() == pytest.approx([0.9] * len(SAMPLES))
    assert table["expected_read_share"].to_list() == pytest.approx([0.72] * len(SAMPLES))
    calls = dict(table.select("sample", "salmon_call").rows())
    assert calls == {**dict.fromkeys(SAMPLES, "reverse"), "se_e": "unstranded"}


def test_only_this_runs_samples_are_read(tmp_path):
    """The other dataset holds a sample this run lacks: its column is left out, never summed.

    Mutation caught: the counts kept whole (an extra sample column reaches the share)."""
    root, genes = planted(tmp_path)
    tables, problems = context.read(root, COMPARE, ["se_a", "se_e"], REGISTRY)
    assert problems == [] and tables is not None
    assert tables.counts.columns == ["gene_id", "se_a", "se_e"]
    assert context.expected(tables, genes, CONTIGS, REGISTRY)["sample"].to_list() == [
        "se_a",
        "se_e",
    ]


@pytest.mark.parametrize("dropped", ["salmon-strand_mapping_bias", "salmon-percent_mapped"])
def test_a_column_the_metrics_lack_is_a_problem_not_a_zero(tmp_path, dropped):
    """The qc metrics table without one Salmon column: a problem naming it, no table.

    Mutation caught: the column check dropped (a missing column read as null or zero)."""
    root, _ = planted(tmp_path)
    path = root / "outputs" / "02_qc" / f"{DAY}_qc_{COMPARE}_metrics_v1.parquet"
    pl.read_parquet(path).drop(dropped).write_parquet(path)
    tables, problems = context.read(root, COMPARE, list(SAMPLES), REGISTRY)
    assert tables is None and len(problems) == 1 and dropped in problems[0]


@pytest.mark.parametrize(
    ("where", "said"), [("counts", "has no ['se_e']"), ("metrics", "no qc metrics for ['se_e']")]
)
def test_a_sample_the_other_dataset_lacks_is_a_problem(tmp_path, where, said):
    """se_e taken out of one of the other dataset's tables: a problem naming it, no table.

    Mutations caught: the counts' sample check dropped (se_e's share read from no column); the
    metrics' sample check dropped (se_e left out of the table silently)."""
    root, _ = planted(tmp_path)
    folder, module = ("01_gather", "gather") if where == "counts" else ("02_qc", "qc")
    path = root / "outputs" / folder / f"{DAY}_{module}_{COMPARE}_{where}_v1.parquet"
    frame = pl.read_parquet(path)
    if where == "counts":
        frame.drop("se_e").write_parquet(path)
    else:
        frame.filter(pl.col("sample") != "se_e").write_parquet(path)
    tables, problems = context.read(root, COMPARE, list(SAMPLES), REGISTRY)
    assert tables is None and len(problems) == 1 and said in problems[0]


def test_a_missing_dataset_names_both_tables(tmp_path):
    """No compare outputs at all: both tables named, each with gather's "run the earlier module".

    Mutation caught: ``read`` returning at the first missing table (one problem)."""
    tables, problems = context.read(tmp_path, "nothing", list(SAMPLES), REGISTRY)
    assert tables is None and len(problems) == 2
    assert all("run the earlier module first" in p for p in problems)


def test_a_gene_the_annotation_lacks_stops_the_share(tmp_path):
    """The compare counts name a gene the annotation does not hold: its contig is unknown.

    Mutation caught: an unknown gene counted as off the contigs (the share silently 0.9)."""
    root, _ = planted(tmp_path)
    tables, _ = context.read(root, COMPARE, list(SAMPLES), REGISTRY)
    assert tables is not None
    genes = gene_types(GTF)  # the slice alone: OFF_CONTIG's gene is not in it
    with pytest.raises(ProcessError, match=OFF_CONTIG[0]):
        context.expected(tables, genes, CONTIGS, REGISTRY)


def test_an_unknown_gene_is_found_without_reading_the_share(tmp_path):
    """``unknown_genes`` names the gene the annotation lacks, for the door to list beside the
    other problems, and finds none in an annotation that holds every counted gene.

    Mutation caught: ``unknown_genes`` returning nothing (the door lets the gene through, and
    only ``expected`` refuses it, after every alignment)."""
    root, genes = planted(tmp_path)
    tables, _ = context.read(root, COMPARE, list(SAMPLES), REGISTRY)
    assert tables is not None
    assert context.unknown_genes(tables, genes) == []
    (found,) = context.unknown_genes(tables, gene_types(GTF))
    assert OFF_CONTIG[0] in found and "1 counted genes are not in the annotation" in found


def test_the_share_is_named_for_the_counts_gather_keeps():
    """The expected share is Salmon's length-scaled counts, the file ``ctx gather`` reads, not
    its reads: the column meaning names that file's kind (D-039).

    Mutation caught: the meaning reverted to "Salmon counts" while gather reads
    gene_counts_length_scaled."""
    from ctx.gather.rnaseq import MATRIX_FILE  # noqa: PLC0415
    from ctx.report import text  # noqa: PLC0415

    kind = Path(MATRIX_FILE).name.removeprefix("salmon.merged.").removesuffix(".tsv")
    assert kind == "gene_counts_length_scaled"
    assert kind in text.word("process", "tables.context.columns.expected_share")
    assert "length-scaled" in text.raw("process", "tables.context.page.note")
