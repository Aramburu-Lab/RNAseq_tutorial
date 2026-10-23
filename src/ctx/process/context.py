"""What another dataset's Salmon run says about the same samples (``--compare-dataset``).

A partial reference cannot be judged by a genome-wide band: reads from the missing contigs have
no locus of their own there. Beside it stands the same samples' full-transcriptome Salmon run:
the share of its length-scaled gene counts (nf-core's ``gene_counts_length_scaled``, what
``ctx gather`` keeps) on genes of the genome's contigs. That share is context, not a ceiling:
reads from absent contigs can align to similar genes on the contigs present and be counted
there, so the share assigned here can exceed it, as it did in all 8 runs of the teaching subset
(D-039), since removed (D-115).
Salmon's strand call stands beside featureCounts' (D-039). Both are read from the other
dataset's gather counts and ``ctx qc`` metrics, never recomputed, and every column is checked
before it is read: a zero from a column the table does not carry is not a result.
"""

from collections.abc import Collection, Sequence
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ctx.dataset import DatasetError
from ctx.de.results import newest
from ctx.io.naming import MODULES
from ctx.process import ProcessError
from ctx.qc.registry import Registry

PERCENT, STRAND = "salmon_percent_mapped", "salmon_strand_bias"
SHOWN = 5


@dataclass(frozen=True)
class Tables:
    """The other dataset's counts and per-sample metrics, cut to this run's samples."""

    dataset: str
    counts: pl.DataFrame
    metrics: pl.DataFrame
    sources: tuple[Path, Path]


def _folder(root: Path, module: str) -> Path:
    """A module's output folder under a ctx root."""
    return root / "outputs" / f"{dict(MODULES)[module]:02d}_{module}"


def read(root: Path, dataset: str, samples: Sequence[str],
         registry: Registry) -> tuple[Tables | None, list[str]]:  # fmt: skip
    """The other dataset's newest gather counts and qc metrics, or every reason they cannot serve.

    Args:
        root: The ctx root holding both datasets' outputs.
        dataset: The other dataset's slug.
        samples: This run's samples; each must be in both tables.
        registry: Names the metrics table's Salmon columns.

    Returns:
        The tables cut to ``samples``, or None, and the problems found, all of them.
    """
    paths, problems = {}, []
    for module, what in (("gather", "counts"), ("qc", "metrics")):
        try:
            paths[what] = newest(_folder(root, module), dataset, what)
        except DatasetError as error:
            problems.append(f"--compare-dataset {dataset}: {error}")
    if problems:
        return None, problems
    counts, metrics = (pl.read_parquet(paths[w]) for w in ("counts", "metrics"))
    columns = {registry.metrics[m].column for m in (PERCENT, STRAND)}
    for name, frame, needed in (("counts", counts, {"gene_id", *samples}),
                                ("metrics", metrics, {"sample", *columns})):  # fmt: skip
        if missing := sorted(needed - set(frame.columns)):
            problems.append(f"--compare-dataset {dataset}: its {name} table "
                            f"{paths[name].name} has no {missing[:SHOWN]}"
                            f"{' and more' if len(missing) > SHOWN else ''}")  # fmt: skip
    if not problems and (absent := sorted(set(samples) - set(metrics["sample"]))):
        problems.append(f"--compare-dataset {dataset}: no qc metrics for {absent[:SHOWN]}")
    if problems:
        return None, problems
    kept = metrics.filter(pl.col("sample").is_in(samples)).select("sample", *sorted(columns))
    return Tables(dataset, counts.select("gene_id", *samples), kept,
                  (paths["counts"], paths["metrics"])), []  # fmt: skip


def unknown_genes(tables: Tables, genes: pl.DataFrame) -> list[str]:
    """The problem a counted gene the annotation lacks makes: its contig is unknown.

    Checked at the door (``run.resolve``), with the other input problems, before anything is
    aligned; ``expected`` refuses the same genes, so neither can read past them.

    Args:
        tables: ``read``'s tables.
        genes: The full annotation's gene records (``checkpoint.gene_types``): gene_id, chrom.

    Returns:
        One problem naming the first genes, or none.
    """
    known = set(genes["gene_id"])
    if unknown := [g for g in tables.counts["gene_id"] if g not in known]:
        return [f"--compare-dataset {tables.dataset}: {len(unknown)} counted genes are not in "
                f"the annotation, so their contig is unknown: {unknown[:SHOWN]}"]  # fmt: skip
    return []


def expected(tables: Tables, genes: pl.DataFrame, contigs: Collection[str],
             registry: Registry) -> pl.DataFrame:  # fmt: skip
    """Per sample: Salmon's share of counts on the genome's contigs, its mapping and strand.

    Args:
        tables: ``read``'s tables.
        genes: The full annotation's gene records (``checkpoint.gene_types``): gene_id, chrom.
        contigs: The genome's contig names.
        registry: Its ``salmon_strand_bias`` bands call Salmon's strand.

    Returns:
        sample, salmon_percent_mapped, expected_share (length-scaled counts on genes of the
        contigs over all of them), expected_read_share (that share times Salmon's mapped
        fraction: a share of every read Salmon saw, if length scaling moved no count between
        contigs), salmon_strand_bias and salmon_call.

    Raises:
        ProcessError: If a counted gene is not in the annotation (``unknown_genes``).
    """
    if problems := unknown_genes(tables, genes):
        raise ProcessError(problems[0])
    on = dict(genes.select("gene_id", pl.col("chrom").is_in(list(contigs))).rows())
    counts = tables.counts.with_columns(on=pl.col("gene_id").replace_strict(on))
    samples = counts.columns[1:-1]
    totals = counts.select(pl.col(samples).sum())
    inside = counts.filter("on").select(pl.col(samples).sum())
    share = {s: float(inside[s][0] / totals[s][0]) if totals[s][0] else None for s in samples}
    percent = registry.metrics[PERCENT].column
    bias = registry.metrics[STRAND].column
    metric = registry.metrics[STRAND]
    frame = tables.metrics.select(
        "sample",
        salmon_percent_mapped=pl.col(percent).cast(pl.Float64),
        salmon_strand_bias=pl.col(bias).cast(pl.Float64),
    ).with_columns(expected_share=pl.col("sample").replace_strict(share, return_dtype=pl.Float64))
    calls = [metric.strandedness(b) if b is not None else "undetermined"
             for b in frame["salmon_strand_bias"]]  # fmt: skip
    return frame.with_columns(
        expected_read_share=pl.col("expected_share") * pl.col("salmon_percent_mapped") / 100,
        salmon_call=pl.Series(calls, dtype=pl.String),
    ).select("sample", "salmon_percent_mapped", "expected_share", "expected_read_share",
             "salmon_strand_bias", "salmon_call")  # fmt: skip
