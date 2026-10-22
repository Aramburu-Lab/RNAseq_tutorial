"""What `ctx process` was asked, what it found, and what each stage hands the next.

The run's dataclasses, apart from ``ctx.process.run``, so the two routes (aligning here, or
reading a cluster's alignment through ``ctx.process.cluster``) share one set of types.
"""

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from ctx.dataset import Dataset
from ctx.process import (
    DEFAULT_SETTINGS,
    Aligner,
    Strandedness,
    align,
    cluster,
    context,
    count,
    reference,
)
from ctx.qc import Fastp
from ctx.qc import multiqc_io as mqc_io
from ctx.qc import registry as reg


@dataclass(frozen=True)
class Options:
    """What `ctx process` was asked to do, beside the dataset.

    Attributes:
        dataset_file: The dataset file; its ``annotation_gtf`` resource is the annotation.
        genome: The genome FASTA, plain or gzipped.
        aligner: STAR or HISAT2 (D-036).
        strandedness: The strand to count at; auto calls it per sample from the reads.
        multimap: STAR's --outFilterMultimapNmax; None keeps STAR's 10.
        settings: The featured sample's extra alignments, as tokens of ``SETTINGS``.
        featured: The featured sample; None takes the sheet's first.
        reads_per_run: Reads kept per run before fastp; None keeps every read.
        compare: A dataset whose Salmon counts and metrics cover the same samples.
        fetchngs_dir: The outdir the sheet's FASTQ paths are relative to; None is
            ``<root>/data/raw/<dataset>``.
        registry: The QC registry; None resolves it.
        fastp: fastp's filter, M2's own, so M2 and M3 trim alike.
        align_dir: A cluster's alignment as ED.6 block D copies it, read instead of aligning
            here; never the default (``ctx.process.cluster``). ``genome`` is then optional.
    """

    dataset_file: Path
    genome: Path | None
    aligner: Aligner = Aligner.STAR
    strandedness: Strandedness = Strandedness.AUTO
    multimap: int | None = None
    settings: tuple[str, ...] = (DEFAULT_SETTINGS,)
    featured: str | None = None
    reads_per_run: int | None = None
    compare: str | None = None
    fetchngs_dir: Path | None = None
    registry: Path | None = None
    fastp: Fastp = field(default_factory=Fastp)
    align_dir: Path | None = None


@dataclass(frozen=True)
class Sheet:
    """The gather sheet, where it is, its readable runs' files and the runs left unread."""

    sheet: pl.DataFrame
    path: Path
    files: pl.DataFrame
    skipped: pl.DataFrame


@dataclass(frozen=True)
class Inputs:
    """Every input, checked: the dataset, its sheet and reads, the annotation and the reference.

    Attributes:
        genes: The whole annotation's gene records (``checkpoint.gene_types``), read once.
        genome_name: The genome FASTA's file name, here or on the cluster.
        source: The cluster's alignment, or None when this run aligns.
    """

    spec: Dataset
    sheet: pl.DataFrame
    sheet_path: Path
    files: pl.DataFrame
    skipped: pl.DataFrame
    gtf: Path
    genes: pl.DataFrame
    registry: reg.Registry
    featured: str
    compare: context.Tables | None
    reference: reference.Reference
    genome_name: str
    source: cluster.Source | None = None


@dataclass(frozen=True)
class Built:
    """A reference and its index, and how long the index took (0 when reused)."""

    reference: reference.Reference
    index: Path
    built: bool
    seconds: float


@dataclass(frozen=True)
class Aligned:
    """The trimmed runs, the main alignment, the featured sample's settings and their reports."""

    trimmed: align.Trimmed
    main: align.Alignment
    featured: dict[str, align.Alignment]
    reports: dict[str, mqc_io.MultiqcRun]


@dataclass(frozen=True)
class Counted:
    """Every count folder, its MultiQC runs, the strand calls and the strand each sample used."""

    plain: count.Counts
    variants: count.Counts
    cut: count.Counts | None
    runs: dict[str, mqc_io.MultiqcRun]
    counted: pl.DataFrame
    calls: pl.DataFrame
    used: dict[str, Strandedness]


# What counts an alignment at some settings: featureCounts here, or a cluster's files joined.
# The str names the reference counted on: "main", or "cut" for the featured sample's half.
Counter = Callable[[str, align.Alignment, Mapping[str, count.FcOptions]], count.Counts]
