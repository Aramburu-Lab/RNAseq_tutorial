"""Read one nf-core/rnaseq outdir, Salmon-only or aligned with STAR, in the pinned layout.

File names are from the 3.27.0 salmon test snapshot (tests/salmon.nf.test.snap) and are
confirmed by tests/fixtures/rnaseq, a real Salmon-only run, and tests/fixtures/rnaseq_star, a
real run with `aligner: star_salmon` (D-116). Nothing here writes into the outdir.
"""

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

# The one definition. An upgrade is one deliberate change: this, a new fixture, and the
# cluster commands in TODO.md, which a test holds to this value.
RNASEQ_VERSION = "3.27.0"
VERSIONS_FILE = "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml"
MATRIX_FILE = "salmon/salmon.merged.gene_counts_length_scaled.tsv"
# What ctx reads. pydeseq2 needs integers and takes one offset per sample, so the gene matrix
# is the length-scaled one (D-007).
READ_FILES = (MATRIX_FILE, VERSIONS_FILE, "pipeline_info/params_*.json")
# What leaves the cluster, the fixture and the bundle: READ_FILES plus per-sample library
# format and mapping facts, the Trim Galore JSON reports (the read fates of `ctx qc`, D-022; JSON
# because MultiQC 1.35 reads Trim Galore 2.x from nothing else, D-025) and the trace (measured
# time and memory). Transcript matrices, .rds files, MultiQC and the HTML
# reports stay on the cluster.
COPIED_FILES = (
    *READ_FILES,
    "salmon/*/lib_format_counts.json",
    "salmon/*/aux_info/meta_info.json",
    "trimgalore/*_trimming_report.json",
    "pipeline_info/execution_trace_*.txt",
)
# What an outdir aligned with STAR adds (`aligner: star_salmon`, D-113): Salmon on STAR's
# alignment and the tables nf-core's own MultiQC wrote, which M2's read stage and M3 read in
# place of running any tool again (D-116). Never a BAM: M3 reads no alignment itself.
ALIGNED_MATRIX = "star_salmon/salmon.merged.gene_counts_length_scaled.tsv"
MULTIQC_DATA = "multiqc/star_salmon/multiqc_report_data"
MULTIQC_TABLES = (
    # The MultiQC version, which names the tables' columns.
    "multiqc.log",
    # M2: Trim Galore, Salmon alone, FastQC on the raw reads and its plots.
    "multiqc_general_stats.txt", "multiqc_trim_galore.txt", "multiqc_salmon.txt",
    "multiqc_fastqc_fastqc_raw.txt", "fastqc_raw_per_base_sequence_quality_plot.txt",
    "fastqc_raw_adapter_content_plot.txt", "fastqc_raw_per_base_n_content_plot.txt",
    "fastqc_raw_per_sequence_gc_content_plot_Percentages.txt",
    "fastqc_raw_overrepresented_sequences_plot.txt",
    # M3: STAR, samtools flagstat, RSeQC, Qualimap, dupRadar, Preseq, featureCounts biotypes.
    "multiqc_star.txt", "multiqc_samtools_flagstat.txt", "rseqc_infer_experiment_plot.txt",
    "multiqc_strand_check_summary_table.txt", "rseqc_read_distribution_plot.txt",
    "multiqc_rseqc_junction_annotation.txt", "qualimap_gene_coverage_profile_Normalised.txt",
    "multiqc_dupradar.txt", "preseq_complexity_plot_molecules.txt",
    "multiqc_featurecounts_biotype_plot.txt",
)  # fmt: skip
ALIGNED_FILES = (ALIGNED_MATRIX, *(f"{MULTIQC_DATA}/{name}" for name in MULTIQC_TABLES))
# Copied when present, never required: MultiQC writes this table only when some sample has an
# overrepresented sequence (BRCA-immune-cells-macrophages, 2026-10-24, has none), and its reader,
# ctx.qc.fastqc, draws nothing without it.
OPTIONAL_FILES = (f"{MULTIQC_DATA}/fastqc_raw_overrepresented_sequences_plot.txt",)
ALIGNER = "star_salmon"
# What leaves Kebnekaise, which aligns (D-113): the Salmon-only list and the aligned files.
CLUSTER_FILES = (*COPIED_FILES, *ALIGNED_FILES)
# Every params file of one outdir must agree on these; -resume adds a file per launch.
AGREE = ("aligner", "pseudo_aligner", "skip_alignment", "gencode", "fasta", "gtf")
INT32_MAX = 2**31 - 1


class RnaseqError(ValueError):
    """An rnaseq outdir that cannot be used as it stands; the message lists every offender."""


@dataclass(frozen=True)
class Counts:
    """The gene matrix of one rnaseq outdir, checked and made integer.

    Attributes:
        matrix: ``gene_id`` then one integer column per sample, Int32 where every value fits.
        genes: ``gene_id`` and ``gene_name``.
        fractional: Values that were not whole numbers before rounding.
        provenance: Release, Salmon version, reference files, and the SHA-256 of each file read.
    """

    matrix: pl.DataFrame
    genes: pl.DataFrame
    fractional: int
    provenance: dict[str, str]


def _sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _versions(outdir: Path) -> dict[str, str]:
    """Refuse any release but the pinned one; record it, Nextflow and Salmon."""
    path = outdir / VERSIONS_FILE
    if not path.is_file():
        raise RnaseqError(f"{path} is missing: not an nf-core/rnaseq {RNASEQ_VERSION} outdir")
    text = path.read_text()
    found = re.search(r"nf-core/rnaseq: v(\d+\.\d+\.\d+)(?:-g([0-9a-f]+))?", text)
    if found is None or found.group(1) != RNASEQ_VERSION:
        seen = found.group(1) if found else "unreadable"
        raise RnaseqError(f"nf-core/rnaseq {seen} wrote {outdir}; ctx reads {RNASEQ_VERSION} only")
    salmon = sorted(set(re.findall(r"^\s+salmon: (\S+)$", text, re.MULTILINE)))
    nextflow = re.search(r"Nextflow: (\S+)", text)
    return {
        "nf-core/rnaseq": found.group(1),
        "nf-core/rnaseq commit": found.group(2) or "not recorded",
        "salmon": ", ".join(salmon) or "not recorded",
        "nextflow (rnaseq)": nextflow.group(1) if nextflow else "not recorded",
        f"sha256:{VERSIONS_FILE}": _sha256(path),
    }


def _params(outdir: Path) -> dict[str, str]:
    """The params files must agree on the settings that shape the matrix."""
    found = sorted((outdir / "pipeline_info").glob("params_*.json"))
    if not found:
        raise RnaseqError(f"{outdir}: no pipeline_info/params_*.json")
    loaded = {p.name: json.loads(p.read_text()) for p in found}
    first = loaded[found[0].name]
    differ = [k for k in AGREE if any(p.get(k) != first.get(k) for p in loaded.values())]
    if differ:
        raise RnaseqError(f"{outdir}: the {len(found)} params files disagree on {differ}")
    if first.get("pseudo_aligner") != "salmon":
        raise RnaseqError(f"{outdir}: pseudo_aligner {first.get('pseudo_aligner')!r}; ctx reads "
                          f"Salmon's {MATRIX_FILE}")  # fmt: skip
    if first.get("skip_alignment") is False:
        _aligned(outdir, first.get("aligner"))
    return {
        **{f"params {k}": str(first.get(k)) for k in AGREE},
        **{
            f"sha256:pipeline_info/{name}": _sha256(outdir / "pipeline_info" / name)
            for name in loaded
        },  # fmt: skip
    }


def _aligned(outdir: Path, aligner: object) -> None:
    """An outdir aligned with STAR holds every file ``ALIGNED_FILES`` names, or is refused.

    A copy that missed one stops gather, before M2 or M3 would stop on it.
    """
    if aligner != ALIGNER:
        raise RnaseqError(f"{outdir}: aligner {aligner!r}; ctx reads {ALIGNER}'s folders only")
    required = [p for p in ALIGNED_FILES if p not in OPTIONAL_FILES]
    if missing := [p for p in required if not (outdir / p).is_file()]:
        raise RnaseqError(f"{outdir} was aligned ({ALIGNER}) but holds no {missing}; copy it "
                          "with the include list of rnaseq.ALIGNED_FILES")  # fmt: skip


def first_params(outdir: Path) -> dict[str, object]:
    """The outdir's first params file, by name; ``_params`` holds the rest to it on ``AGREE``.

    Raises:
        RnaseqError: If it has no params file.
    """
    found = sorted((outdir / "pipeline_info").glob("params_*.json"))
    if not found:
        raise RnaseqError(f"{outdir}: no pipeline_info/params_*.json")
    loaded: dict[str, object] = json.loads(found[0].read_text())
    return loaded


def provenance(outdir: Path) -> dict[str, str]:
    """The release's and the params files' checks, and what they record (``read``'s first two).

    Raises:
        RnaseqError: If the release is not the pinned one or the params files disagree, or an
            aligned outdir lacks a file ``ALIGNED_FILES`` names.
    """
    return _versions(outdir) | _params(outdir)


def aligned(outdir: Path) -> bool:
    """Whether an outdir was aligned with STAR: its params say ``skip_alignment`` false.

    Raises:
        RnaseqError: If it has no params file, or its params files disagree.
    """
    _params(outdir)
    return first_params(outdir).get("skip_alignment") is False


def to_counts(values: np.ndarray) -> tuple[np.ndarray, int]:
    """Round half to even, then choose the narrowest integer type the range allows.

    Args:
        values: Length-scaled counts, genes by samples.

    Returns:
        The integer counts and how many values were not whole before rounding.

    Raises:
        RnaseqError: On a missing or negative value.
    """
    if not np.isfinite(values).all() or (values < 0).any():
        raise RnaseqError("the matrix holds missing, infinite or negative values")
    rounded = np.rint(values)
    fractional = int((rounded != values).sum())
    kind = np.int32 if rounded.max(initial=0) <= INT32_MAX else np.int64
    return rounded.astype(kind), fractional


def read(outdir: Path, samples: set[str]) -> Counts:
    """Check and read one rnaseq outdir, in a fixed order.

    Order: release; params agreement; a lazy scan of the length-scaled matrix; its samples
    against the sheet's; rounding, then the range check, then the integer type.

    Args:
        outdir: An nf-core/rnaseq outdir, as the pipeline wrote it.
        samples: The ``sample`` values of the harmonised sheet.

    Returns:
        The checked counts.

    Raises:
        RnaseqError: If the outdir cannot be used as it stands.
    """
    provenance = _versions(outdir) | _params(outdir)
    path = outdir / MATRIX_FILE
    if not path.is_file():
        raise RnaseqError(f"{path} is missing")
    header = pl.scan_csv(path, separator="\t", n_rows=0).collect_schema().names()
    in_matrix = header[2:]
    problems = [
        f"sample {s}: in the sheet, not in the matrix" for s in sorted(samples - set(in_matrix))
    ]
    problems += [
        f"sample {s}: in the matrix, not in the sheet" for s in sorted(set(in_matrix) - samples)
    ]
    if header[:2] != ["gene_id", "gene_name"]:
        problems.insert(0, f"{path.name}: header starts {header[:2]}, not gene_id, gene_name")
    if len(in_matrix) != len(set(in_matrix)):
        problems.append(f"{path.name}: a sample column appears twice")
    if problems:
        raise RnaseqError(f"{outdir}:\n  " + "\n  ".join(problems))
    schema = {"gene_id": pl.String, "gene_name": pl.String, **dict.fromkeys(in_matrix, pl.Float64)}
    frame = pl.scan_csv(path, separator="\t", schema=schema).collect()
    if frame["gene_id"].is_duplicated().any():
        raise RnaseqError(f"{path.name}: a gene_id appears twice")
    integers, fractional = to_counts(frame.select(in_matrix).to_numpy())
    matrix = pl.DataFrame(integers, schema=in_matrix, orient="row")
    matrix = matrix.insert_column(0, frame["gene_id"])
    provenance[f"sha256:{MATRIX_FILE}"] = _sha256(path)
    if (outdir / ALIGNED_MATRIX).is_file():
        provenance[f"sha256:{ALIGNED_MATRIX}"] = _sha256(outdir / ALIGNED_MATRIX)
    return Counts(matrix, frame.select("gene_id", "gene_name"), fractional, provenance)
