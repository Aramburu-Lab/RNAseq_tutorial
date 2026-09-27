"""Read one nf-core/rnaseq outdir (Salmon-only), in the layout of the pinned release.

File names are from the 3.27.0 salmon test snapshot (tests/salmon.nf.test.snap) and are
confirmed by tests/fixtures/rnaseq, a real run. Nothing here writes into the outdir.
"""

# The one definition. An upgrade is one deliberate change: this, a new fixture, and the
# cluster commands in TODO.md, which a test holds to this value.
RNASEQ_VERSION = "3.27.0"
VERSIONS_FILE = "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml"
# What ctx reads. pydeseq2 needs integers and takes one offset per sample, so the gene matrix
# is the length-scaled one (D-007).
READ_FILES = (
    "salmon/salmon.merged.gene_counts_length_scaled.tsv",
    VERSIONS_FILE,
    "pipeline_info/params_*.json",
)
# What leaves the cluster, the fixture and the bundle: READ_FILES plus per-sample library
# format and mapping facts (Phase 3 QC) and the trace (measured time and memory). Transcript
# matrices, .rds files, MultiQC and the HTML reports stay on the cluster.
COPIED_FILES = (
    *READ_FILES,
    "salmon/*/lib_format_counts.json",
    "salmon/*/aux_info/meta_info.json",
    "pipeline_info/execution_trace_*.txt",
)
