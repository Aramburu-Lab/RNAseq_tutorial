"""M3: trim, align and count reads on a partial reference (`ctx process`).

Dependency-free, so the command line reads its choices without importing the analysis stack
(D-034's rule for `ctx qc`).
"""

from enum import StrEnum


class ProcessError(ValueError):
    """An input `ctx process` cannot use, or a check its outputs failed; the message says which."""


class Aligner(StrEnum):
    """The aligner: STAR on Linux, HISAT2 on Apple silicon.

    The pinned osx-arm64 STAR builds indices but reads no reads (D-036).
    """

    STAR = "star"
    HISAT2 = "hisat2"


class Strandedness(StrEnum):
    """The library's strand: called from the reads (auto) or given."""

    AUTO = "auto"
    UNSTRANDED = "unstranded"
    FORWARD = "forward"
    REVERSE = "reverse"


# The one map of a strand to STAR's ReadsPerGene column and featureCounts' -s: column 2
# counts either strand, column 3 reads whose read 1 lies on the gene's strand, column 4 reads
# whose read 2 does (STAR manual 2.7.11b, section 8); -s 0, 1, 2 likewise (Subread Users Guide,
# 6.2.10).
STRAND = {
    Strandedness.UNSTRANDED: (2, 0),
    Strandedness.FORWARD: (3, 1),
    Strandedness.REVERSE: (4, 2),
}
