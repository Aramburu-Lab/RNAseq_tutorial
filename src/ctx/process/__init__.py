"""M3 (`ctx process`): nf-core/rnaseq's alignment read from its outputs, or reads aligned here.

Two routes. With ``--rnaseq-dir`` the page reads nf-core/rnaseq 3.27.0's STAR alignment and QC
from its own tables (``ctx.process.nfcore``, D-116); with ``--genome`` the laptop route trims,
aligns and counts on a partial reference here. Dependency-free, so the command line reads its
choices without importing the analysis stack (D-034's rule for `ctx qc`).
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


# The featured sample's extra alignments a reader can ask for (--settings); under HISAT2 an
# unspliced alignment is added (D-036: STAR has no switch that turns splicing off).
SETTINGS = ("cut", "untrimmed")
DEFAULT_SETTINGS = "cut"
# --multimap-max: STAR's --outFilterMultimapNmax must stay within --winAnchorMultimapNmax,
# whose default is 50 (STAR manual 2.7.11b, 4.1).
MULTIMAP_RANGE = (1, 50)
# Each aligner's outcomes for a read, the order the mapping figure stacks them in.
CATEGORIES = {
    Aligner.STAR: ("unique", "multi", "too_many", "mismatches", "too_short", "other"),
    Aligner.HISAT2: ("unique", "multi", "unmapped"),
}

# The one map of a strand to STAR's ReadsPerGene column and featureCounts' -s: column 2
# counts either strand, column 3 reads whose read 1 lies on the gene's strand, column 4 reads
# whose read 2 does (STAR manual 2.7.11b, section 8); -s 0, 1, 2 likewise (Subread Users Guide,
# 6.2.10).
STRAND = {
    Strandedness.UNSTRANDED: (2, 0),
    Strandedness.FORWARD: (3, 1),
    Strandedness.REVERSE: (4, 2),
}
