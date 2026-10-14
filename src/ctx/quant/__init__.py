"""M4 (`ctx quant`): what Salmon did with each sample's reads, and how genes get counts.

Dependency-free, so the command line reads its choices without importing the analysis stack
(D-034's rule for `ctx qc`).
"""

from dataclasses import dataclass
from enum import StrEnum

# quant.sf prints NumReads to three decimals: a count of reads carries no finer precision.
NUM_READS_DECIMALS = 3


class QuantError(ValueError):
    """An input `ctx quant` cannot use, or a check its outputs failed; the message says which."""


class CountsFrom(StrEnum):
    """How gene counts are made from transcript estimates: tximport's countsFromAbundance."""

    NO = "no"
    SCALED_TPM = "scaled_tpm"
    LENGTH_SCALED_TPM = "length_scaled_tpm"


@dataclass(frozen=True)
class SalmonOptions:
    """One quantification setting: Salmon's library type, its bias models and the decoys.

    Attributes:
        libtype: Salmon's ``-l``; A detects it per sample.
        gc_bias: Salmon's ``--gcBias``.
        seq_bias: Salmon's ``--seqBias``.
        decoys: Index the genome's contigs as decoys, as nf-core's index does.
    """

    libtype: str = "A"
    gc_bias: bool = False
    seq_bias: bool = False
    decoys: bool = True
