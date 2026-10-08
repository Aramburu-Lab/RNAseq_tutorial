"""M4: what Salmon did with each sample's reads, and how transcripts become gene counts (`ctx quant`).

Dependency-free, so the command line reads its choices without importing the analysis stack
(D-034's rule for `ctx qc`).
"""

from enum import StrEnum


class QuantError(ValueError):
    """An input `ctx quant` cannot use, or a check its outputs failed; the message says which."""


class CountsFrom(StrEnum):
    """How gene counts are made from transcript estimates: tximport's countsFromAbundance."""

    NO = "no"
    SCALED_TPM = "scaled_tpm"
    LENGTH_SCALED_TPM = "length_scaled_tpm"
