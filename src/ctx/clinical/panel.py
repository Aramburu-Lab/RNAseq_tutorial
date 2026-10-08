"""The genes a page carries for live scoring, and the block that carries them.

Decision 2 of the M10 plan (2026-10-08): until Step-3 DE, the panel is the most variable of the
expressed protein-coding genes, as many as fit ``PANEL_BYTES``, the build list always in. Each
gene's stored int16 values over the patients go into the page as two byte planes, all high bytes
then all low bytes, gzip-compressed and in base64: the page inflates them with the browser's
DecompressionStream and ``ctxSurv.planes`` rebuilds the int16 values, the build's own.
"""

import base64
import gzip
from dataclasses import dataclass

import numpy as np

# The block's share of the 15,000,000-byte page cap. On the breast cohort (2026-10-08) the page was
# 14,509,635 bytes at this budget, 8,822,282 of them outside the panel's settings: plotly.js
# 4,815,817, the export's libraries 526,581 and fonts 219,314, the figures and their rows the
# rest. At 6,500,000 it was 15,593,683 and refused; the GS plan's 6,729,268 for the block
# assumed the preview's lighter page.
PANEL_BYTES = 5_600_000
# gzip's mtime is fixed so the same values give the same bytes.
MTIME = 0


@dataclass(frozen=True)
class Panel:
    """The genes carried and their encoded values.

    Attributes:
        rows: Row of each carried gene in the expressed genes, in carried order.
        block: The base64 text of the gzip-compressed byte planes.
        genes: Genes carried.
        patients: Patients per gene.
    """

    rows: np.ndarray
    block: str
    genes: int
    patients: int


def encode(values: np.ndarray) -> str:
    """int16 genes by patients as two byte planes, gzip-compressed, in base64.

    Each value is written as an unsigned 16-bit two's complement: the high bytes of every
    value, then the low bytes, gene after gene, so each plane compresses on its own.
    """
    unsigned = np.asarray(values, np.int16).view(np.uint16).ravel()
    planes = np.concatenate([(unsigned >> 8).astype(np.uint8), (unsigned & 0xFF).astype(np.uint8)])
    return base64.b64encode(gzip.compress(planes.tobytes(), mtime=MTIME)).decode()


def decode(block: str, genes: int, patients: int) -> np.ndarray:
    """``encode``'s inverse, as the page decodes it."""
    planes = np.frombuffer(gzip.decompress(base64.b64decode(block)), np.uint8).astype(np.uint16)
    half = genes * patients
    return ((planes[:half] << 8) | planes[half:]).view(np.int16).reshape(genes, patients)


def choose(values: np.ndarray, coding: np.ndarray, keep: np.ndarray,
           budget: int = PANEL_BYTES) -> Panel:  # fmt: skip
    """The most variable coding genes whose block fits the budget, the ``keep`` rows first.

    Args:
        values: Stored int16 values of the expressed genes, genes by patients.
        coding: Per gene, whether it is protein-coding.
        keep: Rows always carried: the build list.
        budget: Bytes the base64 block may take.

    Returns:
        The panel: by bisection, the largest count of genes whose block fits.
    """
    variance = (np.asarray(values, float) / 1000).var(axis=1)
    ranked = [int(i) for i in np.argsort(-variance, kind="stable") if coding[i]]
    first = [int(i) for i in keep]
    order = np.array([*first, *(i for i in ranked if i not in set(first))])

    def fits(n: int) -> bool:
        """Whether the first n genes' block fits the budget."""
        return len(encode(values[order[:n]])) <= budget

    # Bisect around the count the first 500 genes' bytes per gene predict, each bound checked.
    probe = min(len(order), 500)
    guess = int(budget * probe / max(len(encode(values[order[:probe]])), 1))
    low = max(len(first), min(len(order), int(guess * 0.8)))
    high = max(low, min(len(order), int(guess * 1.25) + 1))
    if not fits(low):
        low = len(first)
    if fits(high):
        low, high = high, len(order)
    while low < high:
        mid = (low + high + 1) // 2
        if len(encode(values[order[:mid]])) <= budget:
            low = mid
        else:
            high = mid - 1
    rows = order[:low]
    return Panel(rows, encode(values[rows]), len(rows), values.shape[1])
