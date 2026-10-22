"""M8 Non-coding RNA: biotype groups, lncRNA neighbours and the page (``ctx ncrna``).

Reads gather's counts and sample sheet, M5's results of the featured contrast and the
annotation's gene records; with a small-RNA library declared, the miRNA files and targets too.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ctx.dataset import Contrast, Dataset


class NcrnaError(ValueError):
    """The inputs cannot make the M8 page; the message names what is missing."""


def small_rna(spec: "Dataset") -> "Contrast | None":
    """The miRNA contrast of the dataset's small-RNA libraries, None where it declares none.

    The one door to every miRNA panel, table and claim of M8: standard RNA-seq holds no
    mature miRNAs (CONVENTIONS.md, Statistics rules), so a file without ``[ncrna.mirna]``
    gets the biotype half alone, whatever its route.
    """
    return spec.ncrna.mirna if spec.ncrna is not None else None
