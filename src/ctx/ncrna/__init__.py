"""M8 Non-coding RNA: biotype groups, lncRNA neighbours and the page (``ctx ncrna``).

Reads gather's counts and sample sheet, M5's results of the featured contrast and the
annotation's gene records; the TCGA half (miRNA files, targets) is not built yet.
"""


class NcrnaError(ValueError):
    """The inputs cannot make the M8 page; the message names what is missing."""
