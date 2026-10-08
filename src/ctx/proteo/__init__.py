"""M7 RNA versus protein: the protein readers, the statistics and the page of ``ctx proteo``.

Two designs (``ctx.dataset.PROTEO_DESIGNS``): one cohort measured at both layers, compared
across its samples, or RNA and protein from two studies, compared as fold changes.
"""


class ProteoError(ValueError):
    """The inputs cannot make the M7 analysis; the message names what is missing."""
