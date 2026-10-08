"""M10 Clinical and survival: the patient table, survival statistics, the score and the page.

Reads gather's GDC tables only (``outputs/01_gather/``), files through ``gdc.chosen()``.
"""


class ClinicalError(ValueError):
    """Gather's tables cannot make the clinical table; the message names what is missing."""
