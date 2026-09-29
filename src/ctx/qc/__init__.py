"""M2: read QC and the checkpoint on counts (`ctx qc`)."""


class QcError(ValueError):
    """A QC input that cannot be used; the message names every problem found."""
