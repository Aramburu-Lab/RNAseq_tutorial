"""M9 Gene cards: a gene's result in this dataset beside GTEx, HPA and Ensembl reference data.

GTEx V11 is read from its release files (``gtex``), HPA and Ensembl from per-gene responses
cached by ``cache``; ``summary`` holds the one quantile rule and the one tau.
"""


class GenecardError(ValueError):
    """A gene card cannot be built; the message names the gene, file or response."""
