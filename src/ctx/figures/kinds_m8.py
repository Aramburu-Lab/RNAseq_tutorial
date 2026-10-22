"""M8's figure kinds: drawn from the tables ``ctx ncrna`` wrote, gather's, M5's and the targets.

The Step-1 half reads the run's options and the annotation's gene records M8 wrote, with
gather's sheet and counts and M5's results; the miRNA half of a GDC dataset reads M8's miRNA
fit, gather's miRNA tables and the four target files ``ctx setup`` fetched. Nothing is fitted
again: the miRNA volcano reads the fit M8 wrote.
"""

from typing import Any

from ctx.dataset import Dataset
from ctx.figures.kind import DATASET_FILE, Inputs, Kind, Plain, m5_results, resource

SHEET, COUNTS = ("gather", "samplesheet"), ("gather", "counts")
GENES, OPTIONS = ("ncrna", "gene_biotypes"), ("ncrna", "figure_options")


def _target_files() -> tuple[str, ...]:
    """The target resources' names, read when asked: their module reads the files."""
    from ctx.ncrna.targets import RESOURCES  # noqa: PLC0415

    return RESOURCES


def _results(spec: Dataset) -> tuple[str, str]:
    """M5's results of the contrast ``[ncrna]`` features."""
    return m5_results(spec.ncrna.featured if spec.ncrna else "")


def _with(*names: str, targets: bool = False) -> Any:
    """These inputs by name, the featured results and the dataset file, the targets if asked."""
    tables = {"samplesheet": SHEET, "counts": COUNTS, "gene_biotypes": GENES,
              "figure_options": OPTIONS, "mirna_results": ("ncrna", "mirna_results"),
              "mirna_samplesheet": ("gather", "mirna_samplesheet"),
              "mirna_mature_counts": ("gather", "mirna_mature_counts")}  # fmt: skip

    def inputs(_: Any, spec: Dataset) -> Inputs:
        """The kind's inputs for one dataset file."""
        found = {n: _results(spec) if n == "results" else tables[n] for n in names}
        found |= {n: resource(n) for n in _target_files()} if targets else {}
        return {**found, "dataset_file": DATASET_FILE}

    return inputs


def _step1(spec: Dataset) -> list[Plain]:
    """The one figure of a kind every dataset file with ``[ncrna]`` makes."""
    return [Plain(spec.slug)] if spec.ncrna else []


def _mirna(spec: Dataset) -> list[Plain]:
    """The one figure of a kind of the miRNA half, which a GDC dataset alone draws."""
    return [Plain(spec.slug)] if spec.ncrna and spec.ncrna.mirna and spec.gdc else []


DRAW = "ctx.figures.draw_m8"
LIBRARY = _with("samplesheet", "counts", "gene_biotypes", "figure_options", "results")
KINDS = (
    Kind("ncrna", "biotype_shares", Plain, "ctx.ncrna.page:library", f"{DRAW}:biotype_shares",
         _step1, lambda _: "biotype_shares", LIBRARY),
    Kind("ncrna", "filter_by_biotype", Plain, "ctx.ncrna.page:library",
         f"{DRAW}:filter_by_biotype", _step1, lambda _: "filter_by_biotype", LIBRARY),
    Kind("ncrna", "volcano_by_biotype", Plain, "ctx.ncrna.page:de", f"{DRAW}:volcano_by_biotype",
         _step1, lambda _: "volcano_by_biotype",
         _with("samplesheet", "gene_biotypes", "figure_options", "results")),
    Kind("ncrna", "lncrna_neighbour", Plain, "ctx.ncrna.page:lncrna",
         f"{DRAW}:lncrna_neighbour", _step1, lambda _: "lncrna_neighbour",
         _with("samplesheet", "counts", "gene_biotypes", "results")),
    Kind("ncrna", "mirna_volcano", Plain, "ctx.ncrna.mirna_page:mirna_volcano",
         f"{DRAW}:mirna_volcano", _mirna, lambda _: "mirna_volcano",
         _with("mirna_results", "mirna_samplesheet")),
    Kind("ncrna", "target_shift", Plain, "ctx.ncrna.mirna_page:shift", f"{DRAW}:target_shift",
         _mirna, lambda _: "target_shift",
         _with("gene_biotypes", "results", "mirna_results", targets=True)),
    Kind("ncrna", "mirna_target_pairs", Plain, "ctx.ncrna.mirna_page:pairs",
         f"{DRAW}:mirna_target_pairs", _mirna, lambda _: "mirna_target_pairs",
         _with("counts", "results", "mirna_samplesheet", "mirna_mature_counts", targets=True)),
    Kind("ncrna", "library_types", Plain, "ctx.ncrna.libraries_page:group",
         f"{DRAW}:library_types", _mirna, lambda _: "library_types",
         _with("samplesheet", "counts", "gene_biotypes")),
)  # fmt: skip
