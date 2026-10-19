"""M4's figure kinds: drawn from the tables ``ctx quant`` wrote.

Two of its six have a spec. The correlation heatmap reads the integer count matrix, which M4
does not write; the detection curves, the units scatter and the isoform bars share their stem
with the table M4 writes after them, whose CSV and README replace the figure's on disk.
"""

from ctx.figures.kind import DATASET_FILE, Kind, OneContrast, Plain, featured, once

DRAW = "ctx.figures.draw_m4"
KINDS = (
    Kind("quant", "fates", Plain, "ctx.quant.figures:fate_bars", f"{DRAW}:fates", once,
         lambda _: "fates", lambda _, __: {"fates": ("quant", "fates"),
                                           "figure_sheet": ("quant", "figure_sheet"),
                                           "dataset_file": DATASET_FILE}),
    Kind("quant", "estimates", OneContrast, "ctx.quant.figures:estimates_scatter",
         f"{DRAW}:estimates", featured, lambda _: "estimates",
         lambda _, __: {"fold_changes": ("quant", "fold_changes"),
                        "contrast_samples": ("quant", "contrast_samples"),
                        "dataset_file": DATASET_FILE}),
)  # fmt: skip
