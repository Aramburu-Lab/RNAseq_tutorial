"""A box figure keeps its feature labels apart at a phone's width (M10, the box labels at 375 px).

Measured on the TCGA kidney page at 375 px in Chromium (2026-10-09): four genes in 63 px slots,
NDUFA4L2's tick label 101 px wide over VEGFA's, the p labels over each other.
"""

import polars as pl

from ctx.viz.box import box_data, least_width
from ctx.viz.metrics import text_width
from ctx.viz.style import BOX_GAP_PX, PHONE_PLOT_PX, SIZES, on_page

GENES = ("CA9", "NDUFA4L2", "VEGFA", "MKI67")


def kidney_like() -> pl.DataFrame:
    """Four genes, two groups of six, values apart enough for small p labels."""
    rows = [{"feature": g, "group": grp, "sample": f"{grp}{i}",
             "value": i + (3.0 if grp == "b" else 0)}
            for g in GENES for grp in ("a", "b") for i in range(6)]  # fmt: skip
    return pl.DataFrame(rows)


def test_the_box_is_never_drawn_narrower_than_its_labels_need():
    """Each feature's slot holds its widest label, tick or p and r line, with the gap; so the
    four kidney genes need more than a phone's plot, and the page scrolls the box. Mutations
    caught: the least width left off the layout, and the p labels left out of the measure."""
    data = box_data(kidney_like(), columns={"value": ("log2 CPM", "x")}, feature_label="Gene",
                    value_label="Expression", groups=("a", "b"))  # fmt: skip
    meta = data.interactive().to_plotly_json()["layout"]["meta"]
    ((scaled, fixed),) = least_width(data, (0.0, 10.0))
    tick, note = on_page(SIZES["tick"]), on_page(SIZES["annotation"])
    longest_note = max(text_width(line, note) for g in GENES
                       for line in data.caption.label(g).split("\n"))  # fmt: skip
    widest = max(text_width("NDUFA4L2", tick), longest_note)
    assert scaled >= len(GENES) * (widest + BOX_GAP_PX)
    assert scaled + fixed > PHONE_PLOT_PX
    short = box_data(kidney_like().with_columns(pl.col("feature").str.slice(0, 1)),
                     columns={"value": ("log2 CPM", "x")}, feature_label="Gene",
                     value_label="Expression", groups=("a", "b"))  # fmt: skip
    ((one_letter, _),) = least_width(short, (0.0, 10.0))
    assert one_letter >= len(set(short.frame["feature"])) * (longest_note + BOX_GAP_PX)
    assert meta["narrow"]["below"] > PHONE_PLOT_PX


def test_a_phone_turns_the_labels_and_the_box_fits_its_plot():
    """C3.4: the kidney page's stage box scrolled inside its card at 375 px. Below the upright
    layout's least width a phone's turns each feature's name, p and r into one quarter-turned
    tick label and hides the labels over the boxes, so four genes fit a phone's plot.

    Mutations caught: no narrow layout; the turned labels without the p and r lines; the
    least width left at the upright layout's, which scrolls the phone again."""
    data = box_data(kidney_like(), columns={"value": ("log2 CPM", "x")}, feature_label="Gene",
                    value_label="Expression", groups=("a", "b"))  # fmt: skip
    meta = data.interactive().to_plotly_json()["layout"]["meta"]
    ((scaled, fixed),) = meta["least_width"]
    assert scaled + fixed <= PHONE_PLOT_PX
    relayout = meta["narrow"]["relayout"]
    assert relayout["xaxis.tickangle"] == -90  # noqa: PLR2004 - a quarter turn
    for gene, tick in zip(GENES, relayout["xaxis.ticktext"], strict=True):
        assert tick.split("<br>") == [gene, data.caption.label(gene, "label_line")]
    hidden = [k for k, v in relayout.items() if k.startswith("annotations") and v is False]
    assert len(hidden) == len(GENES)
