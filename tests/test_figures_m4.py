"""F.1: M4's figures, each redrawn from its spec's inputs alone in a clean folder.

Every figure with a spec writes M4's data table, PDF and SVG byte for byte; the figures whose
stem a later table takes, and the heatmap of the count matrix, have none. Each test names its
mutation.
"""

import polars as pl
from test_figures import drawn_alike, redrawn_each, specs

from ctx.figures.spec import parse

QUANT = ("04_quant",)


def test_fates_and_estimates_redraw_byte_for_byte(quant_built, tmp_path):
    """The fates bars and the estimates scatter, redrawn from M4's tables. Catches the sheet
    read without its fill words, the contrast's levels counted the other way round, and the
    fold changes written after their figure, which leaves the estimates without a spec."""
    root, spec, _ = quant_built
    found = specs(root, QUANT, spec)
    assert sorted(parse(s.figure_id)[2] for s in found) == ["estimates", "fates"]
    for one, csv in zip(found, redrawn_each(found, root, tmp_path), strict=True):
        module_csv = next(root.glob(f"outputs/04_quant/{one.stem}.csv"))
        assert drawn_alike(csv, module_csv) == [], one.stem


def test_the_contrast_samples_are_m5s(quant_built):
    """The contrast's samples M4 wrote are the featured contrast's, as M5 selects them. Catches
    the samples of the whole sheet written instead."""
    from ctx.dataset import load  # noqa: PLC0415
    from ctx.de.deseq import select_samples  # noqa: PLC0415

    root, spec_path, _ = quant_built
    spec = load(spec_path)
    out = root / "outputs"
    written = pl.read_parquet(next(out.glob("04_quant/*_contrast_samples_v1.parquet")))
    sheet = pl.read_parquet(next(out.glob("01_gather/*_samplesheet_v1.parquet")))
    used, _ = select_samples(sheet.unique("sample", keep="first", maintain_order=True),
                             spec.contrasts[0], spec.order)  # fmt: skip
    assert written.equals(used.select("sample", "group"))
