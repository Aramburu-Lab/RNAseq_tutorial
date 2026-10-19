"""M4.14: the quantification page's views in Chromium on the quant fixture, at 1,280 px light
and 375 px dark: F2 by detection threshold, F3 by sample order, F4 by tximport's counts mode,
F6 by gene type. Each view is read back from plotly and held to the tables beside the page and
to numbers computed here. Each test names the mutations it catches.
"""

import re

import numpy as np
import polars as pl
import pytest
from test_export_browser import CHROMIUM, _browse
from views_browser import walk

from ctx.quant import tables

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
SIZES = (("wide_light", 1280, "light"), ("phone_dark", 375, "dark"))


@pytest.fixture(scope="module")
def m4(quant_built, tmp_path_factory):
    """The M4 folder and every view walked at both sizes."""
    root, _, written = quant_built
    out = root / "outputs" / "04_quant"
    return out, walk(written[0], tmp_path_factory.mktemp("walk"), SIZES)


def _table(out, what: str) -> pl.DataFrame:
    (path,) = out.glob(f"*_{what}_v1.parquet")
    return pl.read_parquet(path)


def _view(seen: dict, size: str, label: str) -> dict:
    (view,) = [v for v in seen[size]["acted"]["views"] if v["label"] == label]
    return view


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_each_threshold_draws_its_own_count_of_the_same_draws(m4, size):
    """Mutations caught: every view drawn at 10 reads; a threshold drawn from other draws."""
    out, seen = m4
    assert seen[size]["errors"] == []
    page = seen[size]["acted"]
    assert page["scroll"] <= page["box"]
    view = _view(seen, size, "Detected at")
    table = _table(out, "detection")
    assert [o["line"].split("count of ")[1].split(" or more")[0] if o["line"] else "10"
            for o in view["options"]] == ["10", "1", "100"]  # fmt: skip
    for option, column in zip(view["options"], ("genes", "genes_1", "genes_100"), strict=True):
        drawn = sorted(y for t in option["traces"] if t["name"] != "Depth every sample reaches"
                       for y in t["y"])  # fmt: skip
        assert drawn == sorted(table[column].to_list()), column
    assert (table["genes_1"] >= table["genes"]).all()
    assert (table["genes"] >= table["genes_100"]).all() and table["genes_1"].max() > 0


def test_thinning_counts_every_threshold_on_the_same_draws():
    """Mutations caught: an extra threshold drawn with the generator again, which moves the
    build's own curve; a threshold compared with > rather than >=."""
    rng = np.random.default_rng(3)
    counts = rng.poisson(rng.gamma(0.4, 60, (500, 1)), (500, 4)).astype(float)
    plain, common = tables.thinned(counts, ["a", "b", "c", "d"])
    extra, again = tables.thinned(counts, ["a", "b", "c", "d"], (1, 100))
    assert common == again and plain.equals(extra.select(plain.columns))
    whole = extra.filter(pl.col("depth") == pl.col("depth").max().over("sample"))
    for j, sample in enumerate("abcd"):
        row = whole.filter(pl.col("sample") == sample).row(0, named=True)
        if row["depth"] == counts[:, j].round().sum():
            assert row["genes_1"] == int((np.rint(counts[:, j]) >= 1).sum())


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_counts_modes_redraw_the_counts_and_the_types_their_genes(m4, size):
    """Mutations caught: every mode drawn from the estimates; a type's view drawing every
    gene; the estimates view built when a fixture holds one type."""
    out, seen = m4
    units = _view(seen, size, "Gene counts as")
    ys = [sorted(y for t in o["traces"] for y in (t["y"] or [])) for o in units["options"]]
    assert len(ys) == 3 and ys[0] != ys[1] and ys[1] != ys[2]
    estimates = [v for v in seen[size]["acted"]["views"] if v["label"] == "Genes"]
    frame = _table(out, "fold_changes")
    for view in estimates:
        for option in view["options"][1:]:
            gene_type = option["line"].split("shows ")[1].split(";")[0]
            names = {
                t["name"] for t in option["traces"] if t["name"] != "Equal fold changes" and t["x"]
            }
            assert names <= {gene_type}
            n = frame.filter(pl.col("gene_type") == gene_type).height
            assert f"of {n:,} genes" in option["lead"]


ORDER = "Order samples by"
CELLS = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('FIGURE form.view');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  form.elements.view.value = 'OPTION';
  form.elements.view.dispatchEvent(new Event('change', {bubbles: true}));
  await wait(200);
  const o = cfg.options[OPTION];
  return {z: o.figure.data.filter(t => t.type === 'heatmap').map(t => t.z),
          cells: o.export.cells, caption: o.export.caption, label: o.label};
})()"""


def _runs(z: list[list[float | None]]) -> dict[int, int]:
    """Rectangles per z value: each run of equal z along a row is one, a missing z none."""
    found: dict[int, int] = {}
    for row in z:
        for i, v in enumerate(row):
            if v is not None and (i == 0 or row[i - 1] != v):
                found[int(v)] = found.get(int(v), 0) + 1
    return found


def _rgb(colour: str) -> str:
    """A colour as plotly writes a fill: rgb() of a hex, rgba() kept."""
    if colour.startswith("#"):
        r, g, b = (int(colour[i : i + 2], 16) for i in (1, 3, 5))
        return f"rgb({r}, {g}, {b})"
    return colour


def test_each_order_groups_its_samples_by_its_column(m4):
    """Mutations caught: a grouping view drawn in the clustering's order; a column with one
    level offered; the strips of the build kept where the grouping's column has none."""
    out, seen = m4
    (view,) = [v for v in seen["wide_light"]["acted"]["views"] if v["label"] == ORDER]
    assert view["options"][0]["line"] == "" and all(o["moved"] for o in view["options"][1:])
    assert len(view["options"]) >= 2  # noqa: PLR2004 - the build and one grouping
    for option in view["options"][1:]:
        assert "grouped by" in option["caption"]


def test_a_moved_grouping_view_exports_vector_cells_and_no_raster(m4, tmp_path):
    """HeatmapData's export layers (V2): plotly's SVG draws a heatmap as one raster image; the
    view's PDF and SVG draw each run of equal cells as a filled path in the files' colours.

    Mutations caught: the heatmap traces exported as plotly draws them, a raster; a run split
    or merged across a level change; a track's cells painted in the cells' colours; the export
    caption left on the build's while the view shows another order."""
    out, seen = m4
    (view,) = [v for v in seen["wide_light"]["acted"]["views"] if v["label"] == ORDER]
    (page,) = out.glob("*_quant_*_overview_v1.html")
    figure = f"#{view['figure']}"
    act = CELLS.replace("FIGURE", figure).replace("OPTION", "1")
    (result,) = _browse([{"name": "cells", "url": page.as_uri(), "width": 1280, "figure": figure,
                          "act": act}], tmp_path).values()  # fmt: skip
    acted = result["acted"]
    svg = (tmp_path / result["files"]["svg"]["guid"]).read_text()
    assert "<image" not in svg
    fills = re.findall(r'<path class="js-fill" d="([^"]+)" style="[^"]*?fill: ([^;]+);', svg)
    drawn: dict[str, int] = {}
    for path, colour in fills:
        drawn[colour] = drawn.get(colour, 0) + path.count("Z")
    expected: dict[str, int] = {}
    for z, colours in zip(acted["z"], acted["cells"], strict=True):
        for v, n in _runs(z).items():
            key = _rgb(colours[v])
            expected[key] = expected.get(key, 0) + n
    assert drawn == expected
    assert "grouped by" in acted["caption"]
    assert (tmp_path / result["files"]["pdf"]["guid"]).read_bytes().startswith(b"%PDF")
