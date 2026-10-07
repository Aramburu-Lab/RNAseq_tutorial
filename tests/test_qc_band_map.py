"""The QC band map names what a column is and each row's unit and cut-offs (Rytis 2026-10-01).

Its x-axis read "Sample, as nf-core/rnaseq names it." and its rows gave a title alone, so a
reader could tell neither what a column was nor a metric's unit or band. Both renderers draw
the words built here from conf/qc_metrics.toml and content/qc.toml.
"""

from pathlib import Path

import polars as pl
import pytest

from ctx.qc import figures
from ctx.qc.registry import load
from ctx.viz.style import THEMES

ROOT = Path(__file__).resolve().parents[1]
REG = load(ROOT / "conf" / "qc_metrics.toml")
# Every kind of row: a single cut-off up, two cut-offs down, two up, and a FastQC verdict.
ROWS = {
    "tg_total_reads": "Processed reads (reads; good from 10,000)",
    "percent_duplicates": "Duplicated reads (%; good to 30, usable to 50)",
    "reference_start_codons":
        "Start codons that read ATG (fraction; good from 0.99, usable from 0.9)",
    "fastqc_adapter": "Adapter content (status; pass good, warn usable, fail bad)",
}  # fmt: skip
# Rows in registry order.
LABELS = [ROWS[m] for m in REG.metrics if m in ROWS]
RUNS = ("SRR3", "SRR1", "SRR2")
X = "Each column one run: 3, in name order; hover names it"


@pytest.fixture(scope="module")
def band_map():
    values = pl.DataFrame({
        "run_accession": list(RUNS),
        **{REG.metrics[m].column: v for m, v in (
            ("tg_total_reads", [5_000.0, 20_000.0, 10_000.0]),
            ("percent_duplicates", [10.0, 40.0, 60.0]),
            ("reference_start_codons", [0.995, 0.95, 0.5]),
            ("fastqc_adapter", ["pass", "warn", "fail"]))},
    })  # fmt: skip
    data = figures.band_map(values, REG, "run_accession")
    assert data is not None
    return data


def test_the_processed_reads_metric_is_named_so():
    """Its title only: what it means and how it is computed stay as they were. Mutation: the
    registry title put back to "Reads written by Trim Galore"."""
    metric = REG.metrics["tg_total_reads"]
    assert metric.title == "Processed reads"
    assert (metric.meaning, metric.computed) == (
        "How many reads Trim Galore processed.", "Trim Galore total_reads")  # fmt: skip


OLD_TITLE = "Reads written by Trim Galore"


def test_the_old_title_is_nowhere_a_reader_sees_it():
    """Not in the registry, a page file, the glossary, a template or a fixture; nor on a QC
    page the reads fixture builds. Mutation: the old title appended to a band-map reading
    line in qc.toml, which every other test passed (review, 2026-10-01)."""
    src = ROOT / "src" / "ctx" / "report"
    places = [ROOT / "conf" / "qc_metrics.toml", src / "glossary.yaml",
              *(src / "content").glob("*.toml"), *(src / "templates").glob("*.j2"),
              *(p for p in (ROOT / "tests" / "fixtures").rglob("*") if p.is_file()
                and p.suffix in {".toml", ".tsv", ".csv", ".json", ".md", ".txt"})]  # fmt: skip
    found = [p for p in places if OLD_TITLE in p.read_text(errors="replace")]
    assert found == []


def test_rows_carry_unit_and_cut_offs_and_columns_their_count_and_order(band_map):
    """Mutations: the unit dropped from the row, one cut-off printed where two differ, the
    columns left in arrival order."""
    assert band_map.rows_in_order() == LABELS
    assert band_map.columns_in_order() == sorted(RUNS)
    assert band_map.column_title == X
    bands = dict(band_map.frame.filter(pl.col("run_accession") == "SRR3")
                 .select("metric", "band").rows())  # fmt: skip
    assert bands == dict(zip(ROWS.values(), ("bad", "good", "good", "good"), strict=True))


def _two_lines(label: str, joint: str) -> str:
    """A row label as an axis prints it: the metric, then its unit and cut-offs under it."""
    title, rest = label.split(" (", 1)
    return f"{title}{joint}({rest}"


def test_both_renderers_print_the_axis_title_and_row_labels(band_map):
    """The PDF and SVG, and the page; each row two lines high, the column labels upright, and
    neither the title nor a label off the page. On the page, since the merge with the heatmap
    fix (2026-10-07), rows are bands of a numeric axis and hover is each cell's own text.
    Mutations: either renderer back on the hover label; a row one line high under a two-line
    label; the page's rows named by their two-line labels, which hover printed broken; the
    phone's title left at the desktop's lines."""
    from ctx.viz.save import layout_problems  # noqa: PLC0415

    fig = band_map.static(THEMES["white"])
    (ax,) = [a for a in fig.axes if a.get_xlabel()]
    assert ax.get_xlabel().replace("\n", " ") == X
    assert [t.get_text() for t in ax.get_yticklabels()] == [_two_lines(r, "\n") for r in LABELS]
    assert layout_problems(fig) == []
    page = band_map.interactive()
    # Right-aligned under the upright column labels, wrapped to the plot at the narrowest
    # width each layout is drawn at: a desktop's, then a phone's.
    (title,) = [a for a in page.layout.annotations if a.text.replace("<br>", " ") == X]
    assert (title.x, title.xanchor, page.layout.xaxis.tickangle) == (1, "right", -90)
    from ctx.viz.metrics import text_width  # noqa: PLC0415
    from ctx.viz.style import PHONE_PLOT_PX, PLOTLY_MARGIN_PX, SIZES, on_page  # noqa: PLC0415

    phone = page.layout.meta["narrow"]["relayout"]
    assert phone["annotations[0].text"].replace("<br>", " ") == X
    room = PHONE_PLOT_PX - PLOTLY_MARGIN_PX["l"] - PLOTLY_MARGIN_PX["r"]
    assert all(text_width(line, on_page(SIZES["axis_label"])) <= room
               for line in phone["annotations[0].text"].split("<br>"))  # fmt: skip
    assert phone["annotations[0].text"].count("<br>") >= title.text.count("<br>")
    assert page.layout.margin.b > -title.yshift and phone["margin.b"] >= page.layout.margin.b
    # The axis prints each row on two lines; hover names it by its one-line name.
    assert list(page.layout.yaxis.ticktext) == [_two_lines(r, "<br>") for r in LABELS]
    for row, cells in zip(LABELS, page.data[0].text, strict=True):
        assert all(t.replace("<br>", " ").startswith(f"Metric: {row} Run: ") for t in cells)


def test_hover_names_the_band_and_the_legend_keys_have_a_size(band_map):
    """Hover's last line is the cell's band, where it read "None" on every cell; the legend's
    swatches have a size, so opened alone they double: the static's patch at the page's size
    (``PLOTLY_LEGEND_SWATCH_PX``, the heatmap fix's, where this branch had the page's point
    size; merged 2026-10-07). Mutations: the value column back in hover; the swatches left at
    plotly's own 6 px."""
    from ctx.viz.style import PLOTLY_LEGEND_SWATCH_PX  # noqa: PLC0415

    page = band_map.interactive()
    heat = page.data[0]
    rows, columns = band_map.rows_in_order(), band_map.columns_in_order()
    cells = {(r, c): t.split("<br>")[-1] for r, row in zip(rows, heat.text, strict=True)
             for c, t in zip(columns, row, strict=True)}  # fmt: skip
    bands = {(r, c): b for c, r, b in band_map.frame.select("run_accession", "metric", "band")
             .rows()}  # fmt: skip
    assert cells == bands
    assert not [t for row in heat.text for t in row if "None" in t]
    swatches = [t for t in page.data if t.type == "scatter" and t.showlegend is not False]
    assert [t.name for t in swatches] == ["good", "usable", "bad", "no band"]
    assert {t.marker.size for t in swatches} == {PLOTLY_LEGEND_SWATCH_PX}


def test_a_narrow_plot_wraps_its_row_labels_to_half_a_phone(band_map):
    """Narrower than the labels and a line per column need, the page swaps in labels wrapped
    to half a phone's plot and the height they take; every line fits, no word is split,
    and the rows keep their names. Mutation: the narrow labels left at two lines, which
    started up to 132 px left of a 375 px plot (review, 2026-10-01)."""
    from ctx.viz.metrics import text_width  # noqa: PLC0415
    from ctx.viz.style import PHONE_PLOT_PX, SIZES, on_page  # noqa: PLC0415

    page = band_map.interactive()
    narrow = page.layout.meta["narrow"]
    labels = narrow["relayout"]["yaxis.ticktext"]
    tick = on_page(SIZES["tick"])
    assert all(text_width(line, tick) <= PHONE_PLOT_PX / 2
               for label in labels for line in label.split("<br>"))  # fmt: skip
    assert [" ".join(label.split("<br>")) for label in labels] == LABELS
    assert max(label.count("<br>") for label in labels) >= 2
    assert narrow["relayout"]["height"] > page.layout.height
    assert PHONE_PLOT_PX < narrow["below"] < 686


def test_every_registry_row_fits_its_two_lines():
    """All the registry's banded and FastQC rows, one run each, in the figure files: each row
    two lines high, no label on another. Mutation: a row one line high under two lines."""
    from ctx.viz.save import layout_problems  # noqa: PLC0415

    metrics = [m for m in REG.metrics.values() if m.banded or m.kind == "status"]
    values = pl.DataFrame({"run_accession": ["SRR1"],
                           **{m.column: ["pass" if m.kind == "status" else m.good]
                              for m in metrics}})  # fmt: skip
    data = figures.band_map(values, REG, "run_accession")
    assert data is not None and len(data.rows_in_order()) == len(metrics) >= 10
    assert layout_problems(data.static(THEMES["white"])) == []


def test_a_row_is_as_tall_as_its_label_lines():
    """Thirty two-line rows stay apart in the figure files: each row takes two lines' height.
    Mutation: a row one line high under a two-line label."""
    from ctx.viz.heatmap import HeatmapData, Swatch  # noqa: PLC0415
    from ctx.viz.save import layout_problems  # noqa: PLC0415

    rows = [f"Metric {i} (reads; good from {i},000)" for i in range(30)]
    frame = pl.DataFrame({"row": rows, "col": ["S1"] * 30, "band": ["good"] * 30})
    frame = frame.with_columns(value=pl.lit(None, pl.Float64))
    data = HeatmapData(
        frame=frame, row="row", column="col", value="value", level="band",
        columns={c: ("none", c) for c in frame.columns}, row_label="Metric", column_label="Run",
        swatches=(Swatch("good", "#009E73"),), caption="Two-line rows; no test.",
        row_display={r: r.replace(" (", "\n(") for r in rows},
    )  # fmt: skip
    assert layout_problems(data.static(THEMES["white"])) == []
