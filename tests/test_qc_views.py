"""M2's lookup views (work-order item 5, D-073): the sex figure against each sample's own
recorded sex, and the read flow by level. Each view is drawn in Python by the build's own
definition; the page swaps it in, and its PDF and SVG carry its own caption.

Each test names the mutation it catches.
"""

import json
import re
from datetime import date
from html import unescape

import polars as pl
import pytest
from synth import Synthetic, design_counts, write_gather_outputs
from test_export_browser import CHROMIUM, _browse
from test_qc_run import RNASEQ, RNASEQ_SAMPLES, SPEC, ctx

# The sample whose donor sex is planted as the other sex: its cells are its donor's on paper,
# its counts its own, so the method's view flags it and the own-sex view does not.
TRANSPLANT = RNASEQ_SAMPLES[1]


@pytest.fixture(scope="module")
def page(tmp_path_factory):
    """ctx qc on the rnaseq fixture's reports, one sample given a donor of the other sex."""
    root = tmp_path_factory.mktemp("qc_views")
    synthetic = design_counts((("control", "short", 3), ("case", "short", 2)))
    names = dict(zip(synthetic.samples["sample"], RNASEQ_SAMPLES, strict=True))
    renamed = Synthetic(synthetic.counts.rename(names),
                        synthetic.samples.with_columns(pl.col("sample").replace_strict(names)),
                        synthetic.truth)  # fmt: skip
    write_gather_outputs(root, "desynth", renamed, date(2026, 9, 30))
    (sheet_path,) = (root / "outputs" / "01_gather").glob("*_samplesheet_v1.parquet")
    sheet = pl.read_parquet(sheet_path)
    own = sheet.filter(pl.col("sample") == TRANSPLANT)["sex"][0]
    sheet.with_columns(donor_sex=pl.when(pl.col("sample") == TRANSPLANT)
                       .then(pl.lit("M" if own == "F" else "F")).otherwise(pl.lit(""))
                       ).write_parquet(sheet_path)  # fmt: skip
    done = ctx("--root", str(root), "qc", "--dataset", "desynth", "--dataset-file", str(SPEC),
               "--rnaseq-dir", str(RNASEQ), cwd=root)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "02_qc"
    (html,) = out.glob("*_overview_v1.html")
    return root, out, pl.read_parquet(sheet_path), html


def _views(html: str, title: str) -> dict:
    """The controls of the figure titled ``title``: its view options and their settings."""
    figure = next(
        f for f in html.split('<figure class="fig"') if f'class="fig-title">{title}<' in f
    )
    found = re.search(r'id="ctx-controls-\d+">(.*?)</script>', figure, re.S)
    assert found, title
    return json.loads(found.group(1))


def test_the_sex_view_compares_each_sample_with_its_own_record(page):
    """Mutation caught: the own-sex view built with the donor's sex, as the method is."""
    _, out, sheet, html = page
    views = pl.read_csv(next(out.glob("*_sex_views_v1.csv")))
    by = {v: dict(f.select("sample", "agrees").rows()) for (v,), f in views.group_by("view")}
    donor, own = (
        by["the donor's sex for a transplanted sample"],
        by["each sample's own recorded sex"],
    )
    assert donor[TRANSPLANT] is False and own[TRANSPLANT] is True
    recorded = dict(sheet.select("sample", "sex").rows())
    called = dict(views.select("sample", "expression_sex").rows())
    assert own == {s: recorded[s] == called[s] for s in own}
    cfg = _views(html.read_text(), "Recorded sex against expression")
    assert [o["label"] for o in cfg["options"]] == [
        "the donor's sex for a transplanted sample",
        "each sample's own recorded sex",
    ]
    first, second = cfg["options"]
    assert "compared with its donor" in first["export"]["caption"]
    assert "its own recorded sex" in second["export"]["caption"]
    assert "1 of 5 samples disagree" in unescape(first["lead"])
    assert "0 of 5 samples disagree" in unescape(second["lead"])


def test_each_flow_view_sums_its_own_samples(page):
    """Mutation caught: a level's flow summing every sample."""
    _, out, sheet, html = page
    fates = pl.read_csv(next(out.glob("*_desynth_fates_v1.csv")))
    links = pl.read_csv(next(out.glob("*_flow_views_v1.csv")))
    cfg = _views(html.read_text(), "Where every read goes")
    labels = [o["label"] for o in cfg["options"]]
    assert labels[0] == "all samples" and labels[1:] == list(dict.fromkeys(links["view"]))[1:]
    for condition in sheet["condition"].unique().sort():
        label = f"the samples whose condition is {condition}"
        assert label in labels
        kept = set(sheet.filter(pl.col("condition") == condition)["sample"])
        sequenced = fates.filter(pl.col("sample").is_in(kept))["total_reads"].sum()
        mine = links.filter((pl.col("view") == label) & (pl.col("source") == "total_reads"))
        assert mine["value"].sum() == pytest.approx(sequenced)
        option = cfg["options"][labels.index(label)]
        assert f"of the {len(kept)} samples" in option["export"]["caption"]


SWAP = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(2000);
  const fig = Array.from(document.querySelectorAll('figure.fig')).find(f =>
    (f.querySelector('.fig-title') || {}).textContent === %s);
  fig.setAttribute('data-probe', '1');
  const form = fig.querySelector('form.view');
  form.elements.view.value = String(form.elements.view.options.length - 1);
  form.dispatchEvent(new Event('change'));
  await wait(400);
  return {line: fig.querySelector('.view-summary').textContent,
          caption: fig.querySelector('figcaption').textContent};
})()"""


@pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
@pytest.mark.parametrize("title", ["Recorded sex against expression", "Where every read goes"])
def test_a_view_s_files_carry_its_own_caption(page, title, tmp_path):
    """In Chromium at 375 px: the last view shown, its SVG and PDF drawn with that view's
    caption and the live line. Mutation caught: the page exporting the build's caption."""
    _, _, _, html = page
    task = {"name": "v", "url": html.as_uri(), "width": 375, "height": 800,
            "act": SWAP % json.dumps(title), "figure": "[data-probe]"}  # fmt: skip
    seen = _browse([task], tmp_path)["v"]
    assert seen["errors"] == [] and seen["acted"]["line"].startswith("This view shows ")
    svg = (tmp_path / seen["files"]["svg"]["guid"]).read_text()
    assert "<image" not in svg
    words = " ".join(re.findall(r"<text[^>]*>([^<]*)</text>", svg))
    shown = " ".join(seen["acted"]["caption"].split())
    assert shown in unescape(words), (shown, words)
    assert (tmp_path / seen["files"]["pdf"]["guid"]).read_bytes().startswith(b"%PDF")
