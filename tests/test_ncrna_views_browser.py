"""M8's controls in Chromium (WO10): F3's volcano by biotype at the reader's cut-offs, and F7's
targets by context++ cutoff and class, on the made-up GDC root, at 1,280 px light and 375 px
dark. Each figure is read back from plotly and checked against Python's own numbers. Each test
names the mutations it catches.
"""

import re
import subprocess
import sys
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from test_export_browser import CHROMIUM, _browse
from views_browser import walk

from ctx.de.stats import classify
from ctx.ncrna import mirna, targets

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
CTX = Path(sys.executable).parent / "ctx"
SIZES = (("wide_light", 1280, "light"), ("phone_dark", 375, "dark"))
SETTINGS = [("padj_within", "0.05", "1"), ("pvalue", "0.01", "0.5"), ("padj_within", "0.1", "0"),
            ("pvalue", "0.001", "2")]  # fmt: skip
STRATA = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.strata'), fig = form.closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const read = () => {
    const out = {summary: fig.querySelector('.strata-summary').textContent,
                 moved: !!fig.querySelector('a[data-export][role=button]'), traces: {}};
    Object.keys(cfg.traces).forEach(l => {
      const t = d.data[cfg.traces[l]]; out.traces[l] = (t.customdata || []).map(c => c[0]);
    });
    out.lines = cfg.lines.map(i => d.data[i].x[0]);
    return out;
  };
  const out = {loaded: read(), settings: []};
  for (const [basis, cutoff, fc] of SETTINGS) {
    form.elements.basis.value = basis; form.elements.cutoff.value = cutoff;
    form.elements.fc.value = fc;
    form.elements.fc.dispatchEvent(new Event('change', {bubbles: true}));
    await wait(150); out.settings.push(read());
  }
  form.querySelector('button[type=reset]').click();
  await wait(200); out.reset = read();
  out.box = document.documentElement.clientWidth; out.scroll = document.documentElement.scrollWidth;
  return out;
})()"""


@pytest.fixture(scope="module")
def m8(tmp_path_factory):
    """The made-up GDC root's M8 folder, its target files, and what each act read."""
    import ncrna_tcga_root  # noqa: PLC0415 - imports synth, which only tests use

    path, spec = ncrna_tcga_root.make(tmp_path_factory.mktemp("tcga"), date.today())
    done = subprocess.run([str(CTX), "--root", str(path), "ncrna", "--dataset",
                           ncrna_tcga_root.DATASET, "--dataset-file", str(spec)],
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = path / "outputs" / "08_ncrna"
    (page,) = out.glob("*_overview_v1.html")
    act = STRATA.replace("SETTINGS", repr([list(s) for s in SETTINGS]))
    tasks = [{"name": name, "url": page.as_uri(), "width": width, "height": 900, "scheme": scheme,
              "act": act} for name, width, scheme in SIZES]  # fmt: skip
    seen = _browse(tasks, tmp_path_factory.mktemp("downloads"))
    views = walk(page, tmp_path_factory.mktemp("walk"), SIZES)
    folder = path / ncrna_tcga_root.TARGETS
    found = targets.read({"targetscan_predictions": folder / "predictions.zip",
                          "targetscan_families": folder / "families.zip",
                          "targetscan_genes": folder / "genes.zip",
                          "mirtarbase_strong": folder / "mirtarbase.csv"})  # fmt: skip
    return path, out, seen, views, found


def _table(out: Path, what: str) -> pl.DataFrame:
    return pl.read_parquet(next(out.glob(f"*_{what}_v1.parquet")))


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_volcano_calls_each_biotype_at_the_readers_cutoffs(m8, size):
    """Mutations caught: the pooled padj read for the stratum's; a called gene left grey; the
    fold-change lines left at the build's; reset leaving the view moved."""
    _, out, seen, _, _ = m8
    result = seen[size]
    assert result["errors"] == []
    acted = result["acted"]
    assert acted["scroll"] <= acted["box"]
    frame = _table(out, "results_by_biotype")
    levels = [b for b in acted["loaded"]["traces"] if b != "ns"]
    for (basis, cutoff, fc), drawn in zip(SETTINGS, acted["settings"], strict=True):
        calls = frame.select("gene_name", "biotype", call=classify(
            pl.col(basis), pl.col("log2FoldChange"), float(cutoff), float(fc)))  # fmt: skip
        for b in levels:
            want = calls.filter((pl.col("biotype") == b) & (pl.col("call") != "ns"))
            assert sorted(drawn["traces"][b]) == sorted(want["gene_name"].to_list()), (b, basis)
        assert drawn["lines"] == [-float(fc), float(fc)]
        assert drawn["moved"] == ((basis, cutoff, fc) != ("padj_within", "0.05", "1"))
    assert not acted["reset"]["moved"] and acted["reset"]["lines"] == [-1.0, 1.0]


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_targets_view_draws_each_cutoff_and_class(m8, size):
    """Mutations caught: every cutoff drawn at the default's; the classes drawn whatever is
    chosen; the build not first."""
    path, out, _, views, found = m8
    acted = views[size]["acted"]
    assert views[size]["errors"] == []
    (view,) = [v for v in acted["views"] if v["label"] == "Targets"]
    assert len(view["options"]) == 9  # noqa: PLR2004 - 3 cutoffs by 3 class choices
    first = view["options"][0]["title"].split(" of ")[1].split(" targets")[0]
    de = path / "outputs" / "05_de"
    results = pl.read_parquet(next(de.glob("*_results_v1.parquet")))
    tested = results.filter(pl.col("padj").is_not_null()).select(
        "gene_id", "gene_name", "log2FoldChange"
    )
    for option in view["options"]:
        line = option["line"] or "both classes apart, context++ at most -0.2;"
        preset = float(line.split("context++ at most ")[1].split(";")[0])
        curves = mirna.curves(found, first, tested, preset)
        names, present = {t["name"] for t in option["traces"]}, set(curves["target"])
        if "validated targets," in line:
            assert names == present - {targets.PREDICTED}
        elif "predicted-only targets," in line:
            assert names == present - {targets.VALIDATED}
        else:
            assert names == present
        for t in option["traces"]:
            full = curves.filter(pl.col("target") == t["name"])
            pairs = set(zip(full["log2FoldChange"].round(9), full["share"].round(9), strict=True))
            drawn = set(zip(np.round(t["x"], 9), np.round(t["y"], 9), strict=True))
            assert drawn <= pairs, (preset, t["name"])


def test_moved_cutoffs_export_the_page_and_the_genes_no_cutoff_calls(m8, tmp_path):
    """Mutations caught: the export switch left off; the files missing the genes no cut-off
    calls, or drawing the build's thinned ones over the moved calls."""
    path, out, _, _, _ = m8
    (page,) = out.glob("*_overview_v1.html")
    act = """(async () => {
      await new Promise(r => setTimeout(r, 1500));
      const form = document.querySelector('form.strata');
      form.elements.basis.value = 'pvalue'; form.elements.fc.value = '0.5';
      form.elements.fc.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 200));
      return form.closest('.fig').id;
    })()"""
    probe = _browse([{"name": "p", "url": page.as_uri(), "act": act}], tmp_path)["p"]
    (result,) = _browse([{"name": "e", "url": page.as_uri(), "act": act,
                          "figure": f"#{probe['acted']}"}], tmp_path).values()  # fmt: skip
    svg = (tmp_path / result["files"]["svg"]["guid"]).read_text()
    frame = _table(out, "results_by_biotype")
    # Every tested gene is a point in the file: a path element per marker plotly draws.
    assert len(re.findall(r'<path class="point', svg)) == frame.height
    assert "raw p below 0.05" in svg and "Moved, the page draws" in svg
