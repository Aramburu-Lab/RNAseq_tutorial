"""M3's controls in Chromium (WO8) on the chr22_se fixture's page, at 1,280 px light and 375 px
dark: the counter scatter by gene type and sample, the agreement plot by limits and gene type,
and the forward share's strand calls at the reader's thresholds. Each figure is read back from
plotly and held to numbers computed here from the tables beside the page. Each test names the
mutations it catches.
"""

from dataclasses import replace
from fractions import Fraction

import numpy as np
import polars as pl
import pytest
from test_export_browser import CHROMIUM, _browse
from views_browser import walk

from ctx import concordance
from ctx.process import tables
from ctx.process import views as views_module
from ctx.qc import registry

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
SIZES = (("wide_light", 1280, "light"), ("phone_dark", 375, "dark"))
SETTINGS = [["0.6", "0.3"], ["0.95", "0.05"], ["0.8", "0.1"]]
STRAND = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.strand'), fig = form.closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const read = () => {
    const out = {summary: fig.querySelector('.strand-summary').textContent, calls: {},
                 moved: !!fig.querySelector('a[data-export][role=button]'),
                 lines: cfg.lines.map(i => d.data[i].y[0])};
    cfg.calls.forEach(c => {
      out.calls[c] = (d.data[cfg.traces[c]].customdata || []).map(r => r[0]);
    });
    return out;
  };
  const out = {loaded: read(), settings: []};
  for (const [st, un] of SETTINGS) {
    form.elements.stranded.value = st; form.elements.unstranded.value = un;
    form.elements.unstranded.dispatchEvent(new Event('change', {bubbles: true}));
    await wait(150); out.settings.push(read());
  }
  out.box = document.documentElement.clientWidth; out.scroll = document.documentElement.scrollWidth;
  return out;
})()"""


@pytest.fixture(scope="module")
def m3(process_root, tmp_path_factory):
    """The M3 folder, every view walked and the strand act at both sizes."""
    out = process_root[0] / "outputs" / "03_process"
    (page,) = out.glob("*_overview_v1.html")
    act = STRAND.replace("SETTINGS", repr(SETTINGS))
    tasks = [{"name": n, "url": page.as_uri(), "width": w, "height": 900, "scheme": s,
              "act": act} for n, w, s in SIZES]  # fmt: skip
    return out, walk(page, tmp_path_factory.mktemp("walk"), SIZES), _browse(
        tasks, tmp_path_factory.mktemp("strand"))  # fmt: skip


def _table(out, what: str) -> pl.DataFrame:
    (path,) = out.glob(f"*_{what}_v1.parquet")
    return pl.read_parquet(path)


def _typed(out) -> pl.DataFrame:
    long = _table(out, "gene_comparison")
    return long.with_columns(type_group=pl.when(pl.col("gene_type").is_in(tables.GENE_TYPES))
                             .then("gene_type").otherwise(pl.lit(tables.OTHER)))  # fmt: skip


TYPE_WORDS = {"Every gene": None, "Protein-coding": "protein_coding", "lncRNA": "lncRNA",
              "Other": "other"}  # fmt: skip


def _genes(part: pl.DataFrame, line: str) -> tuple[pl.DataFrame, str | None, str]:
    """A view's genes from its live line: the type and the sample it names."""
    named = (line or "Every gene; every sample summed;").split("shows ")[-1].split(";")
    gene_type, unit = TYPE_WORDS[named[0].strip()], named[1].strip()
    genes = part if gene_type is None else part.filter(pl.col("type_group") == gene_type)
    return genes, gene_type, unit


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_scatter_draws_each_type_and_sample(m3, size):
    """Mutations caught: every view drawn over every gene; a sample's view drawn from the sums;
    the sample axis offered past its cap."""
    out, seen, _ = m3
    assert seen[size]["errors"] == []
    acted = seen[size]["acted"]
    assert acted["scroll"] <= acted["box"]
    long = _typed(out)
    samples = long["sample"].unique().to_list()
    views = [v for v in acted["views"] if v["options"][0]["title"].startswith("Gene counts")]
    assert views
    for view in views:
        n_units = 1 + (len(samples) if len(samples) <= 8 else 0)  # noqa: PLR2004
        name = view["options"][0]["title"].split(": ")[1]
        part = long.filter(pl.col("comparison") == _comparison(long, name))
        units = [concordance.ALL, *(samples if len(samples) <= 8 else [])]  # noqa: PLR2004
        assert len(view["options"]) == len(views_module.types(part, units)) * n_units
        for option in view["options"]:
            genes, _, unit = _genes(part, option["line"])
            counts = (genes.group_by("gene_id").agg(pl.col("a").sum(), pl.col("b").sum())
                      if unit.startswith("every") else
                      genes.filter(pl.col("sample") == unit.removeprefix("sample ")))  # fmt: skip
            drawn = sum(len(t["x"] or []) for t in option["traces"] if t["name"] != "Equal counts")
            assert drawn == counts.filter(pl.col("a") + pl.col("b") > 0).height, option["line"]


def _comparison(long: pl.DataFrame, name: str) -> str:
    """A comparison's identifier from its words in a title."""
    from ctx.process.figures import comparison_words  # noqa: PLC0415

    return next(c for c in long["comparison"].unique() if comparison_words(c)["name"] == name)


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_agreement_plot_counts_beyond_each_limit(m3, size):
    """Mutations caught: the empirical view drawn at 1.96 SD; a type's limits taken from every
    gene."""
    out, seen, _ = m3
    long = _typed(out)
    views = [v for v in seen[size]["acted"]["views"]
             if v["options"][0]["title"].startswith("How far")]  # fmt: skip
    assert views
    for view in views:
        name = view["options"][0]["title"].split(": ")[1]
        part = long.filter(pl.col("comparison") == _comparison(long, name))
        for option in view["options"]:
            line = option["line"] or "Every gene; limits bias ± 1.96 SD;"
            genes, _, _ = _genes(part, line.replace("limits ", ""))
            sums = genes.group_by("gene_id").agg(pl.col("a").sum(), pl.col("b").sum())
            found = concordance.agreement(sums["a"], sums["b"])
            kept = sums.filter(pl.col("a") + pl.col("b") > 0)
            d = np.log2(kept["b"].to_numpy() + 1) - np.log2(kept["a"].to_numpy() + 1)
            low, high = ((found.d_p2_5, found.d_p97_5) if "percentiles" in line
                         else (found.loa_low, found.loa_high))  # fmt: skip
            if low is None:
                continue
            beyond = int(((d < low) | (d > high)).sum())
            assert f"{beyond:,} of {kept.height:,} genes lie beyond" in option["lead"], line


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_the_strand_calls_follow_the_thresholds_as_the_registry_decides(m3, size):
    """Mutations caught: a sample left in its build trace; the lines left at the registry's;
    the live line not refilled."""
    out, _, strand = m3
    assert strand[size]["errors"] == []
    acted = strand[size]["acted"]
    base = registry.load(registry.resolve(None)).metrics["gene_strand_bias"]
    totals = _table(out, "strandedness").select("sample", "assigned_s1", "assigned_s2")
    for (st, un), drawn in zip(SETTINGS, acted["settings"], strict=True):
        metric = replace(base, stranded=float(st), unstranded=float(un))
        want: dict[str, list[str]] = {}
        for sample, s1, s2 in totals.iter_rows():
            if s1 + s2:
                call = metric.strandedness(Fraction(s1) / (Fraction(s1) + Fraction(s2)))
                want.setdefault(call, []).append(sample)
        assert {c: sorted(v) for c, v in drawn["calls"].items() if v} == {
            c: sorted(v) for c, v in want.items()
        }
        s, half = float(st), float(un) / 2
        assert np.allclose(drawn["lines"], [s, 1 - s, 0.5 - half, 0.5 + half])
        assert drawn["summary"].startswith(f"At stranded {st} and unstranded {un}:")
        assert drawn["moved"] == ((st, un) != ("0.8", "0.1"))
