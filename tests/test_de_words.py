"""M5's and M6's words after P2S.4: contrasts by plain labels, the VST's trend named, the GDC
acknowledged on M6 as on M5. Each test names the mutation it was seen to fail on (2026-10-14).
"""

import dataclasses
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
from test_function_words import SPEC, _enrichment, _gsea, _sets

from ctx.dataset import PcaSpec, load
from ctx.de.deseq import Fit
from ctx.de.page import every, level_words, subset_words, words
from ctx.de.run import Writer, _contrasts_table, _pca_figure, p_histogram
from ctx.function import run as m6
from ctx.gather import gdc
from ctx.report import text

ROOT = Path(__file__).resolve().parents[1]
STEP1 = load(ROOT / "conf" / "datasets" / "Non-cancerous-blood.toml")


def test_every_step1_contrast_reads_in_plain_words():
    """Rytis 2026-10-01: titles in plain words, never a contrast's identifier. Every Step-1
    contrast has its own label, no label holds an underscore, and the levels read as the SDRF
    spells them. Mutation: ``level_words`` returning the curated level (``steady_state``)."""
    said = every(STEP1.contrasts, STEP1.columns)
    assert all(c.label for c in STEP1.contrasts)
    assert len({w.label for w in said.values()}) == len(said)
    assert not [w.label for w in said.values() if "_" in w.label]
    assert level_words(STEP1.columns, "condition", "steady_state") == "steady state"
    assert level_words(STEP1.columns, "condition", "gcsf") == "G-CSF"
    assert said["ndn_gcsf"].tested == "G-CSF" and said["ndn_gcsf"].reference == "steady state"


def test_a_pca_names_its_samples_in_plain_words():
    """The Step-1 page titled its steady-state PCA "condition steady_state", a level's
    identifier, found reading the built page in Chromium (2026-10-14). Mutation: the filter's
    levels joined as curated."""
    (ss_ndn,) = [p for p in STEP1.de.pca if p.name == "ss_ndn"]
    shown = subset_words(ss_ndn.filters, STEP1.columns)
    assert "steady state" in shown and "_" not in shown


def test_contrasts_with_no_label_that_read_alike_are_told_apart():
    """The synthetic file labels none of its five case-control contrasts: each reads as its
    levels, and the identifier tells them apart, so no two bars share a name. Mutation:
    ``every`` returning ``words`` unchanged: five bars named "case against control"."""
    said = every(SPEC.contrasts, SPEC.columns)
    assert len({w.label for w in said.values()}) == len(SPEC.contrasts)
    alike = text.fill("de", "contrast.alike", label=words(SPEC.contrasts[0]).label,
                      name=SPEC.contrasts[0].name)  # fmt: skip
    assert said[SPEC.contrasts[0].name].label == alike


def _pca(counts: np.ndarray):
    """One PCA of made-up counts, two groups of six samples."""
    n = counts.shape[1]
    frame = pl.DataFrame({"gene_id": [f"g{i}" for i in range(counts.shape[0])],
                          **{f"s{j}": counts[:, j] for j in range(n)}})  # fmt: skip
    sheet = pl.DataFrame({"sample": [f"s{j}" for j in range(n)],
                          "donor": [f"d{j}" for j in range(n)],
                          "condition": ["a"] * (n // 2) + ["b"] * (n - n // 2)})  # fmt: skip
    spec = PcaSpec("t", {}, ("condition",), None, None, None)
    return _pca_figure(spec, sheet, frame, threads=1)[0]


def test_the_pca_caption_names_the_dispersion_trend_the_vst_used():
    """P2S.7: Step-1's ss_ndn VST fell back to a mean-based trend and no page said so. Poisson
    counts have a flat dispersion, where pydeseq2's parametric fit does not converge; the
    caption names the mean-based trend, and a negative binomial one names the parametric.
    Mutation: the caption filled with ``parametric`` whatever ``Pca.trend`` says."""
    rng = np.random.default_rng(0)
    means = np.repeat(np.exp(rng.uniform(np.log(5), np.log(5000), 800))[:, None], 12, axis=1)
    flat = _pca(rng.poisson(means))
    assert text.word("de", "figures.pca.trend.mean") in flat.caption
    spread = _pca(rng.negative_binomial(10, 10 / (10 + means)))
    assert text.word("de", "figures.pca.trend.parametric") in spread.caption
    assert text.word("de", "figures.pca.trend.mean") not in spread.caption


def test_m6_of_gdc_data_acknowledges_the_project_and_the_gdc(tmp_path):
    """Every page that draws GDC data acknowledges the TCGA project and the GDC (TODO.md,
    GDC's terms): M6 prints M5's sentences, verbatim from ui.toml, first among its results.
    Mutation: M6's results without ``_acknowledged``."""
    brca = load(ROOT / "conf" / "datasets" / "tcga_brca.toml")
    spec = dataclasses.replace(SPEC, gdc=brca.gdc)
    writer = Writer(tmp_path / "06_function", "desynth", date(2026, 10, 14), {})
    found, _, _ = m6._page(writer, spec, _enrichment(_sets(), _gsea(16, 4)))
    assert brca.gdc is not None
    said = gdc.acknowledgement(brca.gdc)
    assert tuple(found[: len(said)]) == said
    plain, _, _ = m6._page(writer, SPEC, _enrichment(_sets(), _gsea(16, 4)))
    assert not set(said) & {r for r in plain if isinstance(r, str)}


def test_cassetta_contrasts_read_by_their_labels():
    """D1.2: BRCA-immune-cells-macrophages, BRCA-immune-cells-monocytes and Breast-tissue-BRCA
    each had two contrasts of one pair of levels, which D-067 told apart by identifier.
    Labelled, every contrast reads by its own label on the contrasts table and in a caption,
    and no two read alike. Mutation: ``words``
    ignoring ``Contrast.label``: the table and caption print the levels and the identifier."""
    results = pl.DataFrame({"pvalue": [0.01, 0.5], "padj": [0.02, 0.6], "call": ["up", "ns"]})
    empty = pl.DataFrame()
    for slug in (
        "BRCA-immune-cells-macrophages",
        "BRCA-immune-cells-monocytes",
        "Breast-tissue-BRCA",
    ):
        spec = load(ROOT / "conf" / "datasets" / f"{slug}.toml")
        said = every(spec.contrasts, spec.columns)
        assert all(c.label for c in spec.contrasts), slug
        assert {n: w.label for n, w in said.items()} == {c.name: c.label for c in spec.contrasts}
        assert len({w.label for w in said.values()}) == len(said), slug
        fits = {c.name: Fit(c, results, empty, (), empty, empty) for c in spec.contrasts}
        shown = [row[0] for row in _contrasts_table(fits, said).rows]
        assert shown == [c.label for c in spec.contrasts], slug
        for c in spec.contrasts:
            caption = p_histogram(fits[c.name], said[c.name])[0].caption
            assert c.label in caption and c.name not in caption, (slug, c.name)
