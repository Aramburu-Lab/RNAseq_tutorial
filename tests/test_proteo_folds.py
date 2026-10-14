"""F4, F5 and the sample pairing of ctx proteo on planted data; each test names the mutation
it catches. Nothing here comes from a real dataset."""

import numpy as np
import polars as pl
import pytest
from scipy import stats

from ctx.dataset import Contrast
from ctx.de.stats import adjust
from ctx.proteo import ProteoError, figures, folds
from ctx.proteo.annotate import pair_samples
from ctx.viz.style import DATA_COLOURS, SPENT_GREYS, contrast

SHEET = pl.DataFrame({"case": ["C1", "C2", "C3"], "sample_id": ["S1", "S2", "S3"],
                      "sample": ["A1", "A2", "A3"]})  # fmt: skip


def test_pairs_by_case_then_by_sample_id():
    """A column named by its GDC sample ID pairs with that tumour, after the case columns.

    Mutation caught: the sample-ID branch removed (the freeze's two columns named so are left
    out, 119 tumours where 121 pair); the prefix not removed before the case lookup."""
    pairs, by_id = pair_samples(("XC1", "S2", "XC9"), SHEET, "X")
    assert pairs == [(0, "A1"), (1, "A2")] and by_id == 1


def test_two_columns_on_one_tumour_refused():
    """A case column and a sample-ID column of the same tumour would count it twice.

    Mutation caught: the duplicate check removed."""
    with pytest.raises(ProteoError, match="one tumour"):
        pair_samples(("XC1", "S1"), SHEET, "X")
    with pytest.raises(ProteoError, match="no protein sample"):
        pair_samples(("XC9",), SHEET, "X")


def planted_protein(n=8):
    """Four genes over n samples, the first half the reference: up, down, untested, small."""
    rng = np.random.default_rng(3)
    values = rng.normal(0, 0.2, (4, n))
    tested = np.arange(n) >= n // 2
    values[0, tested] += 3.0
    values[1, tested] -= 3.0
    values[2, tested] += 3.0
    values[2, [0, 1, n // 2, n // 2 + 1]] = np.nan  # half missing in each group: below 70 %
    values[3, tested] += 0.5  # significant but under |log2 FC| 1
    return ("UP", "DOWN", "HALF", "SMALL"), values, tested


def test_protein_calls_filter_bh_and_rule():
    """Filter, then Welch, then BH over the filtered proteins, then D-012's rule.

    Mutations caught: the 70 % filter skipped (HALF tested and called); BH over every protein
    (padj differs); classify given lfc_min 0 (SMALL called up); Welch replaced by Student's t
    (p differs)."""
    genes, values, tested = planted_protein()
    frame = folds.protein_calls(genes, values, tested, 0.7, 1.0)
    assert frame["tested"].to_list() == [True, True, False, True]
    assert frame["protein_call"].to_list() == ["up", "down", "ns", "ns"]
    assert frame["protein_padj"][2] is None
    p = [stats.ttest_ind(values[i, tested], values[i, ~tested], equal_var=False).pvalue
         for i in (0, 1, 3)]  # fmt: skip
    assert frame["protein_p"].to_numpy()[[0, 1, 3]] == pytest.approx(p)
    kept = frame.filter("tested")
    assert kept["protein_padj"].to_numpy() == pytest.approx(adjust(np.array(p), "padj"))
    assert frame["protein_padj"][3] < 0.05 and frame["protein_lfc"][3] < 1.0


def test_rna_calls_one_row_per_symbol_by_the_rule():
    """M5's untested genes and shared symbols leave; calls use D-012's rule.

    Mutations caught: the null-padj filter removed (an untested gene enters as ns); the
    shared-symbol filter removed (a symbol twice); the call made at |log2 FC| > 0."""
    results = pl.DataFrame({"gene_name": ["A", "B", "C", "D", "D", None],
                            "log2FoldChange": [1.0, 0.9, -2.0, 3.0, 3.0, 5.0],
                            "padj": [0.01, 0.001, None, 0.01, 0.01, 0.01]})  # fmt: skip
    calls = folds.rna_calls(results, 1.0)
    assert calls["gene"].to_list() == ["A", "B"]
    assert calls["rna_call"].to_list() == ["up", "ns"]


def test_kind_and_opposite():
    """Where a gene was called, and opposite calls counted among those at both layers.

    Mutations caught: the RNA and protein branches of ``kind`` swapped; opposite counted over
    every gene instead of those called at both."""
    genes = pl.DataFrame({"rna_call": ["up", "up", "ns", "ns", "down"],
                          "protein_call": ["up", "ns", "down", "ns", "up"]})  # fmt: skip
    kinds = genes.select(kind=folds.kind(pl.col("rna_call"), pl.col("protein_call")))["kind"]
    assert kinds.to_list() == ["both", "rna", "protein", "neither", "both"]
    made = folds.Folds(Contrast("c", "v", "r", "t", (1, 1)),
                       genes.with_columns(kind=kinds), (1, 1), (None, None, None))  # fmt: skip
    counts = made.counts()
    assert (counts["both"], counts["opposite"], counts["rna_up"], counts["protein_down"]) == (
        2,
        1,
        2,
        1,
    )


def test_contrast_folds_uses_the_tumours_m5_compared():
    """Protein is tested on the paired tumours at the contrast's two levels, by M5's selection.

    Mutations caught: every paired tumour used (the LumB column's planted 50 enters the
    reference mean); groups taken from column order instead of the sheet's variable."""
    sheet = pl.DataFrame(
        {
            "sample": [f"A{i}" for i in range(7)],
            "donor": [f"C{i}" for i in range(7)],
            "run_accession": [f"R{i}" for i in range(7)],
            "pam50": ["Basal", "LumA", "Basal", "LumA", "Basal", "LumA", "LumB"],
        }
    )
    pairs = [(j, f"A{j}") for j in range(7)]
    values = np.tile([0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 50.0], (3, 1))
    values += np.random.default_rng(1).normal(0, 0.01, values.shape)
    results = pl.DataFrame({"gene_name": ["G0", "G1", "G2"], "log2FoldChange": [-2.0, 0.0, 1.0],
                            "padj": [0.001, 0.9, 0.9]})  # fmt: skip
    contrast = Contrast("basal_vs_luma", "pam50", "LumA", "Basal", (3, 3))
    protein = (("G0", "G1", "G2"), values)
    made = folds.contrast_folds(results, sheet, pairs, protein, contrast, order={},
                                limits=(0.7, 1.0))  # fmt: skip
    assert made.n_protein == (3, 3)
    assert made.genes["protein_lfc"].to_numpy() == pytest.approx([-2.0] * 3, abs=0.05)
    assert made.genes["kind"].to_list() == ["both", "protein", "protein"]
    with pytest.raises(ProteoError, match="no column 'er_status'"):
        folds.contrast_folds(results, sheet, pairs, protein,
                             Contrast("er", "er_status", "ER+", "ER-", (3, 3)), order={},
                             limits=(0.7, 1.0))  # fmt: skip


def test_zscores_keep_missing_missing():
    """A row is scaled over its present values; a missing cell stays missing, a flat row too.

    Mutation caught: NaN filled with 0 before scaling (the row's mean moves). A constant row
    is 0 / 0 and stays NaN, never a colour."""
    z = figures.zscores(np.array([[1.0, np.nan, 3.0], [2.0, 2.0, 2.0]]))
    assert z[0, 0] == pytest.approx(-1.0) and z[0, 2] == pytest.approx(1.0)
    assert np.isnan(z[0, 1]) and np.isnan(z[1]).all()


def test_heatmap_genes_named_first_then_both_ends():
    """F5 draws the named tested genes, then the highest and the lowest rho, once each.

    Mutations caught: an untested named gene drawn; the low end taken from the top again;
    a gene drawn twice when the two ends overlap."""
    genes = pl.DataFrame({"gene": ["N", "U", "a", "b", "e", "c", "d"],
                          "tested": [True, False, True, True, True, True, True],
                          "rho": [0.1, 0.9, 0.8, 0.5, 0.3, -0.2, -0.7]})  # fmt: skip
    assert figures.heatmap_genes(genes, ["N", "U"], 2) == ["N", "a", "b", "c", "d"]
    assert figures.heatmap_genes(genes, [], 4) == ["a", "b", "e", "N", "c", "d"]


def test_plex_strip_alternates_two_greys_one_run_per_plex():
    """22 plexes, over the colour cap, take two greys in turn and sort first, so every plex
    boundary changes grey (D-063); five subtypes keep their data colours.

    Mutations caught: the plexes given data colours in turn (as before D-063); the subtype
    sorted first, which breaks plexes into runs whose neighbours can share a grey; the greys
    alternating by anything but display order (plex 2 before plex 10)."""
    plexes = [str(p) for p in range(1, 23)]
    subtypes = ["Basal", "Her2", "LumA", "LumB", "Normal-like"]
    annotation = pl.DataFrame({"PAM50": [subtypes[i % 5] for i in range(44)],
                               "Plex": [plexes[(i * 7) % 22] for i in range(44)]})  # fmt: skip
    tracks = (("PAM50", "PAM50"), ("Plex", "Plex"))
    assert figures.sort_order(annotation, tracks) == [tracks[1], tracks[0]]
    greys = figures.track_colours(figures.by_label(annotation["Plex"]))
    assert greys["2"] == greys["10"] != greys["1"] == greys["3"]
    assert set(greys.values()) == set(SPENT_GREYS)
    by = [pl.col("Plex").str.len_chars(), pl.col("Plex"), pl.col("PAM50")]
    shown = annotation.sort(by)["Plex"].to_list()
    assert all(greys[a] != greys[b] for a, b in zip(shown, shown[1:], strict=False) if a != b)
    assert contrast(*SPENT_GREYS) >= 3  # noqa: PLR2004 - the 3:1 graphics threshold
    assert set(figures.track_colours(subtypes).values()) <= set(DATA_COLOURS)
