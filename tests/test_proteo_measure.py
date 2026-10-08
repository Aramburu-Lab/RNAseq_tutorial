"""ctx.proteo.measure and tables on planted data; each test names the mutation it catches."""

import io
import tarfile

import numpy as np
import pytest
from scipy import stats

from ctx.de.stats import adjust
from ctx.proteo import ProteoError, measure, tables

RNG = np.random.default_rng(11)


def planted(n_genes=40, n_samples=10, missing=0.0, seed=5):
    rng = np.random.default_rng(seed)
    rna = rng.normal(5, 1, (n_genes, n_samples))
    protein = 0.6 * rna + rng.normal(0, 0.5, (n_genes, n_samples))
    protein[rng.random(protein.shape) < missing] = np.nan
    return tuple(f"G{i}" for i in range(n_genes)), rna, protein


def test_missing_is_never_zero():
    """A gene's correlation uses its complete pairs only.

    Mutation caught: missing protein filled with 0 before the correlation."""
    genes, rna, protein = planted(n_genes=1, n_samples=12)
    protein[0, [2, 7]] = np.nan
    frame = measure.per_gene(genes, rna, protein, minimum=0.5)
    keep = ~np.isnan(protein[0])
    assert frame["rho"][0] == pytest.approx(stats.spearmanr(rna[0, keep], protein[0, keep])[0])
    assert frame["r"][0] == pytest.approx(stats.pearsonr(rna[0, keep], protein[0, keep])[0])
    assert frame["rho_p"][0] == pytest.approx(stats.spearmanr(rna[0, keep], protein[0, keep])[1])
    assert frame["n_quantified"][0] == 10


@pytest.mark.parametrize(("present", "total", "minimum", "passes"),
                         [(7, 10, 0.7, True), (6, 10, 0.7, False), (69, 99, 0.7, False),
                          (70, 99, 0.7, True), (55, 100, 0.55, True)])  # fmt: skip
def test_quantified_edge(present, total, minimum, passes):
    """The edge both ways, in exact arithmetic, at the default and at a share a reader may set.

    Mutations caught: ``present >= minimum * total`` (0.55 * 100 is 55.00000000000001, so 55
    of 100 fails; at 0.7 floats happen to land exactly for every total up to 300); a strict
    ``>`` (7 of 10 fails)."""
    assert measure.enough(np.array([present]), total, minimum)[0] == passes


def test_bh_over_tested_genes_only():
    """BH counts the tested genes; a gene below the share is neither tested nor counted.

    Mutation caught: untested genes entered into BH with p = 1, which inflates every padj."""
    genes, rna, protein = planted(n_genes=30, n_samples=10)
    protein[:10, :5] = np.nan
    frame = measure.per_gene(genes, rna, protein)
    tested = frame.filter("tested")
    assert tested.height == 20 and frame.filter(~frame["tested"])["rho"].null_count() == 10
    expected = adjust(tested["rho_p"].to_numpy(), "padj")
    assert tested["rho_padj"].to_numpy() == pytest.approx(expected)


def test_order_one_row_per_gene_before_the_filter():
    """The order selftest: collapse to one row per gene, then filter, then correlate, then BH.

    A gene measured twice keeps its better-quantified row, so it passes the filter.
    Mutation caught: the first row kept (the gene would fall below 70 % and leave the test)."""
    values = np.array([[np.nan] * 5 + [1.0] * 5, [1.0, 2, 3, 4, 5, 6, 7, 8, np.nan, 9]])
    kept, report = tables.one_per_gene(["ERBB2", "ERBB2"], values)
    assert kept == [1] and report.extra_rows == 1 and report.genes_twice == 1
    rna = np.arange(10, dtype=float)[None, :]
    frame = measure.per_gene(("ERBB2",), rna, values[kept])
    assert frame["tested"][0] and frame["rho"][0] == pytest.approx(1.0)
    with pytest.raises(ValueError, match="twice"):
        measure.per_gene(("A", "A"), np.vstack([rna, rna]), values)


def test_one_per_gene_report():
    """Mutations caught: a group naming several genes kept under its first; an empty symbol kept."""
    values = np.ones((5, 3))
    kept, report = tables.one_per_gene(["A", "B;C", None, "", "D"], values)
    assert kept == [0, 4]
    assert (report.rows, report.no_symbol, report.several_genes, report.kept) == (5, 2, 1, 2)


def test_membership_drops_genes_in_two_classes():
    """Mutation caught: a gene in both classes counted in the first (ribosomal proteins sit in
    Reactome's amino acid metabolism, so complex against metabolic would compare them twice)."""
    labels, twice = measure.membership(
        ["RPL3", "IDH1", "X"], {"complex": {"RPL3"}, "metabolic": {"RPL3", "IDH1"}}
    )
    assert labels == ["", "metabolic", "other"] and twice == 1


def test_planted_class_shift_and_dunn():
    """A planted complex class with lower rho than a metabolic class is found, one-sided.

    Mutations caught: the two-sided p; the rank-biserial with U of the other group (sign
    flipped); Dunn's z without the tie correction (here ties are planted)."""
    low = np.round(RNG.normal(0.1, 0.1, 40), 2)
    high = np.round(RNG.normal(0.4, 0.1, 40), 2)
    p, r = measure.shifted(high, low)
    one_sided = stats.mannwhitneyu(high, low, alternative="greater").pvalue
    assert p == pytest.approx(one_sided, rel=1e-9, abs=0)
    assert p < 0.001 and r > 0.8
    test = measure.class_test({"complex": low, "metabolic": high})
    pooled = np.concatenate([low, high])
    ranks = stats.rankdata(pooled)
    n = pooled.size
    _, t = np.unique(pooled, return_counts=True)
    assert (t > 1).any()
    var = (n * (n + 1) / 12 - (t**3 - t).sum() / (12 * (n - 1))) * (2 / 40)
    z = (ranks[:40].mean() - ranks[40:].mean()) / np.sqrt(var)
    assert test.pairs[0][2] == pytest.approx(z)
    assert test.h == pytest.approx(stats.kruskal(low, high)[0])
    assert test.epsilon2 == pytest.approx(test.h / (n - 1))


def test_compressed_ratios_keep_ranks_and_shrink_changes():
    """Ratio compression scales protein changes down: ranks and r stay, fold changes shrink.

    Mutation caught: the protein fold change taken as a ratio of means of unlogged values."""
    genes, rna, protein = planted(n_genes=5, n_samples=12)
    full = measure.per_gene(genes, rna, protein)
    squeezed = measure.per_gene(genes, rna, 0.4 * protein)
    assert squeezed["rho"].to_list() == pytest.approx(full["rho"].to_list())
    assert squeezed["r"].to_list() == pytest.approx(full["r"].to_list())
    group = np.arange(12) < 6
    lfc, _ = measure.welch(protein, group)
    lfc_squeezed, _ = measure.welch(0.4 * protein, group)
    assert lfc_squeezed == pytest.approx(0.4 * lfc)


def test_welch_needs_two_values_per_group():
    """Mutation caught: a gene with one present value in a group given a p."""
    values = np.array([[1.0, np.nan, np.nan, 2.0, 3.0, 4.0], [1.0, 1.5, 2.0, 2.5, 4.0, 5.0]])
    lfc, p = measure.welch(values, np.array([True, True, True, False, False, False]))
    assert np.isnan(p[0]) and not np.isnan(p[1])
    assert lfc[0] == pytest.approx(1.0 - 3.0)
    expected = stats.ttest_ind([1.0, 1.5, 2.0], [2.5, 4.0, 5.0], equal_var=False).pvalue
    assert p[1] == pytest.approx(expected)


def gct_archive(tmp_path, body_rows, n_rows=None):
    lines = ["#1.3", f"{n_rows or len(body_rows)}\t3\t1\t1",
             "id\tgeneSymbol\tXA01\tXA02\tXA03", "TMT.Plex\tna\t1\t1\t2", *body_rows]  # fmt: skip
    data = "\n".join(lines).encode()
    path = tmp_path / "freeze.tar.gz"
    with tarfile.open(path, "w:gz") as tar:
        info = tarfile.TarInfo("freeze/prot.gct")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    return path


def test_read_gct(tmp_path):
    """Mutations caught: "NA" read as 0; the column annotation rows read as genes."""
    path = gct_archive(tmp_path, ["NP_1\tESR1\t1.5\tNA\t-0.5", "NP_2\tESR1\tNA\tNA\t1",
                                  "NP_3\tPGR\t0.1\t0.2\t0.3"])  # fmt: skip
    layer, annotation = tables.read_gct(path, "prot.gct")
    assert layer.genes == ("ESR1", "PGR") and layer.samples == ("XA01", "XA02", "XA03")
    assert np.isnan(layer.values[0, 1]) and layer.values[0, 0] == 1.5
    assert annotation["TMT.Plex"].to_list() == ["1", "1", "2"]
    assert layer.collapse is not None and layer.collapse.extra_rows == 1
    with pytest.raises(ProteoError, match="declares"):
        tables.read_gct(gct_archive(tmp_path, ["NP_3\tPGR\t0.1\t0.2\t0.3"], n_rows=2), "prot.gct")


def test_read_copies(tmp_path):
    """Mutations caught: a 0 copy number kept as a value (log2 of 0); a group's columns missed."""
    openpyxl = pytest.importorskip("openpyxl")
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.append(["Majority protein IDs", "Gene names", "CopyNumber_N_01", "CopyNumber_N_02",
                  "CopyNumber_M_01", "CopyNumber_Other_01"])  # fmt: skip
    sheet.append(["P1", "ELANE", 1024, 0, 4, 9])
    sheet.append(["P2", "A;B", 1, 1, 1, 1])
    sheet.append(["P3", "CD14", None, 8, 64, 9])
    path = tmp_path / "copies.xlsx"
    book.save(path)
    layer, n_tested, n_reference = tables.read_copies(path, "_N_", "_M_")
    assert (n_tested, n_reference, layer.zeros) == (2, 1, 1)
    assert layer.genes == ("ELANE", "CD14")
    assert layer.values[0].tolist()[0] == 10 and np.isnan(layer.values[0, 1])
    assert layer.values[1, 2] == 6
