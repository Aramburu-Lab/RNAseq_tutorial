"""M8's TCGA half: target readers and classes, the shift and enrichment, the pairs within a
condition, the hairpin placement, and the page end to end on a made-up GDC root
(tests/ncrna_tcga_root.py). Each test names the mutation it catches, planted and seen to fail."""

import json
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from ctx.dataset import DatasetError, load
from ctx.ncrna import mirna, targets

CTX = Path(sys.executable).parent / "ctx"


@pytest.fixture(scope="module")
def root(tmp_path_factory) -> tuple[Path, Path]:
    """The made-up GDC root and its dataset file."""
    import ncrna_tcga_root  # noqa: PLC0415 - imports synth, which only tests use

    return ncrna_tcga_root.make(tmp_path_factory.mktemp("tcga"), date.today())


@pytest.fixture(scope="module")
def found(root) -> targets.Targets:
    """The made-up target files, read."""
    import ncrna_tcga_root  # noqa: PLC0415

    folder = root[0] / ncrna_tcga_root.TARGETS
    return targets.read({"targetscan_predictions": folder / "predictions.zip",
                         "targetscan_families": folder / "families.zip",
                         "targetscan_genes": folder / "genes.zip",
                         "mirtarbase_strong": folder / "mirtarbase.csv"})  # fmt: skip


def test_readers_keep_human_rows_and_sum_a_genes_sites(found):
    """Catches: a mouse or macaque row kept, a gene's two sites kept as two predictions (or
    one of them dropped), a version left on a gene ID, and miRTarBase's byte-order mark read
    into its first column name."""
    p = found.predicted
    assert set(p["mirna"]) == {"hsa-miR-1-3p", "hsa-miR-2-5p"}
    two = p.filter(pl.col("gene_id") == "g00011")
    assert two.rows() == [("hsa-miR-1-3p", "g00011", "G00011", pytest.approx(-0.3), 2)]
    assert found.names == {"MIMAT9000001": "hsa-miR-1-3p", "MIMAT9000002": "hsa-miR-2-5p"}
    assert "g00010" in found.universe and len(found.universe) == 200
    v = found.validated
    assert v.select("mirna", "gene_name", "references").rows() == [
        ("hsa-miR-1-3p", "G00002", 1),
        ("hsa-miR-1-3p", "G00010", 2),
    ]


def test_validated_first_and_a_weak_prediction_is_no_non_target(found):
    """A gene validated and predicted is validated once; one predicted below the cutoff is
    neither predicted nor a non-target. Catches: validated and predicted merged into one class
    (G00010 counted in both), and non-targets taken as every gene not predicted at the cutoff
    (G00012's weak site would join the non-targets)."""
    tested = pl.DataFrame({"gene_id": [f"g{i:05d}" for i in range(1, 30)],
                           "gene_name": [f"G{i:05d}" for i in range(1, 30)]})  # fmt: skip
    classed = targets.classes(found, "hsa-miR-1-3p", tested, targets.DEFAULT_PRESET)
    by_gene = dict(classed.select("gene_id", "target").rows())
    assert by_gene["g00010"] == targets.VALIDATED
    assert by_gene["g00011"] == targets.PREDICTED
    assert by_gene["g00012"] is None
    assert by_gene["g00002"] == targets.VALIDATED
    assert by_gene["g00003"] == targets.NON_TARGET
    stricter = targets.classes(found, "hsa-miR-1-3p", tested, -0.4)
    assert dict(stricter.select("gene_id", "target").rows())["g00011"] is None
    kept = targets.derived(found, ["hsa-miR-1-3p"], tested)
    both = kept.filter(pl.col("gene_id") == "g00010")
    assert sorted(both["source"]) == ["TargetScan", "miRTarBase"]


def test_a_class_lower_than_the_non_targets_has_a_negative_shift():
    """Catches: the shift taken non-targets minus class, a flipped sign that would turn the
    repression a rising miRNA predicts into a rise."""
    got = targets.shift(np.array([-1.0, -0.5, 0.0]), np.array([0.0, 0.5, 1.0, 1.5]))
    assert got["shift"] == pytest.approx(-0.5 - 0.75)
    assert got["n"] == 3 and 0 < got["d"] <= 1


def test_enrichment_looks_for_targets_among_genes_called_against_the_mirna(found):
    """A rising miRNA's targets are sought among the down genes. Catches: the direction
    flipped (the up genes searched), and a background beyond TargetScan's universe."""
    genes = [f"g{i:05d}" for i in range(1, 260)]
    calls = ["down" if 20 <= i < 60 else "ns" for i in range(1, 260)]
    results = pl.DataFrame({"gene_id": genes, "gene_name": [g.upper() for g in genes],
                            "log2FoldChange": 0.0, "call": calls})  # fmt: skip
    rows = mirna.enrichment(found, {"hsa-miR-1-3p": 2.0}, results)
    predicted = rows.filter(pl.col("target") == targets.PREDICTED).row(0, named=True)
    assert predicted["genes"] == "down" and predicted["hits"] == 40
    assert predicted["expected"] == pytest.approx(41 * 40 / 200)
    flipped = mirna.enrichment(found, {"hsa-miR-1-3p": -2.0}, results)
    assert flipped.filter(pl.col("target") == targets.PREDICTED)["hits"][0] == 0


def test_pairs_correlate_within_each_condition_never_pooled():
    """The tumours' values rise together and so do the normals', but the normals sit at a
    higher miRNA and a lower target, so pooled rho is negative. Catches: rho over all samples
    (negative), which is the tumour-normal gap, not the pair (plan decision 4)."""
    n = 12
    x = np.r_[np.arange(n), np.arange(n) + 100.0]
    y = np.r_[np.arange(n) + 100.0, np.arange(n)]
    samples = pl.DataFrame({"mirna_sample": [f"m{i}" for i in range(2 * n)],
                            "rna_sample": [f"r{i}" for i in range(2 * n)],
                            "case": [f"c{i % n}" for i in range(2 * n)],
                            "condition": ["tumour"] * n + ["normal"] * n})  # fmt: skip
    both = mirna.Layers(samples, x[None, :], ("M1",), y[None, :], ("G1",))
    pairs = pl.DataFrame({"mimat": ["M1"], "mirna": ["miR-x"], "gene_id": ["G1"],
                          "gene_name": ["G1"], "source": ["validated"]})  # fmt: skip
    row = mirna.pair_table(both, pairs).row(0, named=True)
    assert row["rho_tumour"] == pytest.approx(1.0)
    assert row["rho_normal"] == pytest.approx(1.0)
    assert row["n_tumour"] == row["n_normal"] == n
    assert np.corrcoef(x, y)[0, 1] < 0


def test_hairpins_are_placed_one_to_one_on_their_strand():
    """Catches: a hairpin over two genes kept (its STAR count is ambiguous), and a gene on
    the other strand matched."""
    loci = pl.DataFrame({"miRNA_ID": ["h1", "h2", "h3"], "chrom": ["chr1", "chr1", "chr2"],
                         "start": [100, 1000, 100], "end": [180, 1080, 180],
                         "strand": ["+", "+", "+"]})  # fmt: skip
    types = pl.DataFrame({"gene_id": ["a", "b", "c", "d"], "gene_type": "miRNA",
                          "chrom": ["chr1", "chr1", "chr1", "chr2"],
                          "start": [90, 1000, 1050, 90], "end": [200, 1060, 1100, 200],
                          "strand": ["+", "+", "+", "-"]})  # fmt: skip
    assert mirna.placed(loci, types).rows() == [("h1", "a")]


def test_mirnas_without_their_contrast_are_refused(tmp_path):
    """Catches: featured miRNAs accepted with no [ncrna.mirna] to test them in."""
    path = tmp_path / "x.toml"
    path.write_text('[dataset]\nslug = "x"\n\n[[contrasts]]\nname = "c"\nvariable = "v"\n'
                    'reference = "a"\ntested = "b"\nn = [1, 1]\n\n[ncrna]\nfeatured = "c"\n'
                    'pairs = 1\nmirnas = ["hsa-miR-1-3p"]\n')  # fmt: skip
    with pytest.raises(DatasetError, match=r"needs \[ncrna.mirna\]"):
        load(path)


@pytest.fixture(scope="module")
def built(root):
    """ctx ncrna on the made-up GDC root, once for the page tests."""
    import ncrna_tcga_root  # noqa: PLC0415

    path, spec = root
    done = subprocess.run([str(CTX), "--root", str(path), "ncrna", "--dataset",
                           ncrna_tcga_root.DATASET, "--dataset-file", str(spec)],
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    return path / "outputs" / "08_ncrna", done.stdout.split()


def table(out: Path, what: str) -> pl.DataFrame:
    """One table of the run."""
    return pl.read_parquet(next(out.glob(f"*_ncrna_fakegdc_{what}_v1.parquet")))


def test_the_tcga_page_end_to_end(built):
    """The page, its miRNA, target and pair tables and the rules. Catches: the planted miRNA
    not called up, the planted target's anticorrelation lost, a page without the GDC
    acknowledgement, a rule on an unmeasured metric shown with no note, and a page over the
    15 MB cap."""
    out, written = built
    page = Path(written[0])
    assert page.stat().st_size < 15 * 10**6
    html = page.read_text()
    assert "distributed by the NCI Genomic Data Commons" in html
    assert "could not measure: library.histone_pairs" in html
    fit = table(out, "mirna_results").filter(pl.col("gene_name") == "hsa-miR-1-3p")
    assert fit["call"][0] == "up"
    verdicts = dict(table(out, "sanity").select("rule", "rule_verdict").unique().rows())
    assert verdicts == {"up": "PASS", "pair": "PASS", "absent": "PASS", "unmeasured": "FAIL"}
    pairs = table(out, "mirna_target_pairs")
    assert set(pairs["source"]) <= {targets.VALIDATED, targets.PREDICTED}
    planted = pairs.filter(pl.col("gene_name") == "G00010").row(0, named=True)
    assert planted["rho_tumour"] < -0.9 and planted["n_tumour"] == 12
    shifts = table(out, "target_shift")
    assert shifts.height == 2 * (len(targets.PRESETS) + 1)
    blocks = re.findall(r'<script type="application/json" id="ctx-controls-\d+">(.*?)</script>',
                        html, flags=re.S)  # fmt: skip
    kinds = [c["kind"] for c in map(json.loads, blocks)]
    assert kinds.count("pairs") == 2


def test_the_contrast_takes_each_patients_chosen_tumour_alone(built):
    """The root holds a second tumour aliquot of one patient that gather's rule did not take,
    and that sorts first by run accession. Catches: the sheet read without gdc.chosen, where
    rank_replicates would take that aliquot in place of the chosen one."""
    import ncrna_tcga_root  # noqa: PLC0415

    out, _ = built
    shares = table(out, "sample_shares")
    assert ncrna_tcga_root.EXTRA not in set(shares["sample"])
    values = table(out, "filter_values")
    assert values.height > 0


def test_m5_takes_a_gdc_sheet_through_chosen(root, tmp_path):
    """ctx de on the made-up GDC root, the one door its sheet passes. Catches: M5 reading every
    file of the sheet, where rank_replicates drops the patient's chosen tumour for the extra
    aliquot that sorts first; and a GDC sheet refused for its missing cell type."""
    import ncrna_tcga_root  # noqa: PLC0415

    path, spec = root
    done = subprocess.run([str(CTX), "--root", str(path), "de", "--dataset",
                           ncrna_tcga_root.DATASET, "--dataset-file", str(spec),
                           "--no-lfc-shrink"], capture_output=True, text=True,
                          check=False)  # fmt: skip
    out = path / "outputs" / "05_de"
    contrasts = pl.read_parquet(next(out.glob("*_de_fakegdc_contrasts_v1.parquet")))
    assert done.returncode == 0, done.stderr[-3000:]
    row = contrasts.row(0, named=True)
    assert (row["n_reference"], row["n_tested"], row["dropped_runs"]) == (12, 12, "")
