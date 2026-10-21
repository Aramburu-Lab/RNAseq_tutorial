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


@pytest.fixture(scope="session")
def built(ncrna_tcga_run):
    """ctx ncrna on the made-up GDC root, the session's run: (output folder, paths printed)."""
    root, _, printed = ncrna_tcga_run
    return root / "outputs" / "08_ncrna", printed


def table(out: Path, what: str) -> pl.DataFrame:
    """One table of the run."""
    return pl.read_parquet(next(out.glob(f"*_ncrna_fakegdc_{what}_v1.parquet")))


def test_the_tcga_page_end_to_end(built):
    """The page, its miRNA, target and pair tables and the rules. Catches: the planted miRNA
    not called up, the planted target's anticorrelation lost, a page without the GDC
    acknowledgement, a rule on an unmeasured metric called FAIL or shown without its reason, the
    miRNA volcano's live line naming genes, and a page over the 15 MB cap."""
    out, written = built
    page = Path(written[0])
    assert page.stat().st_size < 15 * 10**6
    html = page.read_text()
    assert "were retrieved from the NCI Genomic Data Commons" in html
    assert "holds no read-group record" in html
    fit = table(out, "mirna_results").filter(pl.col("gene_name") == "hsa-miR-1-3p")
    assert fit["call"][0] == "up"
    verdicts = dict(table(out, "sanity").select("rule", "rule_verdict").unique().rows())
    assert verdicts == {"up": "PASS", "pair": "PASS", "absent": "PASS",
                        "unmeasured": "NOT MEASURED"}  # fmt: skip
    pairs = table(out, "mirna_target_pairs")
    assert set(pairs["source"]) <= {targets.VALIDATED, targets.PREDICTED}
    planted = pairs.filter(pl.col("gene_name") == "G00010").row(0, named=True)
    assert planted["rho_tumour"] < -0.9 and planted["n_tumour"] == 12
    shifts = table(out, "target_shift")
    assert shifts.height == 2 * (len(targets.PRESETS) + 1)
    blocks = re.findall(r'<script type="application/json" id="ctx-controls-\d+">(.*?)</script>',
                        html, flags=re.S)  # fmt: skip
    controls = list(map(json.loads, blocks))
    assert [c["kind"] for c in controls].count("pairs") == 2
    (calls,) = [c for c in controls if c["kind"] == "calls"]
    assert calls["unit"] == "mature miRNAs" and "mature miRNAs" in calls["default_summary"]
    assert "genes" not in calls["default_summary"]


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
    # Mutation caught: the M5 page of GDC data without the GDC acknowledgement.
    (page,) = out.glob("*_de_fakegdc_overview_v1.html")
    assert "were retrieved from the NCI Genomic Data Commons" in page.read_text()


TUMOUR = "Primary Tumor"


def f5_files(**changes: list) -> pl.DataFrame:
    """Gather's sheet as the rule marks it: two tumours sequenced both ways, each rRNA-depleted
    twin dropped as another protocol; B also has an FFPE library and A a redacted one, which the
    rule records as protocol, as the reasons' order puts it first."""
    rows = {
        "sample": ["a1", "a2", "a3", "b1", "b2", "b3", "n1", "n2"],
        "sample_id": ["A", "A", "A", "B", "B", "B", "N", "N"],
        "sample_type": [TUMOUR] * 6 + ["Solid Tissue Normal"] * 2,
        "condition": ["tumour"] * 6 + ["normal"] * 2,
        "drop_reason": ["", "protocol", "protocol", "", "protocol", "ffpe", "", "protocol"],
        "preservation": ["OCT", "OCT", "OCT", "OCT", "OCT", "FFPE", "OCT", "OCT"],
        "redacted": [False, False, True, False, False, False, False, False],
        "library_selection": ["Poly-T", "rRNA Depletion", "rRNA Depletion", "Poly-T",
                              "rRNA Depletion", "rRNA Depletion", "Poly-T", "rRNA Depletion"],
    }  # fmt: skip
    return pl.DataFrame(rows | changes)


F5_GENES = pl.DataFrame({"gene_id": ["g1", "g2", "g3"], "gene_name": ["H2BC1", "H2AZ1", "ACTB"],
                         "biotype": ["protein_coding"] * 3})  # fmt: skip
F5_COUNTS = pl.DataFrame({"gene_id": ["g1", "g2", "g3"], "a1": [1, 50, 99], "a2": [10, 0, 90],
                          "a3": [90, 0, 10], "b1": [2, 0, 98], "b2": [20, 0, 80],
                          "b3": [0, 0, 100], "n1": [5, 5, 5], "n2": [9, 1, 5]})  # fmt: skip


def test_f5_reads_the_protocol_twin_and_never_an_ffpe_or_redacted_file():
    """Each rRNA-depleted twin is dropped as another protocol, as gather drops it on every TCGA
    project (LUSC: 9 tumours). Mutations caught: F5 reading only the chosen and tie-break files
    (no pair, as kit v8 found), an FFPE library compared, a redacted file compared because the
    rule recorded protocol first, a normal compared, and the histone share taken over H2AZ1,
    which is replication-independent."""
    from ctx.ncrna import libraries  # noqa: PLC0415

    paired = libraries.pairs(f5_files(), F5_COUNTS, F5_GENES, TUMOUR)
    assert sorted(paired["sample"]) == ["a1", "a2", "b1", "b2"]
    found = libraries.compare(paired)["histone"]
    assert (found["pairs"], found["higher"]) == (2, 2)
    assert found["other"] == pytest.approx(np.median([1 / 150 * 100, 2.0]))


def test_f5_without_a_pair_is_not_measured_with_its_reason():
    """No rRNA-depleted twin, as in BRCA, KIRC and LUAD: the pairs are an empty table and the
    histone rule's metrics are absent with the reason, never a crash (kit v8: IndexError at
    libraries.pairs on all four projects). Mutation caught: the empty selection indexed."""
    from types import SimpleNamespace  # noqa: PLC0415

    from ctx.ncrna import libraries, libraries_page  # noqa: PLC0415

    files = f5_files(library_selection=["Poly-T"] * 8)
    paired = libraries.pairs(files, F5_COUNTS, F5_GENES, TUMOUR)
    assert paired.is_empty() and paired.columns == ["sample_id", "sample", "library_selection",
                                                    *libraries.CLASSES]  # fmt: skip
    tcga = SimpleNamespace(files=files, counts=F5_COUNTS)
    shown, metrics, reasons = libraries_page.group(None, tcga, F5_GENES, TUMOUR)  # type: ignore[arg-type]
    assert shown is None and metrics == {}
    assert set(reasons) == set(libraries_page.RULE_METRICS)
    assert all("no tumour" in why for why in reasons.values())


def test_read_groups_give_each_aliquot_its_library_selection(tmp_path):
    """Catches: a read-group record of another project or with no log row read, and an aliquot
    with two selections kept as one of them."""
    from gdc_fixture import PROJECT, RELEASE  # noqa: PLC0415

    from ctx.gather import gdc  # noqa: PLC0415

    folder = tmp_path / "gdc"
    folder.mkdir()

    def hit(aliquot: str, *selections: str) -> dict:
        """One aligned-reads file hit with its read groups."""
        tree = {
            "samples": [{"portions": [{"analytes": [{"aliquots": [{"submitter_id": aliquot}]}]}]}]
        }
        groups = [{"library_selection": s} for s in selections]
        return {"cases": [tree], "analysis": {"metadata": {"read_groups": groups}}}

    hits = [
        hit("X-1", "Poly-T"),
        hit("X-2", "rRNA Depletion"),
        hit("X-3", "Poly-T", "rRNA Depletion"),
    ]
    (folder / gdc.READ_GROUPS).write_text(json.dumps({"data": {"hits": hits}}))
    query = {"filters": {"op": "in", "content": {"field": "cases.project.project_id",
                                                 "value": [PROJECT]}}}  # fmt: skip
    (folder / f"{gdc.READ_GROUPS}{gdc.QUERY}").write_text(json.dumps(query))
    head = "retrieved_utc\tgdc_release\tfile\tsource\n"
    (folder / gdc.RETRIEVAL).write_text(f"{head}2099\t{RELEASE}\t{gdc.READ_GROUPS}\tmade up\n")
    from gdc_fixture import spec_file  # noqa: PLC0415

    gdc_spec = load(spec_file(tmp_path / "x.toml")).gdc
    assert gdc_spec is not None
    found = gdc.library_selection(folder, gdc_spec)
    assert found is not None
    assert dict(found.rows()) == {"X-1": "Poly-T", "X-2": "rRNA Depletion",
                                  "X-3": "Poly-T;rRNA Depletion"}  # fmt: skip
    (folder / gdc.RETRIEVAL).write_text("retrieved_utc\tgdc_release\tfile\tsource\n")
    with pytest.raises(gdc.GdcError, match="no row for"):
        gdc.library_selection(folder, gdc_spec)
    assert gdc.library_selection(tmp_path, gdc_spec) is None
