"""U: RSEM expected counts rounded at gather's door, and NCBI Gene IDs placed through HGNC.

The annotation, GENCODE's Entrez table and the sample table are test_gather_counts's. Each
test names the mutation it catches.
"""

import gzip
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from conftest import ctx_run
from test_gather_counts import GENES, SAMPLES, matrix, raw, root_with, whole

from ctx.dataset import CountsSpec, DatasetError, load
from ctx.gather import counts, counts_rsem, counts_run

RSEM_KIND = '\nkind = "rsem_expected"\n'
# Values whose rounding is known: halves go to the even neighbour, as R's round() takes them.
FRACTIONS = np.array([[2.5, 3.5, 0.4, 100.0], [10.0, 20.0, 30.0, 40.0], [7.6, 0.0, 1.0, 2.0]])
ROUNDED = [[2, 4, 0, 100], [10, 20, 30, 40], [8, 0, 1, 2]]
IDS = ["ENSG00000000001.4", "ENSG00000000002.2", "ENSG00000000003.7"]


def declared(tmp_path: Path, extra: str = RSEM_KIND) -> tuple[Path, Path]:
    """test_gather_counts's root, its dataset file's [counts] given ``extra`` lines."""
    root, spec = root_with(tmp_path, entrez=True)
    line = 'annotation = "the test\'s own annotation"'
    spec.write_text(spec.read_text().replace(line + "\n", line + extra))
    return root, spec


def rsem_file(path: Path, rows: list[tuple[str, float]], tpm: float = 999.0) -> Path:
    """An RSEM genes.results file; TPM and FPKM hold ``tpm`` so a misread column shows."""
    head = "gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM"
    lines = [head] + [f"{g}\tT{i}\t1000\t900\t{v:g}\t{tpm:g}\t{tpm:g}"
                      for i, (g, v) in enumerate(rows)]  # fmt: skip
    data = ("\n".join(lines) + "\n").encode()
    path.write_bytes(gzip.compress(data) if path.suffix == ".gz" else data)
    return path


def test_rsem_matrix_rounded_and_stated(tmp_path):
    """A declared matrix of expected counts enters rounded, halves to even, and the page and
    every README state how many values were fractional and the largest change. Catches
    truncation in place of rounding (7.6 to 7), halves rounded up (2.5 to 3), and the
    rounding left out of a README."""
    root, spec = declared(tmp_path)
    path = matrix(tmp_path / "m.tsv", IDS, FRACTIONS)
    done = ctx_run("--root", str(root), "gather", "--dataset", "mat", "--counts", str(path),
                   "--dataset-file", str(spec), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "01_gather"
    (table,) = out.glob("*_gather_mat_counts_v1.parquet")
    assert pl.read_parquet(table).select(SAMPLES).rows() == [tuple(r) for r in ROUNDED]
    said = "4 of 12 values were fractional, the largest change 0.5"
    readmes = sorted(out.glob("*_README.md"))
    assert len(readmes) == 4
    for readme in readmes:
        assert said in readme.read_text(), readme.name
    page = next(out.glob("*_gather_mat_overview_v1.html")).read_text()
    assert said in page
    assert (
        "RSEM expected count for sample S01, rounded"
        in table.with_name(table.stem + "_README.md").read_text()
    )


def test_undeclared_fractions_still_refused():
    """The same fractions without the declaration are refused, naming the declaration. Catches
    the rounding applied to every matrix."""
    found = counts.value_problems(raw(IDS, FRACTIONS))
    assert [k for _, k, _ in found] == ["fractional"]
    assert 'kind = "rsem_expected"' in found[0][2]


@pytest.mark.parametrize(
    ("values", "said"),
    [
        (np.array([[1.25, 7.5, 12.75, 3.0], [9.5, 2.25, 6.0, 0.5]]), "log2 values"),
        (np.array([[400_000.4] * 4, [599_999.6] * 4]), "TPM or CPM"),
        (np.array([[-1.5, 200.0, 3.0, 4.0], [10.4, 0.0, 7.0, 8.0]]), "below 0"),
    ],
)  # fmt: skip
def test_declared_rsem_still_refuses_log_tpm_and_negatives(values, said):
    """Declared expected counts on a log2 or TPM scale, or negative, are refused before any
    rounding. Catches the declaration passing every fraction, and the rounding moved before
    the check (rounded log2 values are whole numbers)."""
    found = counts.value_problems(raw(["A", "B"], values), counts.RSEM_EXPECTED)
    assert any(said in d for _, _, d in found)


def test_check_rounds_after_the_value_check(tmp_path):
    """Through counts_run.check, a declared log2 matrix is refused. Catches the order of
    operations: rounding first makes its values whole and lets it in."""
    root, spec = declared(tmp_path)
    values = np.array([[1.25, 7.5, 12.75, 3.0], [9.5, 2.25, 6.0, 0.5], [4.0, 4.5, 5.5, 1.0]])
    path = matrix(tmp_path / "m.tsv", IDS, values)
    with pytest.raises(counts.CountsError, match="log2 values"):
        counts_run.check(root, load(spec), path)


def test_rsem_folder_combined_by_gene(tmp_path):
    """A folder of genes.results files becomes one matrix joined on gene ID, the sample from
    each file's name by the pattern, expected_count alone read. Catches files joined by row
    order (the second file lists its genes reversed), TPM read as counts, and the pattern
    ignored."""
    root, spec = declared(tmp_path, RSEM_KIND + 'sample_pattern = "^(S\\\\d+)_"\n')
    folder = tmp_path / "rsem"
    folder.mkdir()
    for j, sample in enumerate(SAMPLES):
        rows = [(g, float(FRACTIONS[i, j])) for i, g in enumerate(IDS)]
        rsem_file(folder / f"{sample}_lib.genes.results{'.gz' if j % 2 else ''}",
                  rows[::-1] if j == 1 else rows)  # fmt: skip
    (folder / "notes.txt").write_text("not an RSEM file\n")
    found = counts_run.check(root, load(spec), folder)
    assert found.raw.format == counts_rsem.FORMAT
    assert found.raw.samples == tuple(SAMPLES)
    assert found.mapped.counts.select(SAMPLES).rows() == [tuple(r) for r in ROUNDED]
    assert "TPM" in found.raw.dropped
    assert found.rounding == counts_rsem.Rounding(12, 4, 0.5)


def test_rsem_folder_problems_listed(tmp_path):
    """A file missing a gene, a name the pattern misses, and an RSEM file undeclared are each
    named. Catches a missing gene filled with zero, an unmatched file skipped silently, and
    RSEM read as raw counts without the declaration."""
    root, spec = declared(tmp_path, RSEM_KIND + 'sample_pattern = "^(S\\\\d+)_"\n')
    folder = tmp_path / "rsem"
    folder.mkdir()
    for sample in SAMPLES[:2]:
        rsem_file(folder / f"{sample}_x.genes.results", [(g, 5.0) for g in IDS])
    rsem_file(folder / "S03_x.genes.results", [(g, 5.0) for g in IDS[:2]])
    rsem_file(folder / "other.genes.results", [(g, 5.0) for g in IDS])
    with pytest.raises(counts.CountsError) as caught:
        counts_run.check(root, load(spec), folder)
    assert [(c, k) for c, k, _ in caught.value.problems] == [("rsem", "other.genes.results")]
    (folder / "other.genes.results").unlink()
    with pytest.raises(counts.CountsError) as caught:
        counts_run.check(root, load(spec), folder)
    absent = [d for c, k, d in caught.value.problems if k == "S03_x.genes.results"]
    assert absent and absent[0].startswith("1 genes of S01_x.genes.results are absent")
    (tmp_path / "plain").mkdir()
    root2, spec2 = root_with(tmp_path / "plain")
    one = rsem_file(tmp_path / "S01.genes.results", [(g, 5.0) for g in IDS])
    assert counts_rsem.read(one, CountsSpec("s")).samples == ("S01",)
    with pytest.raises(counts.CountsError) as caught:
        counts_run.check(root2, load(spec2), one)
    assert ("values", "kind") in {(c, k) for c, k, _ in caught.value.problems}


def test_counts_kind_and_pattern_schema(tmp_path):
    """[counts] kind takes raw or rsem_expected; a sample_pattern needs a group. Catches an
    unknown kind read as raw."""
    path = tmp_path / "x.toml"
    for extra, said in (('kind = "tpm"', "counts.kind"),
                        ('sample_pattern = "^S"', "has no group"),
                        ('sample_pattern = "("', "counts.sample_pattern")):  # fmt: skip
        path.write_text(f'[dataset]\nslug = "x"\n\n[counts]\nsource = "s"\n{extra}\n')
        with pytest.raises(DatasetError, match=said):
            load(path)
    path.write_text('[dataset]\nslug = "x"\n\n[counts]\nsource = "s"\n')
    assert load(path).counts == CountsSpec("s", kind="raw")


# HGNC's columns out of order, ensembl before entrez, so reading by position fails.
HGNC = """hgnc_id\tsymbol\tensembl_gene_id\tlocus_type\tentrez_id
HGNC:1\tDELTA\tENSG00000000004\tgene\t104
HGNC:2\tTWINA\tENSG00000000005\tgene\t105
HGNC:3\tTWINB\tENSG00000000006\tgene\t105
HGNC:4\tGONE\tENSG00000000077\tgene\t106
HGNC:5\tALPHA2\tENSG00000000001\tgene\t110
HGNC:6\tETA\tENSG00000000007\tgene\t107
HGNC:7\tIOTA\tENSG00000000009\tgene\t102
HGNC:8\tNOENS\t\tgene\t111
"""


def test_hgnc_places_what_gencode_does_not_link(tmp_path):
    """HGNC places an NCBI Gene ID only where GENCODE's table links none, and only onto one
    annotation gene, never unseating a GENCODE row. Catches HGNC overriding GENCODE (102,
    107), a two-gene HGNC ID kept (105), an Ensembl ID the annotation lacks placed (106), 110
    kept on 101's gene or taking 101 out with it, and the columns read by position."""
    path = tmp_path / "hgnc.txt"
    path.write_text(HGNC)
    hgnc = counts.hgnc_ensembl(path)
    assert hgnc["105"] == {"ENSG00000000005", "ENSG00000000006"}
    assert "111" not in hgnc
    annotation = pl.DataFrame({"gene_id": [g for g, _, _ in GENES],
                               "gene_name": [n for _, n, _ in GENES]})  # fmt: skip
    entrez = {"101": {GENES[0][0]}, "102": {GENES[1][0]},
              "107": {GENES[6][0], GENES[8][0]}}  # fmt: skip
    ids = ["101", "102", "104", "105", "106", "107", "110", "111"]
    mapped = counts.place(raw(ids, whole(len(ids))), annotation, entrez, hgnc)
    assert mapped.genes.rows() == [(GENES[0][0], "ALPHA", "101"), (GENES[1][0], "BETA", "102"),
                                   (GENES[3][0], "DELTA", "104")]  # fmt: skip
    left = {r[0]: r[2] for r in mapped.left_out.rows()}
    assert left == {"105": "ambiguous", "106": "not_in_annotation", "107": "ambiguous",
                    "110": "shared_gene", "111": "not_in_annotation"}  # fmt: skip
    assert mapped.through_hgnc == 1
    without = counts.place(raw(ids, whole(len(ids))), annotation, entrez)
    assert without.through_hgnc == 0
    assert "104" not in without.genes["source_id"].to_list()


def test_hgnc_missing_column_refused(tmp_path):
    """A file without entrez_id is named, not read as empty. Catches a silent zero map."""
    path = tmp_path / "hgnc.txt"
    path.write_text("hgnc_id\tensembl_gene_id\nHGNC:1\tENSG00000000004\n")
    with pytest.raises(counts.CountsError, match="entrez_id"):
        counts.hgnc_ensembl(path)


def test_gather_reads_the_hgnc_resource(tmp_path):
    """With hgnc_complete declared, gather places 104 through it and the page counts it;
    the README names the file. Catches the resource declared and never read."""
    root, spec = root_with(tmp_path, entrez=True)
    (root / "resources" / "h.txt").write_text(HGNC)
    resource = ('\n[[resources]]\nname = "hgnc_complete"\nurls = ["http://127.0.0.1/h.txt"]\n'
                'bytes = 0\nsha256 = "' + "0" * 64 + '"\ndest = "resources/h.txt"\n')  # fmt: skip
    spec.write_text(spec.read_text().replace("[counts]", resource + "\n[counts]"))
    path = matrix(tmp_path / "m.tsv", ["101", "102", "104"], whole(3))
    done = ctx_run("--root", str(root), "gather", "--dataset", "mat", "--counts", str(path),
                   "--dataset-file", str(spec), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "01_gather"
    (genes,) = out.glob("*_gather_mat_genes_v1.parquet")
    assert pl.read_parquet(genes)["source_id"].to_list() == ["101", "102", "104"]
    page = next(out.glob("*_gather_mat_overview_v1.html")).read_text()
    assert "1 NCBI Gene ID rows that GENCODE&#39;s table links to no gene" in page
    assert "HGNC, h.txt" in next(out.glob("*_counts_v1_README.md")).read_text()


def test_ctx_run_check_takes_an_rsem_folder(tmp_path):
    """ctx run --check accepts a folder for --counts. Catches the file-only check."""
    from ctx import pipeline  # noqa: PLC0415 - the module under test

    root, spec = declared(tmp_path)
    folder = tmp_path / "rsem"
    folder.mkdir()
    for sample in SAMPLES:
        rsem_file(folder / f"{sample}.genes.results", [(g, 55.5) for g in IDS])
    problems, _ = pipeline.check(root, pipeline.Options("mat", spec, counts=folder))
    assert not [p.text for p in problems if "--counts" in p.text or p.text.startswith("counts")]


def test_two_hgnc_rows_on_one_gene_both_leave():
    """Two rows HGNC alone places on one gene both leave, as two GENCODE rows do. Catches the
    fallback rule keeping the first of them."""
    annotation = pl.DataFrame({"gene_id": [GENES[3][0]], "gene_name": ["DELTA"]})
    hgnc = {"104": {"ENSG00000000004"}, "204": {"ENSG00000000004"}}
    mapped = counts.place(raw(["104", "204"], whole(2)), annotation, {}, hgnc)
    assert mapped.genes.height == 0
    assert set(mapped.left_out["reason"]) == {"shared_gene"}
