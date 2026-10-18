"""R2: a count matrix from any repository, checked at gather's door (``ctx gather --counts``).

Every input is written at test time: a GTF of made-up genes, GENCODE's Entrez table for
them, a matrix and the sample table. Each test names the mutation it catches.
"""

import gzip
import hashlib
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from conftest import ctx_run

from ctx.dataset import DatasetError, load
from ctx.gather import counts, counts_run

# (gene_id, gene_name, transcripts) of the made-up annotation. G5 and G6 share a name, G8
# has two Entrez IDs, and G7 and G9 each take one Entrez ID that also maps to G8's transcript.
GENES = [
    ("ENSG00000000001.4", "ALPHA", ["ENST00000000011.1"]),
    ("ENSG00000000002.2", "BETA", ["ENST00000000021.1"]),
    ("ENSG00000000003.7", "GAMMA", ["ENST00000000031.1", "ENST00000000032.1"]),
    ("ENSG00000000004.1", "DELTA", ["ENST00000000041.1"]),
    ("ENSG00000000005.3", "TWIN", ["ENST00000000051.1"]),
    ("ENSG00000000006.3", "TWIN", ["ENST00000000061.1"]),
    ("ENSG00000000007.1", "ETA", ["ENST00000000071.1"]),
    ("ENSG00000000008.1", "THETA", ["ENST00000000081.1", "ENST00000000082.1"]),
    ("ENSG00000000009.2", "IOTA", ["ENST00000000091.1"]),
]
ENTREZ = [
    ("ENST00000000011.1", "101"),
    ("ENST00000000021.1", "102"),
    ("ENST00000000031.1", "103"),
    ("ENST00000000032.1", "103"),
    ("ENST00000000081.1", "108"),
    ("ENST00000000082.1", "109"),
    ("ENST00000000071.1", "107"),
    ("ENST00000000091.1", "107"),
    ("ENST99999999999.1", "999"),
]
SAMPLES = ["S01", "S02", "S03", "S04"]
TABLE = "sample,patient,tissue\nS01,p1,tumour\nS02,p1,normal\nS03,p2,tumour\nS04,p2,normal\n"
SPEC = """[dataset]
slug = "mat"

[guard]
ids = []
patterns = []

[[resources]]
name = "annotation_gtf"
urls = ["http://127.0.0.1/a.gtf.gz"]
bytes = 0
sha256 = "0000000000000000000000000000000000000000000000000000000000000000"
dest = "resources/a.gtf.gz"
{entrez}
[counts]
source = "made-up matrix for a test"
annotation = "the test's own annotation"

[sdrf]
file = "mat.csv"
join = ["sample", "sample"]

[columns.donor]
source = "patient"
meaning = "Patient."

[columns.cell_type]
source = "tissue"
meaning = "Tissue."
values = {{ "tumour" = "larynx", "normal" = "larynx" }}

[columns.condition]
source = "tissue"
meaning = "Tumour or normal."
"""
ENTREZ_RESOURCE = """
[[resources]]
name = "gencode_entrez"
urls = ["http://127.0.0.1/e.gz"]
bytes = 0
sha256 = "0000000000000000000000000000000000000000000000000000000000000000"
dest = "resources/e.gz"
"""


def gtf_text() -> str:
    """A GTF of ``GENES``: a header naming the release, then gene and transcript lines."""
    lines = ["##description: made-up annotation, version 1 (test)", "##provider: test"]
    for i, (gene, name, transcripts) in enumerate(GENES):
        where = f"chr1\ttest\t{{}}\t{100 * i + 1}\t{100 * i + 90}\t.\t+\t."
        attrs = f'gene_id "{gene}"; gene_type "protein_coding"; gene_name "{name}";'
        lines.append(f"{where.format('gene')}\t{attrs}")
        lines += [f'{where.format("transcript")}\t{attrs} transcript_id "{t}";'
                  for t in transcripts]  # fmt: skip
    return "\n".join(lines) + "\n"


def root_with(tmp_path: Path, entrez: bool = False) -> tuple[Path, Path]:
    """A root with the annotation fetched, and the dataset file beside its table."""
    root = tmp_path / "root"
    (root / "resources").mkdir(parents=True)
    (root / "resources" / "a.gtf.gz").write_bytes(gzip.compress(gtf_text().encode()))
    rows = "".join(f"{t}\t{e}\n" for t, e in ENTREZ)
    (root / "resources" / "e.gz").write_bytes(gzip.compress(rows.encode()))
    conf = tmp_path / "conf"
    conf.mkdir()
    (conf / "mat.csv").write_text(TABLE)
    spec = conf / "mat.toml"
    spec.write_text(SPEC.format(entrez=ENTREZ_RESOURCE if entrez else ""))
    return root, spec


def matrix(path: Path, ids: list[str], values: np.ndarray, sep: str = "\t",
           head: str = "gene", samples: list[str] = SAMPLES) -> Path:  # fmt: skip
    """A matrix file: a header, then one row per ID; gzipped when the name ends .gz."""
    lines = [sep.join([head, *samples])]
    lines += [sep.join([i, *(f"{v:g}" for v in row)]) for i, row in zip(ids, values, strict=True)]
    data = ("\n".join(lines) + "\n").encode()
    path.write_bytes(gzip.compress(data) if path.suffix == ".gz" else data)
    return path


def whole(n: int, seed: int = 1) -> np.ndarray:
    """n rows of whole counts, four samples."""
    return np.random.default_rng(seed).integers(0, 500, (n, len(SAMPLES))).astype(float)


def test_gather_counts_end_to_end(tmp_path):
    """ctx gather --counts on a gzipped TSV with versioned Ensembl IDs: the sheet, counts and
    genes tables and their READMEs. Catches counts reordered against their genes, a version
    compared as equal when it differs, and the matrix's SHA-256 left out of the README."""
    root, spec = root_with(tmp_path)
    ids = ["ENSG00000000001.4", "ENSG00000000002.9", "ENSG00000000003", "ENSG00000000077.1"]
    values = whole(4)
    path = matrix(tmp_path / "m.tsv.gz", ids, values)
    done = ctx_run("--root", str(root), "gather", "--dataset", "mat", "--counts", str(path),
                   "--dataset-file", str(spec), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    out = root / "outputs" / "01_gather"
    (sheet,) = out.glob("*_gather_mat_samplesheet_v1.parquet")
    frame = pl.read_parquet(sheet)
    assert frame["sample"].to_list() == SAMPLES
    assert frame["counts_file"].unique().to_list() == ["m.tsv.gz"]
    assert frame["condition"].to_list() == ["tumour", "normal", "tumour", "normal"]
    assert counts.is_sheet(frame)
    (table,) = out.glob("*_gather_mat_counts_v1.parquet")
    got = pl.read_parquet(table)
    assert got["gene_id"].to_list() == [g for g, _, _ in GENES[:3]]
    assert got.select(SAMPLES).to_numpy().tolist() == values[:3].astype(int).tolist()
    assert got.dtypes[1] == pl.Int32
    (genes,) = out.glob("*_gather_mat_genes_v1.parquet")
    assert pl.read_parquet(genes)["source_id"].to_list() == ids[:3]
    (left,) = out.glob("*_gather_mat_genes_left_out_v1.parquet")
    assert pl.read_parquet(left).rows() == [
        ("ENSG00000000077.1", "ensembl", "not_in_annotation", "", float(values[3].sum()))
    ]
    readme = table.with_name(table.stem + "_README.md").read_text()
    assert hashlib.sha256(path.read_bytes()).hexdigest() in readme
    assert "made-up matrix for a test" in readme
    assert "made-up annotation, version 1 (test)" in readme
    page = next(out.glob("*_gather_mat_overview_v1.html")).read_text()
    assert "1 equal the annotation&#39;s, 1 differ, 1 without a version" in page


def raw(ids: list[str], values: np.ndarray, samples: tuple[str, ...] = tuple(SAMPLES)):
    """A Raw as ``counts.read`` returns it."""
    return counts.Raw(tuple(ids), samples, values, (), (), "table", ())


@pytest.mark.parametrize(
    ("values", "check", "said"),
    [
        (np.array([[1.5, 2.25, 3.0, 4.0], [310.0, 0.0, 7.0, 8.0]]), "fractional",
         "FPKM, RPKM or estimated counts"),
        (np.array([[-1.2, 2.5, 3.1, 4.0], [10.4, 0.0, 7.0, 8.0]]), "negative", "log2"),
        (np.array([[400_000.0] * 4, [600_000.0] * 4]), "scaled", "rounded TPM or CPM"),
        (np.array([[400_000.5] * 4, [599_999.5] * 4]), "fractional", "TPM or CPM"),
        (np.array([[np.nan, 2.0, 3.0, 4.0], [1.0, 0.0, 7.0, 8.0]]), "missing", "never missing"),
    ],
)  # fmt: skip
def test_value_problems_name_the_scale(values, check, said):
    """Fractions, negatives, a million per sample and empty cells each refuse, naming how it
    told. Catches the sum rule dropped (rounded TPM passes as counts), the negative check
    reading only fractions, and NaN passing as a whole number."""
    found = counts.value_problems(raw(["A", "B"], values))
    assert check in [k for _, k, _ in found]
    assert any(said in d for _, _, d in found)


def test_scaled_edge_is_one_percent():
    """A sum exactly 1 % off a million is scaled; just past it is not. Catches the tolerance
    compared with < instead of <=, and a tolerance of 10 %."""
    on = np.array([[990_000.0] * 4])
    past = np.array([[989_999.0] * 4])
    assert [k for _, k, _ in counts.value_problems(raw(["A"], on))] == ["scaled"]
    assert counts.value_problems(raw(["A"], past)) == []
    assert counts.value_problems(raw(["A"], np.array([[900_000.0] * 4]))) == []


def test_raw_counts_pass():
    """Whole counts with ordinary library sizes raise nothing."""
    assert counts.value_problems(raw(["A", "B"], whole(2))) == []


def test_every_problem_listed_together(tmp_path):
    """A matrix with fractions, a repeated gene (by its stable ID), a repeated sample, and a
    sample the table lacks: one run names all four, writes the problem table and no sheet.
    Catches a check that stops the run before the others, and the version kept when looking
    for repeats."""
    root, spec = root_with(tmp_path)
    ids = ["ENSG00000000001.4", "ENSG00000000001.5", "ENSG00000000002.2"]
    values = whole(3).copy()
    values[0, 0] = 2.5
    path = matrix(tmp_path / "m.tsv", ids, np.column_stack([values, values[:, :1]]),
                  samples=[*SAMPLES, "S01"])  # fmt: skip
    path2 = matrix(tmp_path / "m2.tsv", ids, np.column_stack([values, values[:, :1]]),
                   samples=[*SAMPLES, "S05"])  # fmt: skip
    with pytest.raises(counts.CountsError) as caught:
        counts_run.check(root, load(spec), path)
    found = {(c, k) for c, k, _ in caught.value.problems}
    assert {("values", "fractional"), ("duplicate gene", "ENSG00000000001"),
            ("duplicate sample", "S01")} <= found  # fmt: skip
    with pytest.raises(counts.CountsError) as caught:
        counts_run.check(root, load(spec), path2)
    found = {(c, k) for c, k, _ in caught.value.problems}
    assert ("join", "S05") in found
    assert ("values", "fractional") in found
    done = ctx_run("--root", str(root), "gather", "--dataset", "mat", "--counts", str(path2),
                   "--dataset-file", str(spec), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 1
    out = root / "outputs" / "01_gather"
    assert list(out.glob("*_counts_problems_v1.csv"))
    assert not list(out.glob("*_samplesheet_v1.*"))


def test_a_table_sample_missing_from_the_matrix(tmp_path):
    """A table row with no matrix column is named. Catches the join checked one way only."""
    root, spec = root_with(tmp_path)
    path = matrix(tmp_path / "m.tsv", ["ALPHA"], whole(1)[:, :3], samples=SAMPLES[:3])
    with pytest.raises(counts.CountsError) as caught:
        counts_run.check(root, load(spec), path)
    assert ("join", "S04") in {(c, k) for c, k, _ in caught.value.problems}


def test_ids_placed_by_kind():
    """Symbols, Entrez IDs and Ensembl IDs, each by its own kind; ambiguous and shared IDs are
    left out with the genes they hit, never guessed. Catches a symbol taken by its first gene,
    an Entrez ID with two genes kept, two rows summed onto one gene, and a row left out
    without its counts."""
    annotation = pl.DataFrame({"gene_id": [g for g, _, _ in GENES],
                               "gene_name": [n for _, n, _ in GENES]})  # fmt: skip
    entrez = {"101": {GENES[0][0]}, "103": {GENES[2][0]}, "107": {GENES[6][0], GENES[8][0]},
              "108": {GENES[7][0]}, "109": {GENES[7][0]}}  # fmt: skip
    ids = ["BETA", "TWIN", "NOPE", "101", "103", "107", "108", "109", "ENSG00000000004"]
    values = whole(len(ids))
    mapped = counts.place(raw(ids, values), annotation, entrez)
    assert mapped.genes.rows() == [
        (GENES[1][0], "BETA", "BETA"), (GENES[0][0], "ALPHA", "101"),
        (GENES[2][0], "GAMMA", "103"), (GENES[3][0], "DELTA", "ENSG00000000004"),
    ]  # fmt: skip
    left = {r[0]: r[2] for r in mapped.left_out.rows()}
    assert left == {"TWIN": "ambiguous", "NOPE": "not_in_annotation", "107": "ambiguous",
                    "108": "shared_gene", "109": "shared_gene"}  # fmt: skip
    totals = dict(mapped.left_out.select("source_id", "total").rows())
    assert totals["108"] == values[6].sum()
    assert mapped.kinds == {"symbol": 3, "entrez": 5, "ensembl": 1}
    assert mapped.versions == {"equal": 0, "differ": 0, "unversioned": 1}


def test_entrez_through_the_gtf(tmp_path):
    """GENCODE's table read through the GTF's transcripts; a transcript the GTF lacks is
    skipped. Catches the gene read from the Entrez table's own column, and the unknown
    transcript kept as a gene."""
    root, _ = root_with(tmp_path)
    found = counts.entrez_genes(root / "resources" / "e.gz", root / "resources" / "a.gtf.gz")
    assert found["103"] == {GENES[2][0]}
    assert found["107"] == {GENES[6][0], GENES[8][0]}
    assert "999" not in found


def test_entrez_ids_need_the_table(tmp_path):
    """A matrix of Entrez IDs with no gencode_entrez resource stops with the fix named, and
    with the resource gather places them. Catches every row left out silently instead."""
    root, spec = root_with(tmp_path)
    path = matrix(tmp_path / "m.tsv", ["101", "102"], whole(2))
    with pytest.raises(counts.CountsError, match="add the gencode_entrez resource"):
        counts_run.check(root, load(spec), path)
    (tmp_path / "with").mkdir()
    root2, spec2 = root_with(tmp_path / "with", entrez=True)
    found = counts_run.check(root2, load(spec2), path)
    assert found.mapped.genes["gene_name"].to_list() == ["ALPHA", "BETA"]


def test_featurecounts_and_summary_rows(tmp_path):
    """featureCounts' comment line and five annotation columns, and HTSeq's and STAR's summary
    rows, leave before any check, and are named. Catches Length read as a sample and
    __no_feature counted as a gene."""
    lines = ["# Program:featureCounts v2.0.6; Command:...",
             "\t".join(["Geneid", *counts.FEATURECOUNTS, *SAMPLES]),
             "\t".join(["BETA", "chr1", "1", "90", "+", "90", "1", "2", "3", "4"]),
             "\t".join(["__no_feature", "", "", "", "", "", "9", "9", "9", "9"]),
             "\t".join(["N_unmapped", "", "", "", "", "", "5", "5", "5", "5"])]  # fmt: skip
    path = tmp_path / "fc.txt"
    path.write_text("\n".join(lines) + "\n")
    found = counts.read(path)
    assert found.format == "featureCounts"
    assert found.samples == tuple(SAMPLES)
    assert found.ids == ("BETA",)
    assert found.values.tolist() == [[1.0, 2.0, 3.0, 4.0]]
    assert found.summary == (("__no_feature", 36.0), ("N_unmapped", 20.0))


def test_csv_by_ending_and_text_cells(tmp_path):
    """A .csv is read with commas; the same text as .txt is refused with how it was read; a
    word in a count cell is named with its column. Catches tabs forced on a CSV, and a text
    cell read as missing."""
    body = "gene,S01,S02\nBETA,1,2\nALPHA,3,NA\n"
    (tmp_path / "m.csv").write_text(body)
    (tmp_path / "m.txt").write_text(body)
    found = counts.read(tmp_path / "m.csv")
    assert found.samples == ("S01", "S02")
    assert found.problems == (("values", "S02", "1 cells are not numbers, the first 'NA'"),)
    with pytest.raises(counts.CountsError, match="read with tabs"):
        counts.read(tmp_path / "m.txt")


def test_read_level_modules_refuse(tmp_path):
    """process, quant and qc --reads refuse a count matrix's sheet with the reason. Catches
    a module that reads the sheet as a fetchngs one and fails on its missing FASTQ paths."""
    root, spec = root_with(tmp_path)
    path = matrix(tmp_path / "m.tsv", ["ALPHA", "BETA"], whole(2))
    done = ctx_run("--root", str(root), "gather", "--dataset", "mat", "--counts", str(path),
                   "--dataset-file", str(spec), cwd=tmp_path)  # fmt: skip
    assert done.returncode == 0, done.stderr[-2000:]
    for command in (
        ["process", "--genome", str(tmp_path)],
        ["quant", "--rnaseq-dir", str(tmp_path)],
        ["qc", "--reads"],
    ):
        done = ctx_run("--root", str(root), command[0], "--dataset", "mat", "--dataset-file",
                       str(spec), *command[1:], cwd=tmp_path)  # fmt: skip
        assert done.returncode != 0, command
        assert "count matrix's sheet" in done.stderr + done.stdout, (command, done.stderr[-800:])


def test_counts_table_schema(tmp_path):
    """[counts] needs its source, and is refused beside [gdc]."""
    path = tmp_path / "x.toml"
    path.write_text('[dataset]\nslug = "x"\n\n[counts]\nurl = "u"\n')
    with pytest.raises(DatasetError, match="counts.source"):
        load(path)
    path.write_text('[dataset]\nslug = "x"\n\n[counts]\nsource = "s"\nfile = "f"\n')
    with pytest.raises(DatasetError, match="counts.file: unknown key"):
        load(path)


DE_SPEC = """[dataset]
slug = "matde"

[guard]
ids = []
patterns = []

[[resources]]
name = "annotation_gtf"
urls = ["http://127.0.0.1/a.gtf.gz"]
bytes = 0
sha256 = "0000000000000000000000000000000000000000000000000000000000000000"
dest = "resources/a.gtf.gz"

[counts]
source = "tests/synth.py counts, seed 42"

[sdrf]
file = "matde.tsv"
join = ["sample", "sample"]

[columns.donor]
source = "donor"
meaning = "Donor."

[columns.cell_type]
source = "cell_type"
meaning = "Cell type."

[columns.condition]
source = "condition"
meaning = "Condition."

[[contrasts]]
name = "a_case"
variable = "condition"
reference = "control"
tested = "case"
n = [6, 6]
filters = { cell_type = ["type_a"] }

[de]
featured = "a_case"
"""


def test_gather_qc_de_run_on_a_matrix(tmp_path):
    """gather --counts, then qc and de unchanged on its tables: de's calls recover the planted
    condition genes. Catches a counts table whose columns are not the sheet's samples, and a
    sheet qc or de cannot read."""
    from synth import counts as synthetic  # noqa: PLC0415 - tests/ helper

    made = synthetic(n_genes=400)
    ids = [f"ENSG{i + 1:011d}.1" for i in range(made.counts.height)]
    root = tmp_path / "root"
    (root / "resources").mkdir(parents=True)
    lines = [f'chr1\tt\tgene\t{i + 1}\t{i + 9}\t.\t+\t.\tgene_id "{g}"; gene_type "protein_coding";'
             f' gene_name "N{i}";' for i, g in enumerate(ids)]  # fmt: skip
    (root / "resources" / "a.gtf.gz").write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
    conf = tmp_path / "conf"
    conf.mkdir()
    made.samples.select("sample", "donor", "cell_type", "condition").write_csv(
        conf / "matde.tsv", separator="\t"
    )
    (conf / "matde.toml").write_text(DE_SPEC)
    path = matrix(tmp_path / "m.csv", ids, made.counts.drop("gene_id").to_numpy().astype(float),
                  sep=",", samples=made.counts.columns[1:])  # fmt: skip
    base = ("--dataset", "matde", "--dataset-file", str(conf / "matde.toml"))
    for command in (("gather", "--counts", str(path)), ("qc",), ("de",)):
        done = ctx_run("--root", str(root), command[0], *base, *command[1:], cwd=tmp_path)
        assert done.returncode == 0, (command, done.stderr[-3000:])
    (results,) = (root / "outputs").glob("*_de/*_de_matde_a_case_results_v1.parquet")
    called = pl.read_parquet(results).filter(pl.col("call") != "ns")["gene_id"]
    planted = made.truth.with_columns(gene_id=pl.Series(ids)).filter(
        pl.col("lfc_condition_type_a") != 0
    )["gene_id"]
    assert len(set(called) & set(planted)) >= 0.8 * planted.len()


def test_ctx_run_plans_the_matrix_route(tmp_path):
    """ctx run --counts gives gather the matrix and no fetchngs outdir, runs qc and de on its
    counts, never quant, and --check lists a matrix problem with the others. Catches the
    matrix dropped from gather's arguments, and de left out as if no counts existed."""
    from ctx import pipeline  # noqa: PLC0415 - the module under test

    root, spec = root_with(tmp_path)
    contrast = (
        '\n[[contrasts]]\nname = "tumour_normal"\nvariable = "condition"\nreference = "normal"\n'
        'tested = "tumour"\nn = [2, 2]\n\n[de]\nfeatured = "tumour_normal"\n'
    )
    spec.write_text(spec.read_text() + contrast)
    path = matrix(tmp_path / "m.tsv", ["ALPHA", "BETA"], whole(2))
    options = pipeline.Options("mat", spec, counts=path)
    planned = dict(pipeline.steps(load(spec), options))
    assert planned["gather"][-2:] == ("--counts", str(path))
    assert "--fetchngs-dir" not in planned["gather"]
    assert {"gather", "qc", "de", "report"} <= set(planned)
    assert "quant" not in planned
    problems, _ = pipeline.check(root, options)
    assert [p.text for p in problems if p.blocks] == []
    values = whole(2).copy()
    values[0, 0] = 0.5
    bad = matrix(tmp_path / "bad.tsv", ["ALPHA", "BETA"], values)
    problems, _ = pipeline.check(root, pipeline.Options("mat", spec, counts=bad))
    assert any(p.text.startswith("counts: values: fractional") for p in problems)
