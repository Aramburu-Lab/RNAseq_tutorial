"""R3: a protein table of any study, its adapters, the identifier mapping, and M7 on it.

Unit tests on hand-built tables, then ``ctx proteo`` end to end on tests/proteo_study_root.py,
once per adapter. Each test names the mutation it catches. No real accession is typed here.
"""

import math
import subprocess
import zipfile
from dataclasses import replace

import numpy as np
import polars as pl
import pytest
from clinical_root import CTX
from proteo_study_root import SAMPLES, SLUG, build, run

from ctx.dataset import DatasetError, ProteoTable, load
from ctx.proteo import ProteoError, adapters, identifiers, quantities

RNA = pl.DataFrame(
    {
        "gene_id": ["ENSG01.3", "ENSG02.1", "ENSG03.2", "ENSG04.1", "ENSG04.1_PAR_Y", "ENSG05.1"],
        "gene_name": ["TP53", "MYC", "DUP", "DUP", "XY", None],
    }
)
DECLARED = ProteoTable("table", "s.csv", "linear", ("NA", "0"), "symbol", "tmt", "plex")


def samples(names: list[str]) -> quantities.SampleTable:
    """A sample table naming these protein columns, each with an RNA sample, one plex."""
    return quantities.SampleTable(
        tuple(names),
        tuple(f"r{n}" for n in names),
        ("p1",) * len(names),
        pl.DataFrame({"sample": names}),
        "x",
    )


def text(rows: dict[str, list[str]]) -> quantities.Text:
    """A reader's text columns."""
    return quantities.Text(pl.DataFrame(rows), (), len(rows) - 1)


def test_mapping_reasons_in_order():
    """Each reason once. Mutations caught: a several-gene group resolved to its first gene; a
    symbol two RNA genes share taken as either; an isoform or FASTA-wrapped accession unmapped;
    a version left on an Ensembl ID; a group with an unmapped member mapped by the other."""
    cells = ["TP53", "TP53;MYC", "DUP", "NOPE", "", "MYC;MYC"]
    found = identifiers.to_genes(cells, "symbol", RNA)
    assert found.reasons == ("mapped", "several", "shared", "absent", "no_id", "mapped")
    assert found.genes == ("TP53", None, None, None, None, "MYC")
    assert found.several == (("TP53;MYC", "MYC;TP53"),)
    up = {"P04637": "TP53", "P01106": "MYC"}
    cells = ["P04637-2", "sp|P01106|MYC_HUMAN", "P04637;Q99999", "Q99999"]
    found = identifiers.to_genes(cells, "uniprot", RNA, up)
    assert found.reasons == ("mapped", "mapped", "no_map", "no_map")
    found = identifiers.to_genes(["ENSG01.7", "ENSG04", "ENSG04.1_PAR_Y", "ENSG09"], "ensembl", RNA)
    assert found.genes == ("TP53", "ENSG04", "XY", None)
    with pytest.raises(ValueError, match="map"):
        identifiers.to_genes(["P1"], "uniprot", RNA)


def test_missing_never_zero_and_undeclared_zero_refused():
    """Mutations caught: a declared 0 read as a value (log2 of 0); an undeclared 0 or text read
    as missing; a log10 table not turned into log2; problems raised one at a time."""
    read = quantities.read(
        text({"id": ["TP53", "MYC"], "a": ["8", "0"], "b": ["NA", "16"]}),
        samples(["a", "b"]),
        DECLARED,
        RNA,
        "id",
    )
    assert read.layer.genes == ("TP53", "MYC")
    assert np.array_equal(read.layer.values, [[3.0, np.nan], [np.nan, 4.0]], equal_nan=True)
    assert read.tokens == {"NA": 1, "0": 1}
    strict = replace(DECLARED, missing=("NA",))
    with pytest.raises(ProteoError) as error:
        quantities.read(
            text({"id": ["TP53", "MYC"], "a": ["0", "x"], "b": ["-1", "4"]}),
            samples(["a", "b"]),
            strict,
            RNA,
            "id",
        )
    message = str(error.value)
    assert "'a': 1 cells neither" in message and "'a': 1 intensities at or below 0" in message
    assert "'b': 1 intensities" in message
    logged = replace(DECLARED, scale="log", log_base=10, missing=("NA",))
    read = quantities.read(
        text({"id": ["TP53"], "a": ["2"], "b": ["0"]}), samples(["a", "b"]), logged, RNA, "id"
    )
    assert (
        read.layer.values[0, 0] == pytest.approx(2 * math.log2(10)) and read.layer.values[0, 1] == 0
    )


def test_collapse_after_mapping_keeps_the_best_quantified_row():
    """Two rows of one gene: the one quantified in more samples stays. Mutations caught: the
    first row kept; rows collapsed by their identifier text before mapping (TP53 and its
    isoform would be two genes)."""
    up = {"P04637": "TP53"}
    rows = {"id": ["P04637", "P04637-2"], "a": ["2", "4"], "b": ["NA", "8"]}
    read = quantities.read(
        text(rows), samples(["a", "b"]), replace(DECLARED, id_kind="uniprot"), RNA, "id", up
    )
    assert read.layer.genes == ("TP53",) and list(read.layer.values[0]) == [2.0, 3.0]
    assert (read.layer.collapse.extra_rows, read.layer.collapse.kept) == (1, 1)


def test_reader_lists_absent_columns_together():
    """Mutations caught: the first absent column alone reported; an absent identifier column
    reported after a value problem."""
    with pytest.raises(ProteoError) as error:
        quantities.read(
            text({"x": ["TP53"], "a": ["2"]}), samples(["a", "b", "c"]), DECLARED, RNA, "id"
        )
    assert "no identifier column 'id'" in str(error.value)
    assert "2 sample table protein samples" in str(error.value)


def test_sample_table_problems_together(tmp_path):
    """Mutations caught: a protein or RNA sample on two rows accepted; a row with no batch
    accepted; a missing column reported alone after the others."""
    path = tmp_path / "s.csv"
    path.write_text("protein_sample,rna_sample,plex\na,r1,p1\na,r1,\nb,,p2\n")
    problems: list[str] = []
    found = quantities.sample_table(path, DECLARED, problems)
    assert found is not None and found.rna == ("r1", "r1", "")
    assert len(problems) == 3 and "line 3 lacks" in problems[2], problems
    path.write_text("protein_sample,plex\na,p1\n")
    problems = []
    assert quantities.sample_table(path, DECLARED, problems) is None
    assert "no column ['rna_sample']" in problems[0]


def maxquant_bytes() -> bytes:
    """A proteinGroups.txt with each flag planted, and one row carrying two flags."""
    head = "\t".join(
        [
            "Gene names",
            "LFQ intensity s1",
            "LFQ intensity s2",
            "Intensity s1",
            "Reverse",
            "Contaminant",
            "Only identified by site",
        ]
    )
    rows = [
        "A\t1\t2\t9\t\t\t",
        "B\t1\t2\t9\t+\t\t",
        "C\t1\t2\t9\t+\t+\t",
        "D\t1\t2\t9\t\t+\t",
        "E\t1\t2\t9\t\t\t+",
        "F\t0\t2\t9\t\t+\t+",
    ]
    return ("\n".join([head, *rows]) + "\n").encode()


def test_maxquant_drops_each_flag_once_in_order():
    """Mutations caught: a row with two flags counted twice; the pre-1.6 Contaminant column
    ignored; Intensity columns taken for LFQ; a flagged row kept."""
    found = adapters.maxquant(maxquant_bytes(), "Gene names", "LFQ intensity")
    assert found.dropped == (("reverse", 2), ("contaminant", 2), ("only_by_site", 1))
    assert found.frame.columns == ["Gene names", "s1", "s2"] and found.columns == 2
    assert found.frame["Gene names"].to_list() == ["A"]
    with pytest.raises(ProteoError, match="no column starts"):
        adapters.maxquant(maxquant_bytes(), "Gene names", "iBAQ")


def test_diann_names_runs_by_file_and_refuses_twins():
    """Mutations caught: a run named by its whole path; two runs of one file name accepted."""
    head = "Protein.Group\tGenes\t/data/a.raw\tD:\\x\\b.d\n"
    found = adapters.diann((head + "P1\tTP53\t5\t\n").encode(), "Genes")
    assert found.frame.columns == ["Genes", "a", "b"] and found.frame["b"][0] == ""
    with pytest.raises(ProteoError, match="share a file name"):
        adapters.diann(b"Genes\t/x/a.raw\t/y/a.raw\n", "Genes")


def test_archive_member_from_zip(tmp_path):
    """Catches a zip read as a tar: PRIDE serves MaxQuant's txt folder zipped."""
    path = tmp_path / "txt.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("combined/txt/proteinGroups.txt", maxquant_bytes())
    assert quantities.raw_bytes(path, "proteinGroups.txt") == maxquant_bytes()
    with pytest.raises(ProteoError, match="0 members"):
        quantities.raw_bytes(path, "evidence.txt")


def test_dataset_table_problems_together(tmp_path):
    """Mutations caught: an unknown format, scale or kind accepted; UniProt without a map
    accepted; a plain table without its identifier column accepted; table.file beside a
    protein resource accepted; one problem reported at a time."""
    path = tmp_path / "d.toml"
    path.write_text(
        '[dataset]\nslug = "d"\n[proteo]\ndesign = "across_samples"\n'
        'protein = "x"\n[proteo.table]\nformat = "excel"\nfile = "t.tsv"\n'
        'samples = "s.xlsx"\nscale = "ln"\nid_kind = "uniprot"\n'
        'acquisition = "tmt"\nbatch = "plex"\n'
    )
    with pytest.raises(DatasetError) as error:
        load(path)
    message = str(error.value)
    for part in (
        "format: 'excel'",
        "samples: 's.xlsx'",
        "scale: 'ln'",
        "id_map",
        "protein (a fetched resource) or table.file",
        "proteo.table.missing",
    ):
        assert part in message, part
    path.write_text(path.read_text().replace('format = "excel"', 'format = "table"'))
    with pytest.raises(DatasetError, match="id_column"):
        load(path)


@pytest.fixture(scope="module", params=["maxquant", "diann"])
def study_root(request, tmp_path_factory):
    """ctx proteo run on the made-up study, once per adapter."""
    return request.param, run(tmp_path_factory.mktemp(f"r3_{request.param}"), request.param)


def output(root, what: str):
    """The one proteo output of this kind."""
    (path,) = (root / "outputs" / "07_proteo").glob(f"*_proteo_{SLUG}_{what}")
    return path


def test_study_end_to_end(study_root):
    """The page is built and the planted rules pass. Mutations caught: the reverse decoy kept
    (G9's rho falls); G10's missing values read as zero (tested); the unpaired protein sample
    correlated; the several-genes group dropped unlisted; the batch design unstated."""
    fmt, root = study_root
    genes = pl.read_parquet(output(root, "genes_v1.parquet"))
    g9 = genes.filter(pl.col("gene") == "G9").row(0, named=True)
    assert g9["rho"] > 0.95 and g9["n_quantified"] == SAMPLES - 1
    g10 = genes.filter(pl.col("gene") == "G10").row(0, named=True)
    assert g10["n_quantified"] == SAMPLES // 2 and not g10["tested"]
    verdicts = pl.read_parquet(output(root, "sanity_v1.parquet"))
    assert set(verdicts["rule_verdict"]) == {"PASS"}, verdicts.select("rule", "left").rows()
    html = output(root, "overview_v1.html").read_text()
    assert f"{SAMPLES} of {SAMPLES + 1} protein samples are paired" in html
    assert "TMT design: 3 plexes of 8 to 9 samples" in html
    batches = pl.read_csv(output(root, "batches_v1.csv"), comment_prefix="#")
    assert batches["samples"].sum() == SAMPLES + 1 and batches["paired"].sum() == SAMPLES
    if fmt == "maxquant":
        assert "2 reverse decoys, 1 contaminants, 0 identified only" in html
        assert "12 rows map to one RNA gene" in html and "1 naming a symbol two RNA" in html
        several = pl.read_csv(
            output(root, "protein_groups_several_genes_v1.csv"), comment_prefix="#"
        )
        assert several.rows() == [("G2;G3", "G2;G3")]


def test_study_rna_sample_not_gathered(tmp_path):
    """Catches a sample table naming an RNA sample gather never wrote, accepted silently."""
    root, spec = build(tmp_path, "diann")
    table = tmp_path / "fake_study_proteo.csv"
    table.write_text(table.read_text().replace("rna_03", "rna_99"))
    done = subprocess.run(
        [str(CTX), "--root", str(root), "proteo", "--dataset", SLUG, "--dataset-file", str(spec)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert done.returncode == 1 and "1 RNA samples are not in gather" in done.stderr


def test_text_frame_keeps_empty_cells():
    """Catches an empty cell read as null: DIA-NN's missing token is the empty text."""
    frame = quantities.text_frame(b"a\tb\n1\t\n", "\t")
    assert frame.rows() == [("1", "")]
    assert quantities.text_frame(b'a,b\n"x,y",2\n', ",").rows() == [("x,y", "2")]
