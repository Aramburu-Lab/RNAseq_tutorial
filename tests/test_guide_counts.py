"""R2: the guide's Part 1c, its counts template and the GEO table writer.

No accession is typed here: series and samples are made up. Each test names the mutation it
catches.
"""

import gzip
import re
import shutil
import subprocess
from pathlib import Path

import pytest
from test_gather_counts import gtf_text, matrix, whole
from test_guide import codes, parse

from ctx.dataset import load
from ctx.gather import counts_run, geo_table
from ctx.report import guide, text

ROOT = Path(__file__).resolve().parents[1]
FOUR = ROOT / "conf" / "datasets" / "four_runs.toml"
R2 = ROOT / "conf" / "datasets" / "r2_gse130605.toml"
TEMPLATE = ROOT / guide.COUNTS_TEMPLATE
# A series matrix as GEO writes it: quoted, tab-separated, one line per field.
SERIES = """!Series_title\t"A made-up study"
!Sample_title\t"T_1"\t"T_2"\t"N_1"\t"N_2"
!Sample_geo_accession\t"X01"\t"X02"\t"X03"\t"X04"
!Sample_source_name_ch1\t"tumour tissue"\t"tumour tissue"\t"normal tissue"\t"normal tissue"
!Sample_characteristics_ch1\t"stage: 2"\t"stage: 3"\t"stage: NA"\t"stage: NA"
!Sample_characteristics_ch1\t"Tobacco Smoking: yes"\t"Tobacco Smoking: no"\t""\t""
!series_matrix_table_begin
"""


def part(blocks) -> dict[str, str]:
    """Part 1c's code blocks, by their key in the page text."""
    shown = codes(blocks)
    keys = ("find", "size", "download", "table", "run")
    return {k: shown[text.word("guide", f"part1c.{k}")].text for k in keys}


@pytest.fixture(scope="module")
def blocks():
    return part(guide.blocks(load(FOUR), FOUR))


def test_part1c_follows_part1b(blocks):
    """Part 1c sits after Part 1b and before Part 2, its blocks in order. Catches the group
    left off the page or moved after the count."""
    groups = [g.title for g in guide.blocks(load(FOUR), FOUR)]
    at = groups.index(text.fill("guide", "part1c.group", n=1))
    assert groups[at - 1] == text.fill("guide", "part1b.group", n=1)
    assert groups[at + 1] == text.fill("guide", "part2.group", n=2)


def test_part1c_ctx_lines_parse(blocks):
    """Each ctx line of Part 1c is one the CLI takes, and ctx run is given the matrix.
    Catches --counts renamed in the CLI alone, and the check run without the matrix."""
    lines = [line for line in blocks["run"].splitlines() if " ctx --root " in line]
    assert [parse(line).info_name for line in lines] == ["setup", "run", "run"]
    assert all('--counts "$RAW/$FILE"' in line for line in lines[1:])
    assert lines[1].endswith(" --check") and not lines[2].endswith(" --check")


def test_part1c_reads_the_size_before_downloading(blocks):
    """The size block downloads nothing, the download block reads no size, and the size comes
    first. Catches the two in one block, which downloads unread when pasted."""
    assert "curl -fsSIL" in blocks["size"] and " -o " not in blocks["size"]
    assert "content-length" not in blocks["download"] and '-o "$FILE"' in blocks["download"]
    order = list(blocks)
    assert order.index("size") < order.index("download")


def test_part1c_template_and_table_are_the_repositorys(blocks):
    """The table block runs the GEO writer and copies the counts template, both present.
    Catches a template path or module renamed in the code and not on the page."""
    assert f"cp {guide.COUNTS_TEMPLATE} " in blocks["table"]
    assert TEMPLATE.is_file()
    assert "python -m ctx.gather.geo_table " in blocks["table"]


def test_part1c_filled_from_a_counts_dataset():
    """The proof's dataset file, which has [counts], fills the series and slug; the worked
    example leaves placeholders. Catches the worked example's slug written into a matrix route
    it does not take."""
    spec = load(R2)
    filled = part(guide.blocks(spec, R2))["find"].splitlines()[0]
    assert f"GSE={spec.accessions[0]} SLUG={spec.slug}" in filled
    four = part(guide.blocks(load(FOUR), FOUR))["find"].splitlines()[0]
    assert f"GSE={guide.SERIES} SLUG={guide.SLUG}" in four


def test_geo_folder_rule():
    """The series folder is the accession with its last three digits as nnn, the folder the
    proof downloaded the proof's series matrix from; bash expands the page's rule to it.
    Catches a cut of another length, and a rule bash does not expand."""
    series = load(R2).accessions[0]
    rule = re.search(r"\$\{GSE%\?+\}nnn", guide.GEO_FTP)
    assert rule is not None
    done = subprocess.run(["bash", "-c", f"GSE={series}; echo {rule.group()}"],
                          capture_output=True, text=True, check=True)  # fmt: skip
    assert done.stdout.strip() == series[:-3] + "nnn"


def test_geo_table_writes_samples(tmp_path):
    """The writer's table: one row per sample, sample then donor, each characteristic a column
    under its key, a missing one empty. Catches a key kept with its spaces and capitals, the
    rows shifted by an empty characteristic, and a donor left unset."""
    series = tmp_path / "s.txt.gz"
    series.write_bytes(gzip.compress(SERIES.encode()))
    out = tmp_path / "t.csv"
    assert geo_table.main([str(series), str(out), "--donor-pattern", r"_(\d+)$"]) == 0
    lines = out.read_text().splitlines()
    assert lines[0] == "sample,donor,title,source_name,stage,tobacco_smoking"
    assert lines[1:] == ["X01,1,T_1,tumour tissue,2,yes", "X02,2,T_2,tumour tissue,3,no",
                         "X03,1,N_1,normal tissue,NA,", "X04,2,N_2,normal tissue,NA,"]  # fmt: skip
    assert geo_table.main([str(series), str(out)]) == 0
    assert out.read_text().splitlines()[1].startswith("X01,X01,")


def test_geo_table_refuses_a_pattern_without_its_donor(tmp_path, capsys):
    """A pattern that misses a title stops the writer, naming it. Catches a donor silently
    left empty, which curation would read as one individual."""
    series = tmp_path / "s.txt"
    series.write_text(SERIES.replace('"N_2"', '"N-two"'))
    assert (
        geo_table.main([str(series), str(tmp_path / "t.csv"), "--donor-pattern", r"_(\d+)$"]) == 1
    )
    assert "'N-two'" in capsys.readouterr().err
    assert not (tmp_path / "t.csv").exists()


def test_counts_template_curates_a_geo_table(tmp_path):
    """The counts template, as the guide copies it, with the writer's table and a matrix of
    its samples, passes gather's checks once its EDIT lines are filled. Catches a join or a
    column source in the template that the writer's table cannot satisfy."""
    root = tmp_path / "root"
    spec_path = tmp_path / "mine.toml"
    shutil.copy(TEMPLATE, spec_path)
    edited = spec_path.read_text().replace(
        '[columns.condition]\nsource = "source_name"',
        '[columns.condition]\nsource = "source_name"\n'
        'values = { "tumour tissue" = "tumour", "normal tissue" = "normal" }',
    )
    spec_path.write_text(edited)
    series = tmp_path / "s.txt"
    series.write_text(SERIES)
    assert geo_table.main([str(series), str(tmp_path / "mine.csv")]) == 0
    spec = load(spec_path)
    gtf = root / spec.resource("annotation_gtf").dest
    gtf.parent.mkdir(parents=True)
    gtf.write_bytes(gzip.compress(gtf_text().encode()))
    entrez = root / spec.resource("gencode_entrez").dest
    entrez.write_bytes(gzip.compress(b"ENST00000000011.1\t101\n"))
    hgnc = root / spec.resource("hgnc_complete").dest
    hgnc.parent.mkdir(parents=True)
    hgnc.write_text("hgnc_id\tentrez_id\tensembl_gene_id\nHGNC:1\t102\tENSG00000000002\n")
    path = matrix(
        tmp_path / "m.tsv", ["101", "102"], whole(2), samples=["X01", "X02", "X03", "X04"]
    )
    found = counts_run.check(root, spec, path)
    assert found.sheet["condition"].to_list() == ["tumour", "tumour", "normal", "normal"]
    assert found.mapped.genes["gene_name"].to_list() == ["ALPHA", "BETA"]
    assert found.mapped.through_hgnc == 1


def test_template_guide_keeps_the_slug_a_placeholder():
    """The guide made from the counts template names SLUG as a placeholder, not the template's
    own slug. Catches the proof's first gap: SLUG=mine, kept by a reader, named every folder
    after the template while the edited dataset file named another slug."""
    filled = part(guide.blocks(load(TEMPLATE), TEMPLATE))["find"].splitlines()[0]
    assert f"SLUG={guide.SLUG}" in filled


def test_part1c_notes_name_what_the_proof_lacked():
    """The proof's other gaps, each now on the page: an example donor pattern, where the
    annotation release is named, and the table's file name in the dataset file. Catches any
    note cut back to its first version."""
    note, where, edit = (text.word("guide", f"part1c.{k}") for k in ("table_note", "where", "edit"))
    assert "--donor-pattern '_(\\d+)$'" in note
    assert "annotation release" in where and "`annotation`" in where
    assert "`$SLUG.csv`" in edit and "same as SLUG" in edit
    paired = text.word("guide", "part1c.paired")
    assert "`paired = true`" in paired and '`design = "~donor + group"`' in paired


def test_geo_table_reads_descriptions(tmp_path):
    """GEO's sample descriptions become columns as characteristics do, so an authors' matrix
    can join on library_name; a description with no key goes to one description column.
    Catches the descriptions skipped (R2.2's featureCounts proof could not join), and a
    keyless line written as an empty column named after its text."""
    lines = ('!Sample_description\t"Library name: T1"\t"Library name: T2"\t"Library name: N1"'
             '\t"Library name: N2"\n!Sample_description\t"free text"\t""\t"x"\t"y"\n')  # fmt: skip
    series = tmp_path / "s.txt"
    series.write_text(SERIES.replace("!series_matrix_table_begin\n", lines))
    assert geo_table.main([str(series), str(tmp_path / "t.csv")]) == 0
    header, first, second = (tmp_path / "t.csv").read_text().splitlines()[:3]
    assert header.split(",")[-2:] == ["library_name", "description"]
    assert first.endswith(",T1,free text")
    assert second.endswith(",T2,")


def test_part1c_names_the_authors_join():
    """The proof's gaps on an authors' matrix, now on the page and in the template: the join
    on library_name and a donor from a subject column. Catches the note dropped."""
    note = text.word("guide", "part1c.authors")
    assert '`join = ["library_name", "sample"]`' in note and "`subject_id`" in note
    assert '["library_name", "sample"]' in TEMPLATE.read_text()
