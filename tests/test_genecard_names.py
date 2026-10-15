"""M9: HPA's cell-line group on the card, and a symbol the genes table lacks found by HPA's name.

Each test names the mutation it catches; every one was planted and seen to fail. The two HPA
search responses are real (``tests/fixtures/genecard/README.md``).
"""

import gzip
import json

import polars as pl
import pytest
from genecard_root import FIXTURES, LINES, NPTX2, SLUG, Stub, make, responses
from test_genecard import _cards, _register, _settings

from ctx.genecard import GenecardError, cards, hpa
from ctx.genecard import run as module
from ctx.genecard.cache import Cache
from ctx.report.builder import unpacked

GENE = "ENSG00000106236"
RENAMED = 185633


def _body(name: str) -> bytes:
    return gzip.decompress((FIXTURES / name).read_bytes())


@pytest.fixture
def root(tmp_path):
    spec = make(tmp_path / "root")
    return tmp_path / "root", spec


def test_cell_lines_of_the_gene_row_in_hpa_order():
    # Mutations: another gene's row read; a misspelt group (HPA drops the column) read as no
    # lines instead of refused; the plan's measured values (93 lines, HL-60 0.0, THP-1 1.9,
    # U-937 0.4) lost.
    # Mutation: the search URL without compress=no (HPA then sends gzip) or its columns.
    assert hpa.search_url("v25", GENE, ("g", "eg", LINES[0])) == (
        "https://v25.proteinatlas.org/api/search_download.php?search=ENSG00000106236"
        "&format=json&columns=g,eg,cell_RNA_leukemia&compress=no"
    )
    body = _body("hpa_v25_search_ENSG00000106236_cell_RNA_leukemia.json.gz")
    lines = hpa.parse_lines(body, GENE, LINES[0])
    values = dict(lines.values)
    assert len(lines.values) == 93  # noqa: PLR2004 - HPA's leukemia group, plan section 1
    assert (values["HL-60"], values["THP-1"], values["U-937"]) == (0.0, 1.9, 0.4)
    assert lines.values[0][0] == "697"
    other = [{"Gene": "X", "Ensembl": "ENSG00000000001", "Cell line RNA - HL-60 [nTPM]": "99"}]
    planted = json.dumps(other + json.loads(body)).encode()
    assert dict(hpa.parse_lines(planted, GENE, LINES[0]).values)["HL-60"] == 0.0
    with pytest.raises(GenecardError, match="no cell line in group 'cell_RNA_leukaemia'"):
        hpa.parse_lines(json.dumps([{"Gene": "NPTX2", "Ensembl": GENE}]).encode(), GENE,
                        "cell_RNA_leukaemia")  # fmt: skip
    with pytest.raises(GenecardError, match="0 rows"):
        hpa.parse_lines(b"[]", GENE, LINES[0])


def test_dataset_file_cell_lines_checked():
    # Mutations: a group outside HPA's cell_RNA_ columns accepted; one line typed as a string
    # read as a list of characters.
    assert hpa.cell_lines({}, "x.toml") is None
    found = hpa.cell_lines({"hpa_cell_lines": {"group": LINES[0], "lines": list(LINES[1])}}, "x")
    assert found == hpa.CellLines(LINES[0], LINES[1])
    with pytest.raises(GenecardError, match="not a cell_RNA_ column"):
        hpa.cell_lines({"hpa_cell_lines": {"group": "leukemia", "lines": ["HL-60"]}}, "x")
    with pytest.raises(GenecardError, match="must list the lines"):
        hpa.cell_lines({"hpa_cell_lines": {"group": LINES[0], "lines": "HL-60"}}, "x")


def _renamer(root, bodies: dict[str, bytes]) -> tuple[cards.Renamer, cards.genes.GeneTable, Stub]:
    """A renamer on its own empty cache, so each served body is the one read."""
    _, table = cards.gene_table(root, SLUG)
    stub = Stub(bodies)
    folder = root / "caches" / str(len(list((root / "caches").glob("*"))))
    return cards.Renamer(Cache(folder, opener=stub), "v25", table), table, stub


def test_previous_symbol_resolves_once_and_exactly(root):
    # Mutations: HPA's fuzzy search taken as a match (a synonym that only contains the name);
    # a gene outside the genes table accepted; two genes resolved to the first; an Ensembl ID
    # sent to HPA's name search.
    url = hpa.search_url("v25", "NDUFA4L2", ("g", "gs", "eg"))
    renamer, table, stub = _renamer(root[0], {url: _body("hpa_v25_search_name_NDUFA4L2.json.gz")})
    assert cards.named("NPTX2 NDUFA4L2", table, "--genes", renamer) == [NPTX2, RENAMED]
    assert renamer.renamed == {RENAMED: "NDUFA4L2"}
    near = [
        {"Gene": "ACTB", "Gene synonym": ["NDUFA4L2-AS1"], "Ensembl": "ENSG00000075624"},
        {"Gene": "NDUFA4L2", "Gene synonym": [], "Ensembl": "ENSG00000999998"},
    ]
    renamer, table, _ = _renamer(root[0], {url: json.dumps(near).encode()})
    with pytest.raises(GenecardError, match="NDUFA4L2 names no gene"):
        cards.named("NDUFA4L2", table, "--genes", renamer)
    two = [
        {"Gene": "NDUFA4L2", "Gene synonym": [], "Ensembl": "ENSG00000075624"},
        {"Gene": "X", "Gene synonym": ["ndufa4l2"], "Ensembl": "ENSG00000185633"},
    ]
    renamer, table, _ = _renamer(root[0], {url: json.dumps(two).encode()})
    with pytest.raises(GenecardError, match="NDUFA4L2 names ENSG00000075624, ENSG00000185633"):
        cards.named("NDUFA4L2", table, "--genes", renamer)
    assert renamer.renamed == {}
    renamer, table, stub = _renamer(root[0], {})
    with pytest.raises(GenecardError, match="ENSG00000999998 names no gene"):
        cards.named("ENSG00000999998", table, "--genes", renamer)
    assert stub.calls == []


def test_offline_rename_names_why(root):
    # Mutation: an uncached name search under --offline reported as a plain unknown name.
    folder, spec = root
    options = module.Options(spec, genes=("NDUFA4L2",), offline=True)
    with pytest.raises(GenecardError, match="NDUFA4L2 names no gene; NDUFA4L2: HPA search: not"):
        module.run(_settings(folder), _cards(), "ctx genecard", SLUG, options)


def test_run_shows_lines_and_rename(root):
    # Mutations: the cell-line row or its table missing; the rename not stated on the card;
    # the renamed gene's card under the typed name rather than the table's.
    folder, spec = root
    _register(folder, spec)
    options = module.Options(spec, genes=("NPTX2", "NDUFA4L2"))
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG, options,
                         Stub(responses()))  # fmt: skip
    page = unpacked(written[0].read_text())
    assert "HL-60 0.0 nTPM; THP-1 1.9 nTPM; U-937 0.4 nTPM" in page
    assert "93 lines in the group, median 0.0 nTPM" in page
    assert "Jin et al., Nat Commun 2023" in page
    assert ("NDUFA4L2 is not in this dataset&#39;s gene table" in page
            or "NDUFA4L2 is not in this dataset's gene table" in page)  # fmt: skip
    assert "COXFA4L2, ENSG00000185633" in page
    lines = pl.read_parquet(next(p for p in written if p.name.endswith("_lines_v1.parquet")))
    assert lines.height == 93 and set(lines["group"]) == {LINES[0]}  # noqa: PLR2004
