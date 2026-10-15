"""M9: GTEx V10 junctions through API v2, and the V11 re-release line in the sources footer.

Each test names the mutation it catches; every one was planted and seen to fail. NPTX2's two
GTEx responses are real (``tests/fixtures/genecard/README.md``).
"""

import gzip
import json
import urllib.parse

import pytest
from genecard_root import FIXTURES, SLUG, Stub, make, responses
from test_genecard import _cards, _register, _settings

from ctx.genecard import GenecardError, annotation, junctions
from ctx.genecard import run as module
from ctx.report.builder import unpacked

GENE = "ENSG00000106236"
V39 = "ENSG00000106236.4"


def _body(name: str) -> bytes:
    return gzip.decompress((FIXTURES / name).read_bytes())


def test_routes_name_dataset_and_gencode():
    # Mutations: datasetId left to the route's default; reference/gene left on its default
    # GENCODE v26; one page of 250 rows taken as the whole gene.
    gene = urllib.parse.parse_qs(urllib.parse.urlsplit(junctions.gene_url(GENE)).query)
    assert gene == {"geneId": [GENE], "gencodeVersion": ["v39"], "genomeBuild": ["GRCh38/hg38"]}
    found = urllib.parse.parse_qs(urllib.parse.urlsplit(junctions.junction_url(V39)).query)
    assert found == {"gencodeId": [V39], "datasetId": ["gtex_v10"], "itemsPerPage": ["100000"]}


def test_versioned_id_and_junctions_parsed():
    # Mutations: our GENCODE 50 version sent where V10 keys by v39's; a row of another
    # dataset or gene accepted; a second page left unread; junction coordinates one base off
    # the MANE introns (all 4 of NPTX2's match, the plan's check).
    assert junctions.versioned(_body("gtex_v10_ENSG00000106236_reference_v39.json.gz"),
                               GENE) == V39  # fmt: skip
    body = _body("gtex_v10_ENSG00000106236_junctions.json.gz")
    table = junctions.parse(body, V39)
    assert table.height == 216 and table["tissue"].n_unique() == 54  # noqa: PLR2004 - V10
    lookup = annotation.rest_mane(_body("ensembl_116_ENSG00000106236_lookup.json.gz"))
    assert lookup is not None
    pairs = set(zip(table["start"].to_list(), table["end"].to_list(), strict=True))
    assert pairs == set(lookup.introns())
    raw = json.loads(body)
    raw["data"][0]["datasetId"] = "gtex_v8"
    with pytest.raises(GenecardError, match="not gtex_v10"):
        junctions.parse(json.dumps(raw).encode(), V39)
    raw = json.loads(body)
    raw["paging_info"]["numberOfPages"] = 2
    with pytest.raises(GenecardError, match="2 pages"):
        junctions.parse(json.dumps(raw).encode(), V39)
    wrong = json.loads(_body("gtex_v10_ENSG00000106236_reference_v39.json.gz"))
    wrong["data"][0]["gencodeVersion"] = "v26"
    with pytest.raises(GenecardError, match="'v26', not v39"):
        junctions.versioned(json.dumps(wrong).encode(), GENE)


def test_mane_introns_from_rest_where_the_gtf_has_none():
    # Mutation: introns read from the GTF alone, so a gene GDC's v36 GTF tags no MANE Select
    # for (NDUFA4L2) shows "0 of 0" although Ensembl REST 116 has one.
    from ctx.genecard.cards import Card  # noqa: PLC0415
    from ctx.genecard.figures import junctions as panel  # noqa: PLC0415

    rest = annotation.rest_mane(_body("ensembl_116_ENSG00000106236_lookup.json.gz"))
    table = junctions.parse(_body("gtex_v10_ENSG00000106236_junctions.json.gz"), V39)
    card = Card(106236, "ENSG00000106236.6", "NPTX2", "named", rest=rest, junctions=table,
                junction_id=V39)  # fmt: skip
    _, numbers = panel(card)
    assert (numbers["matched"], numbers["introns"], numbers["junctions"]) == (4, 4, 4)


def test_page_labels_v10_and_states_the_rerelease(tmp_path):
    # Mutations: the junction panel labelled with V11; the V11 re-release line dropped from
    # the sources footer; the V10 API row missing from the sources.
    spec = make(tmp_path / "root")
    folder = tmp_path / "root"
    _register(folder, spec)
    options = module.Options(spec, genes=("NPTX2",))
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG, options,
                         Stub(responses()))  # fmt: skip
    page = unpacked(written[0].read_text())
    assert "Splice junctions of NPTX2 across tissues, GTEx V10" in page
    assert "4 of 4 MANE Select introns match a GTEx V10 junction" in page
    assert "re-released V11" in page and "2026-05-19" in page
    assert "adult-gtex/o, prefix bulk-gex/v11/rna-seq/" in page
    assert "V10 (gtex_v10), GENCODE v39" in page
    assert any(p.name.endswith("_junctions_v1.parquet") for p in written)


@pytest.mark.parametrize(
    ("dropped", "words"),
    [(1, "1 sample without exon reads left out"), (2, "2 samples without exon reads left out")],
)
def test_exon_caption_counts_samples_in_the_right_number(dropped, words):
    # C3.4: "1 samples without exon reads". Mutation: one template for every count.
    import polars as pl  # noqa: PLC0415

    from ctx.genecard.cards import Card  # noqa: PLC0415
    from ctx.genecard.figures import exons  # noqa: PLC0415

    shares = pl.DataFrame({"exon": ["E1", "E2"] * 2, "tissue": ["A", "A", "B", "B"],
                           "n": [3] * 4, "share": [0.4, 0.6, 0.5, 0.5]})  # fmt: skip
    card = Card(106236, "ENSG00000106236.6", "NPTX2", "named", exons=shares,
                exons_dropped=dropped)  # fmt: skip
    heat, _, _ = exons(card)
    assert words in heat.caption
