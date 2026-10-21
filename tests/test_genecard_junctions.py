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
                junction_id=V39, junction_tpm=_tpm())  # fmt: skip
    drawn = panel(card)
    assert drawn is not None
    numbers = drawn[1]
    assert (numbers["matched"], numbers["introns"], numbers["junctions"]) == (4, 4, 4)


TPM_FIXTURE = "gtex_v10_ENSG00000106236_median_tpm_madeup.json.gz"
# The made-up V10 fixture's tissues below 1 TPM (tests/fixtures/genecard/README.md).
BELOW = {"Whole_Blood", "Cells_EBV-transformed_lymphocytes", "Kidney_Cortex", "Kidney_Medulla",
         "Liver", "Lung", "Brain_Cerebellar_Hemisphere"}  # fmt: skip


def _tpm() -> dict[str, float]:
    return junctions.parse_tpm(_body(TPM_FIXTURE), V39)


def test_floor_route_names_dataset_and_its_rows_are_checked():
    # Mutations: datasetId left to the route's default; one page of 250 rows taken as the
    # whole gene; a row of another dataset, gene or unit accepted; a tissue read twice.
    url = urllib.parse.urlsplit(junctions.tpm_url(V39))
    assert url.path == "/api/v2/expression/medianGeneExpression"
    assert urllib.parse.parse_qs(url.query) == {"gencodeId": [V39], "datasetId": ["gtex_v10"],
                                                "itemsPerPage": ["100000"]}  # fmt: skip
    tpm = _tpm()
    assert len(tpm) == 54 and tpm["Pituitary"] == 240.0  # noqa: PLR2004 - the fixture's V10
    raw = json.loads(_body(TPM_FIXTURE))
    for field, value, message in (("datasetId", "gtex_v8", "not gtex_v10"),
                                  ("gencodeId", "ENSG00000106236.6", "not gtex_v10"),
                                  ("unit", "read count", "not TPM")):  # fmt: skip
        bad = json.loads(json.dumps(raw))
        bad["data"][0][field] = value
        with pytest.raises(GenecardError, match=message):
            junctions.parse_tpm(json.dumps(bad).encode(), V39)
    raw["data"].append(raw["data"][0])
    with pytest.raises(GenecardError, match="twice"):
        junctions.parse_tpm(json.dumps(raw).encode(), V39)


def _panel(table, tpm):  # noqa: ANN001, ANN202
    from ctx.genecard.cards import Card  # noqa: PLC0415
    from ctx.genecard.figures import junctions as panel  # noqa: PLC0415

    card = Card(106236, "ENSG00000106236.6", "NPTX2", "named", junctions=table,
                junction_id=V39, junction_tpm=tpm)  # fmt: skip
    return panel(card)


def test_floor_is_v10_tpm_and_comes_before_the_reference():
    # Mutation (planted, seen to fail): the floor ignored, every tissue drawn and the
    # reference taken over all 54 (7.5, 12.75, 15, 18.5 reads), as before D-106.
    table = junctions.parse(_body("gtex_v10_ENSG00000106236_junctions.json.gz"), V39)
    drawn = _panel(table, _tpm())
    assert drawn is not None
    heat, numbers = drawn
    assert (numbers["tissues"], numbers["left"], numbers["release"]) == (47, 7, "V10")
    assert not BELOW & set(heat.column_order or ())
    assert heat.row_stats == {"chr7:98617888-98619642:+": 8.0, "chr7:98619860-98624921:+": 14.0,
                              "chr7:98625167-98627164:+": 17.0,
                              "chr7:98627345-98628401:+": 20.0}  # fmt: skip
    assert "47 tissues at 1 TPM or more in GTEx V10; 7 tissues below" in heat.caption
    # The 7 below hold every zero of the fixture, so no drawn cell lacks a ratio.
    assert heat.frame["log2_ratio"].null_count() == 0


def test_zero_steps_take_the_outer_levels_without_a_pseudocount():
    # Mutation (planted, seen to fail): the junction panel's old rule, log2 null at zero with
    # no step, so a junction read only where its median is 0 drew blank. D is below the floor.
    import polars as pl  # noqa: PLC0415

    rows = [
        ("chr7:10-20:+", "chr7", 10, 20, t, m)
        for t, m in (("A", 0.0), ("B", 4.0), ("C", 4.0), ("D", 90.0))
    ]
    rows += [("chr7:30-40:+", "chr7", 30, 40, t, m) for t, m in
             (("A", 2.0), ("B", 0.0), ("C", 0.0), ("D", 90.0))]  # fmt: skip
    table = pl.DataFrame(rows, schema=["junction", "chrom", "start", "end", "tissue", "median"],
                         orient="row")  # fmt: skip
    drawn = _panel(table, {"A": 5.0, "B": 5.0, "C": 1.0, "D": 0.99})
    assert drawn is not None
    frame = drawn[0].frame
    cells = {(j, t): (lv, r) for j, t, lv, r in
             frame.select("junction", "tissue", "level", "log2_ratio").iter_rows()}  # fmt: skip
    assert cells[("chr7:10-20:+", "A")] == ("< -2", None)
    assert cells[("chr7:30-40:+", "A")] == ("≥ 2", None)
    assert cells[("chr7:30-40:+", "B")] == (None, None)
    assert cells[("chr7:10-20:+", "B")] == ("-0.5 to 0.5", 0.0)
    assert "D" not in frame["tissue"].to_list()
    assert drawn[1]["unreferenced"] == 1


def test_fewer_than_two_v10_tissues_leave_a_line_not_a_panel():
    # Mutation: the two-tissue minimum dropped, so one tissue draws a panel of ratios of 0.
    from ctx.genecard.cards import Card  # noqa: PLC0415
    from ctx.genecard.page import _junction_panel  # noqa: PLC0415

    table = junctions.parse(_body("gtex_v10_ENSG00000106236_junctions.json.gz"), V39)
    one = dict.fromkeys(_tpm(), 0.5) | {"Pituitary": 240.0}
    assert _panel(table, one) is None
    card = Card(106236, "ENSG00000106236.6", "NPTX2", "named", junctions=table,
                junction_id=V39, junction_tpm=one)  # fmt: skip
    line = _junction_panel(card)
    assert line == ("No junction panel for NPTX2: fewer than two GTEx V10 tissues reach a "
                    "median of 1 TPM, so no junction has a median across tissues.")  # fmt: skip


def _run(folder, spec, opener=None, offline=False):  # noqa: ANN001, ANN202
    options = module.Options(spec, genes=("NPTX2",), offline=offline)
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG, options, opener)
    return written, unpacked(written[0].read_text())


def test_uncached_floor_says_so_and_never_draws_every_tissue(tmp_path):
    # Mutations: the panel drawn over every tissue when the floor is missing; a 5xx on the
    # floor cached; a 404 on it not recorded, so offline says "not cached" where online said
    # "no entry"; the floor read from another release's request.
    from test_genecard_absent import Failing  # noqa: PLC0415

    from ctx.genecard.cache import MANIFEST  # noqa: PLC0415

    folder, spec = tmp_path / "root", make(tmp_path / "root")
    _register(folder, spec)
    url = junctions.tpm_url(V39)
    key = f"gtex/gtex_v10/{GENE}/median_tpm.gz"
    written, page = _run(folder, spec, Failing(responses(), {url: 503}))
    assert (
        "GTEx V10 junctions unavailable for NPTX2: no GTEx V10 median TPM for its tissue "
        "floor, HTTP Error 503: Server Error." in page
    )
    assert "Splice junctions of NPTX2" not in page
    assert not any("_junctions_" in p.name for p in written)
    assert key not in json.loads((folder / MANIFEST).read_text())
    _, page = _run(folder, spec, offline=True)
    assert ("no GTEx V10 median TPM for its tissue floor, not cached, and --offline makes no "
            "request.") in page  # fmt: skip
    assert "Splice junctions of NPTX2" not in page
    bodies = responses()
    del bodies[url]
    _, page = _run(folder, spec, Stub(bodies))
    assert json.loads((folder / MANIFEST).read_text())[key]["status"] == 404  # noqa: PLR2004
    _, again = _run(folder, spec, offline=True)
    said = "no GTEx V10 median TPM for its tissue floor, no entry: HTTP 404 on "
    assert said in page and said in again and "Splice junctions of NPTX2" not in again


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
    spread = pl.DataFrame({"tissue": ["A", "B"], "sex": ["all"] * 2, "n": [3] * 2,
                           "median": [5.0] * 2, "q1": [5.0] * 2, "q3": [5.0] * 2})  # fmt: skip
    card = Card(106236, "ENSG00000106236.6", "NPTX2", "named", tissues=spread, exons=shares,
                exons_dropped=dropped)  # fmt: skip
    drawn = exons(card)
    assert drawn is not None
    heat, _, _ = drawn
    assert words in heat.caption
