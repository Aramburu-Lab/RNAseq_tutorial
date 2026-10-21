"""E.1a: the M5 evidence M9 carries for a list from a link, and where its cards sit.

Each test names the mutation it catches; each was planted and seen to fail.
"""

import json
import re

import polars as pl
import pytest
from genecard_root import SLUG, Stub, make, responses
from test_genecard import _cards, _register, _settings

from ctx.genecard import figures, listing
from ctx.genecard import run as module
from ctx.report import genes
from ctx.report.builder import unpacked
from ctx.viz.heatmap import diverging_levels

TABLE = genes.table(["ENSG00000000001.1", "ENSG00000000002.1", "ENSG00000000003.1"],
                    ["AAA", "", "CCC"])  # fmt: skip


def _results(rows: list[tuple[int, float | None, float | None, str | None]]) -> pl.DataFrame:
    schema = {"number": pl.Int64, "log2FoldChange": pl.Float64, "padj": pl.Float64,
              "call": pl.String}  # fmt: skip
    return pl.DataFrame(rows, schema=schema, orient="row")


def test_every_tested_gene_travels_rounded_with_m5s_own_call():
    """Gene 1 is tested in the first contrast only, 2 in both, 3 in neither. Mutations caught:
    a gene no contrast tested carried; the call re-derived from the rounded padj (0.04996
    rounds to 0.05 and would lose its call); the digits not cut; an untested cell given a
    call."""
    first = _results([(1, 1.23456, 0.04996, "up"), (2, -0.5, 0.5, "ns"), (3, 2.0, None, None)])
    second = _results([(1, 0.7, None, None), (2, -3.21987, 1.234e-7, "down"), (3, 1.0, None, None)])
    data = listing.evidence({"a": first, "b": second}, TABLE)
    assert data["numbers"] == [1, 2] and data["symbols"] == ["AAA", ""]
    assert data["lfc"] == [[1.23, -0.5], [0.7, -3.22]]
    assert data["padj"] == [[0.05, 0.5], [None, 1.2e-07]]
    assert data["calls"] == ["un", "-d"]
    assert data["steps"] == ["5-", "-0"]


def test_a_called_cell_travels_with_the_evidence_views_step_of_its_exact_fold_change():
    """E.1d: 1.99996 rounds to 2.00, a step higher than its own, and -0.999 is called down at
    a shrunk fold change the step does not read. Mutations caught: the step taken from the
    rounded log2 FC; an uncalled cell given a step; the step of a down call read as up; the
    page's heatmap keyed on the view's legend in another order."""
    first = _results([(1, 1.99996, 1e-5, "up"), (2, -0.999, 1e-6, "down"), (3, 2.5, 0.5, "ns")])
    data = listing.evidence({"a": first}, genes.table(
        [f"ENSG{n:011d}.1" for n in (1, 2, 3)], ["A", "B", "C"]))  # fmt: skip
    steps = [s.level for s in diverging_levels([])[1]]
    assert data["lfc"][0][0] == 2.0  # noqa: PLR2004 - the rounded value crosses the edge
    assert [steps[int(c)] if c != "-" else None for c in data["steps"][0]] == [
        "1 to 2", "-1 to -0.5", None]  # fmt: skip
    scale = figures.live_scale(["A over control"])
    assert scale["levels"][:-1] == steps and scale["levels"][-1] == "not called"
    assert [name for name, _ in scale["legend"]] == scale["levels"]
    assert len(scale["ticks"]) == 1 and scale["least_px"] > scale["label_px"] > 0


def test_the_block_reads_back_and_holds_no_timestamp():
    """Mutations caught: gzip's mtime left in, so two builds differ; NaN written as a number."""
    data = {"numbers": [1], "symbols": ["A"], "lfc": [[None]], "padj": [[None]], "calls": ["-"]}
    assert listing.unblock(listing.block(data)) == data
    assert listing.block(data) == listing.block(data)
    assert listing.rounded(float("nan"), 3) is None and listing.rounded(float("inf"), 2) is None


@pytest.fixture
def page(tmp_path):
    folder = tmp_path / "root"
    spec = make(folder)
    _register(folder, spec)
    options = module.Options(spec, genes=("NPTX2", "ACTB"))
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG, options,
                         Stub(responses()))  # fmt: skip
    return unpacked(written[0].read_text())


def test_the_page_carries_the_tested_genes_and_links_each_card_to_its_group(page):
    """The panel's block holds every gene the root's M5 table tested; each card's anchor is
    the id of the group that holds its card. Mutations caught: the anchors counted without the
    list's own group; the command carrying a gene list of its own."""
    settings = json.loads(re.search(r'id="ctx-panel-1">([^<]*)</script>', page)[1])
    data = listing.unblock(settings["block"])
    assert data["numbers"] == sorted([75624, 900011, 900012, 900013, 900014])
    groups = dict(re.findall(r'<div class="group" id="(group-\d+)">\s*<h3>([^<]*)</h3>', page))
    for number, anchor, symbol in settings["cards"]:
        assert symbol in groups[anchor] and f"ENSG{number:011d}" in groups[anchor]
    assert " genecard --dataset gcsynth " in settings["command"]
    assert "--genes" not in settings["command"]
    assert groups["group-1"] == "Genes from a link"
