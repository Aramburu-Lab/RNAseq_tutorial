"""M9's genes-by-contrasts evidence view (E.1, GS.9 tier 1): M5's rows of the whole card list.

Each test names the mutation it catches; each was planted and seen to fail.
"""

import polars as pl
from genecard_root import CALLED, SLUG, TWINS
from test_genecard import _cards, _register, _settings, root  # noqa: F401 - the fixture

from ctx.de.run import newest
from ctx.genecard import figures
from ctx.genecard import run as module
from ctx.report.builder import unpacked

LABELS = {"a": "A over control", "b": "B over control"}


def _evidence(genes: int) -> pl.DataFrame:
    """Gene 1 up in a and untested in b, gene 2 down in a and not called in b, then fillers."""
    rows = [(1, "a", 3.0, 2.5, 1e-6, "up"), (1, "b", None, None, None, None),
            (2, "a", -1.5, -1.2, 1e-4, "down"), (2, "b", 0.2, 0.1, 0.5, "ns")]  # fmt: skip
    rows += [(n, c, 0.1, 0.1, 0.9, "ns") for n in range(3, genes + 1) for c in LABELS]
    return pl.DataFrame(rows, schema={"number": pl.Int64, "contrast": pl.String,
                                      "log2FoldChange": pl.Float64,
                                      "log2FoldChange_shrunk": pl.Float64, "padj": pl.Float64,
                                      "call": pl.String}, orient="row")  # fmt: skip


def test_a_called_cell_takes_its_step_an_uncalled_one_grey_an_untested_one_nothing():
    # Mutations: an uncalled cell drawn in its fold change's step; an untested cell drawn grey;
    # a called cell drawn grey; the call read from the shrunk fold change's sign.
    heat, numbers, rest = figures.evidence_view(_evidence(3), {1: "G1", 2: "G2", 3: "G3"}, LABELS)
    level = {(g, c): v for g, c, v in heat.frame.select("gene", "contrast", "level").iter_rows()}
    assert level[("G1", "A over control")] == "≥ 2"
    assert level[("G2", "A over control")] == "-2 to -1"
    assert level[("G2", "B over control")] == "not called"
    assert level[("G1", "B over control")] is None
    assert (numbers["called"], numbers["cells"], rest) == (2, 5, [])
    assert {s.level for s in heat.swatches} == {"≥ 2", "-2 to -1", "not called"}


def test_at_most_a_hundred_genes_drawn_the_rest_returned_in_list_order():
    # Mutations: the row limit off; the rest lost; the rest reordered.
    names = {n: f"G{n}" for n in range(1, 131)}
    heat, numbers, rest = figures.evidence_view(_evidence(130), names, LABELS)
    assert heat.rows_in_order() == [f"G{n}" for n in range(1, 101)]
    assert rest == list(range(101, 131))
    assert (numbers["genes"], numbers["listed"]) == (100, 130)


def test_the_page_draws_the_view_for_two_contrasts_and_not_for_one(root):  # noqa: F811
    # Mutations: the view left off a dataset with two contrasts; drawn for one; its files not
    # written; a shared symbol giving two genes one row.
    folder, spec = root
    _register(folder, spec)
    one = module.run(_settings(folder), _cards(), "ctx genecard", SLUG,
                     module.Options(spec, genes=("NPTX2",), offline=True))  # fmt: skip
    assert "The listed genes in every contrast" not in unpacked(one[0].read_text())
    for path in one:
        path.rename(path.with_name(f"one_{path.name}"))
    results = newest(folder / "outputs" / "05_de", SLUG, "a_case_results")
    for suffix in (".parquet", ".csv"):
        source = results.with_suffix(suffix)
        source.with_name(source.name.replace("a_case", "b_case")).write_bytes(source.read_bytes())
    spec.write_text(
        spec.read_text() + '\n[[contrasts]]\nname = "b_case"\nvariable = "condition"'
        '\nreference = "control"\ntested = "case"\nn = [3, 3]\nlabel = "B case"\n'
    )
    named = ("NPTX2", *(t.split(".")[0] for t in TWINS), *[v[0] for v in CALLED.values()][:2])
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG,
                         module.Options(spec, genes=named, offline=True))  # fmt: skip
    page = unpacked(written[0].read_text())
    assert "The listed genes in every contrast" in page and "B case" in page
    view = next(p for p in written if p.name.endswith("_evidence_view_v1.csv"))
    rows = pl.read_csv(view)
    assert rows["gene"].n_unique() == 5  # noqa: PLR2004 - NPTX2, two TWINs apart, UPA, UPB
    assert any(p.name.endswith("_evidence_view_v1_white.pdf") for p in written)
