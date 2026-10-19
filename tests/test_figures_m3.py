"""F.1: M3's figures, each redrawn from its spec's inputs alone in a clean folder.

Every figure with a spec writes M3's data table, PDF and SVG byte for byte. Each test names
its mutation.
"""

import polars as pl
from test_figures import copied, drawn_alike, redrawn, specs

from ctx.figures.spec import parse

PROCESS = ("03_process",)
# What the HISAT2 run of the fixture draws: no tool comparison, which is STAR's alone.
DRAWN = {"mapping", "featured_settings", "strand_assigned", "strand_share",
         *(f"{k}:{c}" for k in ("counters", "agreement")
           for c in ("strand", "overlap", "fraction", "multimap"))}  # fmt: skip


def _kinds(found) -> set[str]:
    """The kinds, a comparison's named by it."""
    names = set()
    for spec in found:
        kind = parse(spec.figure_id)[2]
        extra = getattr(spec.params, "comparison", None)
        names.add(f"{kind}:{extra}" if extra else kind)
    return names


def test_every_figure_redraws_byte_for_byte(process_root, tmp_path):
    """All 12 figures of the HISAT2 run, redrawn from M3's tables. Catches a bar's outcomes
    read back in another order, the featured sample's settings taken from the samples' table,
    the gene sums left out of the scatters, and the groups that filter the bars dropped."""
    root, spec_path, _ = process_root
    found = specs(root, PROCESS, spec_path)
    assert _kinds(found) == DRAWN
    for spec in found:
        module_csv = next(root.glob(f"outputs/03_process/{spec.stem}.csv"))
        csv = redrawn(copied(spec, root, tmp_path / spec.stem), spec)
        assert drawn_alike(csv, module_csv) == [], spec.stem


def test_the_figure_tables_hold_what_the_page_drew(process_root):
    """The outcomes table adds up as the aligner's reports did, every bar's outcomes to its
    input reads, and the options name the run's aligner and featured sample. Catches the input
    reads left out of the table, which the bars' percentages divide by."""
    root = process_root[0]
    out = root / "outputs" / "03_process"
    outcomes = pl.read_parquet(next(out.glob("*_outcomes_v1.parquet")))
    summed = outcomes.filter(pl.col("outcome") != "input").group_by("label").agg(pl.sum("reads"))
    given = outcomes.filter(pl.col("outcome") == "input").select("label", given="reads")
    joined = summed.join(given, on="label")
    assert joined.height == outcomes["label"].n_unique() > 0
    assert (joined["reads"] == joined["given"]).all()
    options = pl.read_parquet(next(out.glob("*_figure_options_v1.parquet"))).row(0, named=True)
    assert options["aligner"] == "hisat2" and options["featured"] == "se_a"
