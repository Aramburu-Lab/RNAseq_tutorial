"""ED.5b: a same-day rerun keeps the earlier output it would overwrite (ctx.io.keep).

Each write below that stands for a new run first clears what this process wrote
(``keep._OURS``), as a new process starts empty. Mutations caught: the old file overwritten in
place (no superseded copy); an unchanged rerun moving its files aside anyway; the staged copy
left behind; one run's own rewrite moved aside; a failed write leaving a staged file.
"""

import polars as pl
import pytest

from ctx.io import keep
from ctx.io.readme_companion import write_table


def run(path, text):
    """One run writing ``text`` at ``path``."""
    keep._OURS.clear()
    keep.written(path, lambda p: p.write_text(text))


def test_a_changed_file_moves_aside_and_an_unchanged_one_stays(tmp_path):
    page = tmp_path / "2026-10-23_qc_demo_overview_v1.html"
    run(page, "first")
    run(page, "first")
    assert not list(tmp_path.glob("superseded_*"))
    run(page, "second")
    (aside,) = tmp_path.glob(f"superseded_{keep.RUN_STAMP}")
    assert (aside / page.name).read_text() == "first" and page.read_text() == "second"
    run(page, "third")
    assert {f.name for f in aside.iterdir()} == {page.name, f"{page.stem}.2{page.suffix}"}
    assert not (tmp_path / keep.STAGING).exists()


def test_one_run_rewriting_its_own_file_replaces_it(tmp_path):
    page = tmp_path / "2026-10-23_genecard_demo_overview_v1.html"
    run(page, "draft")
    keep.written(page, lambda p: p.write_text("final"))
    assert page.read_text() == "final" and not list(tmp_path.glob("superseded_*"))


def test_a_failed_write_leaves_nothing_new(tmp_path):
    page = tmp_path / "2026-10-23_qc_demo_overview_v1.html"
    run(page, "kept")

    def fail(p):
        p.write_text("half")
        raise OSError("disk full")

    keep._OURS.clear()
    with pytest.raises(OSError, match="disk full"):
        keep.written(page, fail)
    assert page.read_text() == "kept" and sorted(p.name for p in tmp_path.iterdir()) == [page.name]


def test_a_table_rewritten_by_its_writer_keeps_the_earlier_rows(tmp_path):
    columns = {"gene": ("none", "Gene."), "n": ("count", "Reads.")}
    stem = "2026-10-23_qc_demo_counts_v1"
    for n in (1, 2):
        keep._OURS.clear()
        write_table(pl.DataFrame({"gene": ["A"], "n": [n]}), tmp_path, stem, columns=columns,
                    data_releases={})  # fmt: skip
    (aside,) = tmp_path.glob("superseded_*")
    assert pl.read_csv(aside / f"{stem}.csv")["n"].to_list() == [1]
    assert pl.read_csv(tmp_path / f"{stem}.csv")["n"].to_list() == [2]
