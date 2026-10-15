"""M1's depth figure: its filters name their columns in words (work-order item 2's leftover).

The selects printed CELL_TYPE and CONDITION, the column names in capitals. Mutation caught:
the labels left off the contract, so each select falls back to its column's name.
"""

from types import SimpleNamespace

import polars as pl

from ctx.gather.harmonise import _library_sizes


def test_the_depth_filters_are_labelled_in_words():
    frame = pl.DataFrame({
        "run_accession": ["R1", "R2", "R3"], "read_count": [3e7, 1e7, 2e7],
        "batch": ["", "", ""], "cell_type": ["NDN", "Mono", "NDN"],
        "condition": ["a", "b", "a"], "read_length": [101, 101, 101],
    })  # fmt: skip
    data = _library_sizes(frame, SimpleNamespace(batch_bins={}), (3e7, 2e7))
    assert data.controls()["filter_labels"] == {"cell_type": "Cell type",
                                                "condition": "Condition"}  # fmt: skip
