"""M3's scatter views (WO8) offer a sample axis only up to ``SAMPLES_MAX`` samples, each view
carrying its genes again. The fixture page has fewer, so the cap is held here on a made-up
comparison. Mutation caught: the sample axis offered past its cap.
"""

import polars as pl

from ctx import concordance
from ctx.process import views


def _part(samples: int) -> pl.DataFrame:
    rows = [(f"s{i}", f"g{j}", j + 1, j + 2) for i in range(samples) for j in range(4)]
    frame = pl.DataFrame(rows, schema=["sample", "gene_id", "a", "b"], orient="row")
    return frame.with_columns(type_group=pl.lit("protein_coding"))


def test_the_sample_axis_stops_at_its_cap():
    for n, units in ((views.SAMPLES_MAX, views.SAMPLES_MAX + 1), (views.SAMPLES_MAX + 1, 1)):
        found, axes, at = views.subsets(_part(n))
        assert len(axes[1][1]) == units
        assert {s.sample for s in found} >= {concordance.ALL}
        assert len(found) == len(at) == len(axes[0][1]) * units
