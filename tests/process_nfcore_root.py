"""A ctx root for M3 on nf-core's alignment: gather outputs named after the fixture's samples.

``tests/fixtures/rnaseq_star`` is nf-core/rnaseq 3.27.0's test profile aligned with STAR; its
five yeast samples (nf-core's own test IDs, the D20 rule) are given gather outputs from
``synth``'s generator, two cell types and two conditions, so the bars have groups to filter by.
"""

from datetime import date
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
STAR = ROOT / "tests" / "fixtures" / "rnaseq_star"
SPEC = ROOT / "tests" / "fixtures" / "datasets" / "de_synthetic.toml"
DATASET = "desynth"
SAMPLES = ("RAP1_IAA_30M_REP1", "RAP1_UNINDUCED_REP1", "RAP1_UNINDUCED_REP2", "WT_REP1",
           "WT_REP2")  # fmt: skip
DAY = date(2026, 10, 23)


def plant(root: Path) -> Path:
    """Gather's sheet, counts and genes for the fixture's samples, under ``root``; the root."""
    from synth import Synthetic, design_counts, write_gather_outputs  # noqa: PLC0415

    synthetic = design_counts((("control", "short", 3, "type_a"), ("case", "short", 2, "type_b")))
    names = dict(zip(synthetic.samples["sample"], SAMPLES, strict=True))
    renamed = Synthetic(synthetic.counts.rename(names),
                        synthetic.samples.with_columns(pl.col("sample").replace_strict(names)),
                        synthetic.truth)  # fmt: skip
    write_gather_outputs(root, DATASET, renamed, DAY)
    return root


__all__ = ["DATASET", "SAMPLES", "SPEC", "STAR", "plant"]
