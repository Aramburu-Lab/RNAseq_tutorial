"""Package wo13-e2a: what Step-1's judge panel shows, counted by checkpoint.judge in Python.

Samples flagged on any strip metric at each k and strata grouping; globin's flags by strata;
agreement and distance flags by peers. The page's twin equals checkpoint.judge
(tests/test_qc_judge.py), so these are the page's numbers.
    python judge_numbers.py <qc folder> <gather sheet parquet>
"""

import json
import re
import sys
from dataclasses import replace
from pathlib import Path

import polars as pl

from ctx.qc import checkpoint, judge
from ctx.qc import registry as reg
from ctx.qc.blocks import STRIPPED

out, sheet_path = Path(sys.argv[1]), Path(sys.argv[2])
values = pl.read_csv(next(out.glob("*_tier1_metrics_v1.csv")))
sheet = checkpoint.length_groups(pl.read_parquet(sheet_path), declared=True)
samples = sheet.unique("sample", keep="first", maintain_order=True)
registry = reg.load(reg.resolve(None))
strata_options, peer_options = checkpoint.groupings(("cell_type",), ("cell_type", "condition"),
                                                    sheet.columns)  # fmt: skip
stripped = [m.id for m in registry.metrics.values()
            if m.kind in STRIPPED and m.column in values.columns]  # fmt: skip
print(f"strip metrics: {len(stripped)}")
for strata in strata_options:
    row = []
    for k in judge.K_STEPS:
        _, flags = checkpoint.judge(values, samples, replace(registry, mad_multiplier=k), strata)
        mine = flags.filter(pl.col("metric").is_in(stripped))
        row.append(f"k {k:g}: {mine['sample'].n_unique()}")
    globin = checkpoint.judge(values, samples, registry, strata)[1].filter(
        pl.col("metric") == "globin_fraction").height  # fmt: skip
    print(f"strata {strata or 'none'}: any strip metric {'; '.join(row)}; globin at k 3 {globin}")
# Each peer grouping's values as the page carries them, judged at the run's k and strata.
page = next(out.glob("*_tier1_overview_v1.html")).read_text()
for found in re.finditer(r'id="ctx-controls-\d+">(.*?)</script>', page, re.S):
    cfg = json.loads(found.group(1))
    if cfg.get("kind") != "judge" or "peered" not in cfg:
        continue
    metric = registry.metrics[cfg["metric"]["id"]]
    row = []
    for (label, _), peered in zip(cfg["peers"], cfg["peered"], strict=True):
        frame = pl.DataFrame({"sample": values["sample"], metric.column: peered},
                             schema_overrides={metric.column: pl.Float64})  # fmt: skip
        _, flags = checkpoint.judge(frame, samples, registry, ("cell_type",))
        row.append(f"{label}: {flags.height} flagged, {sum(v is None for v in peered)} no peer")
    print(f"{metric.id} by peers: {'; '.join(row)}")
