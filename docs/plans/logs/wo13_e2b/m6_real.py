"""M6's running score, live ORA by contrast and background view on a real page, in Chromium.

python m6_real.py <root> <dataset> <downloads dir>: at 1,280 and 375 px in both themes, every
running set's curve against lite._running and its ES against the GSEA table, every contrast's
live ORA rows against its ORA table (D-049's 1.52e-10), the page box against the window.
"""

import math
import sys
from pathlib import Path

import numpy as np
import polars as pl

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402
from test_function_views_browser import ORA, RUNNING  # noqa: E402

from ctx.function import lite, running  # noqa: E402

root, dataset, downloads = Path(sys.argv[1]).absolute(), sys.argv[2], Path(sys.argv[3]).absolute()
downloads.mkdir(parents=True, exist_ok=True)
m6, m5 = root / "outputs" / "06_function", root / "outputs" / "05_de"
(page,) = sorted(m6.glob(f"*_function_{dataset}_overview_v1.html"))[-1:]
SIZES = (("1280_light", 1280, "light"), ("1280_dark", 1280, "dark"),
         ("375_light", 375, "light"), ("375_dark", 375, "dark"))  # fmt: skip


def newest(folder: Path, what: str) -> pl.DataFrame:
    """The newest table of one kind."""
    return pl.read_parquet(sorted(folder.glob(f"*_{dataset}_{what}_v1.parquet"))[-1])


tasks = [{"name": f"{what}_{name}", "url": page.as_uri(), "width": width, "height": 900,
          "scheme": scheme, "act": act}
         for what, act in (("running", RUNNING), ("ora", ORA))
         for name, width, scheme in SIZES]  # fmt: skip
seen = _browse(tasks, downloads, timeout=1800)
print(f"page {page.name} {page.stat().st_size:,} B")
for name, _, _ in SIZES:
    run = seen[f"running_{name}"]
    acted = run["acted"]
    drawn_file = next(iter(sorted(m6.glob(f"*_{dataset}_running_*_v1.csv")))).name
    contrast = drawn_file.split("_running_")[1].removesuffix("_v1.csv")
    weights = running.weights(newest(m5, f"{contrast}_results"))
    gsea = newest(m6, f"{contrast}_gsea")
    worst, es_worst, flagged = 0.0, 0.0, 0
    for drawn in acted["sets"]:
        positions = np.asarray(drawn["dots"][0], dtype=int)
        expected = lite._running(positions, weights[positions], len(weights))  # noqa: SLF001
        worst = max(worst, float(np.max(np.abs(np.asarray(drawn["line"]) - expected))))
        flagged += "gives ES" in drawn["line_words"]
    es_table = gsea["es"].to_list()
    es_page = sorted(d["es"] for d in acted["sets"])
    es_worst = max(abs(a - b) for a, b in zip(sorted(es_table), es_page, strict=True))
    print(f"running {name}: {len(acted['sets'])} sets of {gsea.height}, line max |page - lite| "
          f"{worst:.2e}, ES max |page - table| {es_worst:.2e}, flagged {flagged}, box "
          f"{acted['box']} scroll {acted['scroll']}, errors {run['errors']}")  # fmt: skip
    ora = seen[f"ora_{name}"]
    rel, rows_ok = 0.0, 0
    for contrast_, drawn in ora["acted"]["contrasts"]:
        table = newest(m6, f"{contrast_}_ora")
        want = {(r["direction"], r["gene_set"]): r for r in table.iter_rows(named=True)}
        rows = [line.split("\t") for line in drawn["tsv"].splitlines()[1:]]
        rows_ok += len(rows) == table.height
        for direction, gene_set, _, _, _, p, _, _ in rows:
            ref = want[(direction, gene_set)]["p"]
            rel = max(rel, abs(float(p) - ref) / ref if ref else 0.0)
    top = [len(y) for _, y in ora["acted"]["top"]["bars"]]
    print(f"ora {name}: {rows_ok} of {len(ora['acted']['contrasts'])} contrasts row for row, "
          f"max relative p {rel:.2e}, top bars {top}, errors {ora['errors']}")  # fmt: skip
    assert math.isfinite(rel)
