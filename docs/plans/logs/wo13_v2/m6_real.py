"""V2 M6 on a real page in Chromium: the running score's contrast select and ORA's contrast.

python m6_real.py <root> <dataset> <downloads dir>: at 1,280 and 375 px in both themes, every
contrast's running line against lite._running on that contrast's ranking, its ES against its
GSEA table, its set list against running.order and the set kept against running.keep; every
contrast's live ORA rows against live.by_contrast (D-049's 1.52e-10); the page box.
"""

import sys
from pathlib import Path

import numpy as np
import polars as pl

WT = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(WT / "tests")]

from test_export_browser import _browse  # noqa: E402
from test_function_views_browser import CONTRASTS, ORA  # noqa: E402

from ctx.dataset import load  # noqa: E402
from ctx.de.page import every  # noqa: E402
from ctx.function import lite, live, running  # noqa: E402
from ctx.function.page import set_words  # noqa: E402

root, dataset, downloads = Path(sys.argv[1]).absolute(), sys.argv[2], Path(sys.argv[3]).absolute()
downloads.mkdir(parents=True, exist_ok=True)
m6, m5 = root / "outputs" / "06_function", root / "outputs" / "05_de"
(page,) = sorted(m6.glob(f"*_function_{dataset}_overview_v1.html"))[-1:]
spec = load(WT / "conf" / "datasets" / f"{dataset}.toml")
by_label = {w.label: n for n, w in every(spec.contrasts, spec.columns).items()}
build_set = spec.function.running[1]
sets = lite.read_gmt(root / spec.resource(spec.function.gene_sets).dest)
SIZES = (("1280_light", 1280, "light"), ("1280_dark", 1280, "dark"),
         ("375_light", 375, "light"), ("375_dark", 375, "dark"))  # fmt: skip


def newest(folder: Path, what: str) -> pl.DataFrame:
    """The newest table of one kind."""
    return pl.read_parquet(sorted(folder.glob(f"*_{dataset}_{what}_v1.parquet"))[-1])


tasks = [{"name": f"{what}_{name}", "url": page.as_uri(), "width": width, "height": 900,
          "scheme": scheme, "act": act}
         for what, act in (("contrasts", CONTRASTS), ("ora", ORA))
         for name, width, scheme in SIZES]  # fmt: skip
seen = _browse(tasks, downloads, timeout=2400)
print(f"page {page.name} {page.stat().st_size:,} B")
for name, _, _ in SIZES:
    acted = seen[f"contrasts_{name}"]["acted"]
    shown, worst, es_worst, ok_lists, ok_kept, ok_words = build_set, 0.0, 0.0, 0, 0, 0
    for i, drawn in enumerate(acted["contrasts"]):
        contrast = by_label[drawn["contrast"]]
        gsea = {r["gene_set"]: r for r in newest(m6, f"{contrast}_gsea").iter_rows(named=True)}
        order = running.order(gsea, build_set)
        ok_lists += drawn["sets"] == [set_words(s) for s in order]
        kept = order[running.keep(order, shown, build_set)]
        ok_kept += drawn["set"] == set_words(kept)
        shown = kept
        results = newest(m5, f"{contrast}_results")
        weights = running.weights(results)
        positions = np.asarray(drawn["dots"], dtype=int)
        expected = lite._running(positions, weights[positions], len(weights))  # noqa: SLF001
        worst = max(worst, float(np.max(np.abs(np.asarray(drawn["line"]) - expected))))
        es_worst = max(es_worst, abs(drawn["es"] - gsea[kept]["es"]))
        ok_words += (drawn["moved"] == (i != 0) and drawn["export"].startswith(set_words(kept))
                     and drawn["caption"].startswith(f"{set_words(kept)} in "))  # fmt: skip
    n = len(acted["contrasts"])
    back = acted["back"]["line"] == acted["contrasts"][0]["line"] and not acted["back"]["moved"]
    print(f"running {name}: {n} contrasts; lists {ok_lists}/{n}, kept {ok_kept}/{n}, words "
          f"{ok_words}/{n}; line max |page - lite| {worst:.2e}; ES max |page - table| "
          f"{es_worst:.2e}; back to build {back}; box {acted['box']} scroll {acted['scroll']}; "
          f"errors {seen[f'contrasts_{name}']['errors']}")  # fmt: skip
    ora = seen[f"ora_{name}"]
    rel, rows_ok = 0.0, 0
    for contrast_, drawn in ora["acted"]["contrasts"]:
        table = live.by_contrast(newest(m5, f"{contrast_}_results"), sets)
        want = {(r["direction"], r["gene_set"]): r for r in table.iter_rows(named=True)}
        rows = [line.split("\t") for line in drawn["tsv"].splitlines()[1:]]
        rows_ok += len(rows) == table.height and all((d, g) in want for d, g, *_ in rows)
        for direction, gene_set, _, _, _, p, _, _ in rows:
            ref = want[(direction, gene_set)]["p"]
            rel = max(rel, abs(float(p) - ref) / ref if ref else 0.0)
    print(f"ora {name}: {rows_ok} of {len(ora['acted']['contrasts'])} contrasts equal "
          f"live.by_contrast row for row, max relative p {rel:.2e}, box "
          f"{ora['acted'].get('box')}, errors {ora['errors']}")  # fmt: skip
