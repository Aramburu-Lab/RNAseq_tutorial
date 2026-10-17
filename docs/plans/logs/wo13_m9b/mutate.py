# ruff: noqa: E501 - mutations are quoted whole, a line each
"""Package M9b: plant each new rule's mutation, run the tests that name it, and restore.

A mutation is caught when its tests fail; the file is restored byte for byte after each, and
no bytecode is written (lessons: a same-size mutation can outlive its restore in a .pyc).
Run: python mutate.py [name ...]   (all when none is named); a line per mutation to stdout.
"""

import os
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
SAMPLES, BOXES = "src/ctx/genecard/samples.py", "src/ctx/genecard/boxes.py"
RUN, PAGE, LISTING = (
    "src/ctx/genecard/run.py",
    "src/ctx/genecard/page.py",
    "src/ctx/genecard/listing.py",
)
FIGURES, JS = "src/ctx/genecard/figures.py", "src/ctx/report/templates/_page_m9.js.j2"
OWN, LIST_T, M9B = (
    "tests/test_genecard_own.py",
    "tests/test_genecard_listing.py",
    "tests/test_genecard_browser.py",
)
MUTATIONS = {
    "own_factor_fitted": (SAMPLES, "levels = log2_normalised(counts, one, gene_ids)",
                          'levels = log2_normalised(counts, one.with_columns(size_factor=pl.lit(1.0)), gene_ids)', [OWN]),
    "batch_not_joined": (SAMPLES, 'frame = pl.concat(frames).join(batch, on="sample", how="left", maintain_order="left")',
                         'frame = pl.concat(frames).with_columns(batch=pl.lit(""))', [OWN]),
    "min_per_group_one": (SAMPLES, "MIN_PER_GROUP = 2", "MIN_PER_GROUP = 1", [OWN]),
    "tied_group_tested": (SAMPLES, "np.ptp(np.concatenate([ref, tst])) == 0", "False", [OWN]),
    "pairs_by_row_order": (SAMPLES, '.join(rows.filter(pl.col("group") == "tested").select("donor", tst="value"),\n                      on="donor"))',
                           '.hstack(rows.filter(pl.col("group") == "tested").select(tst="value")))', [OWN]),
    "r_sign_flipped": (SAMPLES, "2 * float(u) / (ref.size * tst.size) - 1", "1 - 2 * float(u) / (ref.size * tst.size)", [OWN]),
    "zeros_ranked": (SAMPLES, "d = differences[differences != 0]", "d = differences", [OWN]),
    "bh_none": (SAMPLES, 'p_adj = false_discovery_control([t.p for t in tests.values()], method="bh")',
                "p_adj = [t.p for t in tests.values()]", [OWN]),
    "batches_uncapped": (SAMPLES, "return found if 2 <= len(found) <= most else []", "return found if 2 <= len(found) else []", [OWN]),
    "family_with_untested": (RUN, "[c.gene_id for c in chosen if c.number in tested]", "[c.gene_id for c in chosen]", [OWN]),
    "missing_table_silent": (RUN, '            card.notes["own"] = own\n', "            pass\n", [OWN]),
    "per_figure_seven": (BOXES, "PER_FIGURE = 6", "PER_FIGURE = 7", [OWN]),
    "labels_wider": (BOXES, "LABEL_CHARS = 14", "LABEL_CHARS = 30", [OWN]),
    "f1_bh_per_page": (BOXES, "tests = {c: family.tests[c] for c in drawn}",
                       "tests = samples.adjusted({c: family.tests[c] for c in drawn})", [OWN]),
    "f1_n_dropped": (BOXES, "with_n=True,", "with_n=False,", [OWN]),
    "f1_no_headroom": (BOXES, "fig.update_yaxes(range=[low, high + span * HEADROOM / 2])", "pass", [OWN]),
    "f1_ticks_to_plotly": (BOXES, 'ticktext=[n.replace("\\n", "<br>") for n in names]', "ticktext=names", [OWN]),
    "f1_phone_layout_kept": (BOXES, 'fig.layout.meta = {"least_width": upright}', "pass", [OWN]),
    "f1_one_contrast_view": (BOXES, "if len(tests) < 2:  # noqa: PLR2004 - one contrast is the box",
                             "if len(tests) < 1:", [OWN]),
    "alone_without_controls": (PAGE, "        if controls is not None:\n            # Opened alone", "        if False:\n            # Opened alone", [OWN, M9B]),
    "steps_from_rounded": (LISTING, "levels, _ = diverging_levels([fc for _, fc, _, _ in rows])",
                           "levels, _ = diverging_levels([rounded(fc, LFC_DIGITS) for _, fc, _, _ in rows])", [LIST_T]),
    "uncalled_given_step": (LISTING, "if p is not None and call in CALLED and level is not None:",
                            "if p is not None and level is not None:", [LIST_T]),
    "heat_cap_two": (LISTING, "HEAT_ROWS = 100", "HEAT_ROWS = 2", [M9B]),
    "heat_untested_grey": (JS, 'if (call === "-") { row.push(null); lines.push(""); return; }',
                           'if (call === "-") { row.push(grey); lines.push(""); return; }', [M9B]),
    "heat_key_every_step": (JS, "        if (!used[i]) return;\n", "", [M9B]),
    "heat_rows_resorted": (JS, "genes = shown.genes.slice(0, h.rows)",
                           "genes = shown.genes.slice().sort(function (a, b) { return b - a; }).slice(0, h.rows)", [M9B]),
    "heat_no_least_width": (JS, '      heatDiv.style.minWidth = Math.ceil(h.least_px) + "px";\n', "", [M9B]),
    "live_scale_legend_reordered": (FIGURES, '"legend": [[t.name, t.marker.color] for t in fig.data[1:] if t.name],',
                                    '"legend": [[t.name, t.marker.color] for t in fig.data[:0:-1] if t.name],', [LIST_T]),
}  # fmt: skip


def run(name: str) -> str:
    """Plant one mutation, run its tests, restore; the verdict as one line."""
    path, old, new, tests = MUTATIONS[name]
    target = WT / path
    body = target.read_bytes()
    text = body.decode()
    if text.count(old) != 1:
        return f"{name}\tNOT PLANTED: {text.count(old)} matches"
    target.write_text(text.replace(old, new))
    try:
        env = {**os.environ, "PATH": f"{BIN}:{os.environ['PATH']}", "PYTHONPATH": str(WT / "src"),
               "PYTHONDONTWRITEBYTECODE": "1"}  # fmt: skip
        done = subprocess.run([f"{BIN}/python", "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider",
                               *tests], cwd=WT, env=env, capture_output=True, text=True,
                              check=False)  # fmt: skip
    finally:
        target.write_bytes(body)
    last = (done.stdout.strip().splitlines() or [""])[-1]
    return f"{name}\t{'CAUGHT' if done.returncode else 'SURVIVED'}\t{last}"


if __name__ == "__main__":
    for name in sys.argv[1:] or MUTATIONS:
        print(run(name), flush=True)
