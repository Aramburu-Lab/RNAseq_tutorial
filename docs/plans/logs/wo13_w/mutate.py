# ruff: noqa: E501 - page scripts and mutations are quoted whole, a line each
"""Package W: plant each new rule's mutation, run the tests that name it, and restore.

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
M9JS = "src/ctx/report/templates/_m9.js.j2"
PAGE = "src/ctx/report/templates/_page_m9.js.j2"
CALLS = "src/ctx/report/templates/_calls.js.j2"
CONTROLS = "src/ctx/report/templates/_page_controls.js.j2"
TWINS, M9B = "tests/test_genecard_twins.py", "tests/test_genecard_browser.py"
LISTING_T, DEBOUNCE = "tests/test_genecard_listing.py", "tests/test_calls_debounce.py"
MUTATIONS = {
    "uneven_last_of_ties": (M9JS, "(best === null || v > best.cv)", "(best === null || v >= best.cv)", [TWINS]),
    "f3_trace_left_visible": (PAGE, "{ visible: cfg.traces.map(function (l) { return k.indexOf(l) >= 0; }) }",
                              "{ visible: cfg.traces.map(function (l) { return true; }) }", [M9B]),
    "f3_axis_keeps_cells": (PAGE, 'Plotly.relayout(div, { "yaxis.categoryarray": cats });', "", [M9B]),
    "f3_scale_ignored": (PAGE, 'log = !scale || scale.value === "log";', "log = true;", [M9B]),
    "f4_reference_over_every_tissue": (PAGE, "return [r.exon, r.tissue, r.share]; }), k);",
                                       "return [r.exon, r.tissue, r.share]; }), cfg.options);", [M9B]),
    "f4_cv_axis_as_built": (PAGE, "if (div.layout.yaxis2) {", "if (false) {", [M9B]),
    "list_anchor_without_list_group": ("src/ctx/genecard/page.py", "_list_panel(ctx, drawn, len(results) + 2)",
                                       "_list_panel(ctx, drawn, len(results) + 1)", [LISTING_T]),
    "list_carries_untested": ("src/ctx/genecard/listing.py", 'frame.filter(pl.col("padj").is_not_null())["number"]',
                              'frame["number"]', [LISTING_T]),
    "list_digits_uncut": ("src/ctx/genecard/listing.py", "lfc[i], padj[i] = rounded(fc, LFC_DIGITS), rounded(p, PADJ_DIGITS)",
                          "lfc[i], padj[i] = fc, p", [LISTING_T]),
    "list_call_from_rounded_padj": ("src/ctx/genecard/listing.py", "calls[i] = CALL_CODES[call]",
                                    'calls[i] = CALL_CODES[call] if (padj[i] or 1) < 0.05 else "n"', [LISTING_T]),
    "entry_drops_fragment": ("src/ctx/report/templates/entry.html.j2", 'a.getAttribute("href") + location.hash',
                             'a.getAttribute("href")', [M9B]),
    "list_shows_untested_uncarded": (PAGE, "if (at[n] !== undefined || cards[n]) rows.push(n);", "rows.push(n);", [M9B]),
    "list_drops_carded_untested": (PAGE, "if (at[n] !== undefined || cards[n]) rows.push(n);",
                                   "if (at[n] !== undefined) rows.push(n);", [M9B]),
    "list_command_names_carded": (PAGE, "var other = found.filter(function (n) { return !cards[n]; });",
                                  "var other = found;", [M9B]),
    "debounce_none": (CALLS, "pending = t.set(function () { pending = null; fn(); }, ms);", "fn();", [DEBOUNCE]),
    "debounce_pending_kept": (CALLS, "      if (pending !== null) t.clear(pending);\n      pending = t.set(",
                              "      pending = t.set(", [DEBOUNCE]),
    "debounce_cancel_noop": (CALLS, "      if (pending !== null) t.clear(pending);\n      pending = null;\n",
                             "", [DEBOUNCE]),
    "input_redraws_per_key": (CONTROLS, "if (e.target !== unit) typed(); });", "if (e.target !== unit) update(false); });",
                              ["tests/test_calls_burst_browser.py"]),
    "driver_killed_alone": ("tests/test_export_browser.py", "os.killpg(proc.pid, signal.SIGKILL)", "proc.kill()",
                            ["tests/test_browser_driver.py"]),
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
