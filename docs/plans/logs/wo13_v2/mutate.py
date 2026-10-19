# ruff: noqa: E501 - page scripts and mutations are quoted whole, a line each
"""Package wo13-v2: plant each new rule's mutation, run the tests that name it, and restore.

A mutation is caught when its tests fail; the file is restored byte for byte after each, and
no bytecode is written. Chromium is named, so a browser test runs and never skips.
Run: python mutate.py <group> [name ...]; a line per mutation to stdout.
"""

import os
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
CHROMIUM = "/Applications/Chromium.app/Contents/MacOS/Chromium"
BOX, BOXES, DEPAGE, GCPAGE = ("src/ctx/viz/box.py", "src/ctx/genecard/boxes.py",
                              "src/ctx/de/page.py", "src/ctx/genecard/page.py")  # fmt: skip
T_PAIRS, T_OWN, T_DEBOX = ("tests/test_box_pairs.py", "tests/test_genecard_own.py",
                           "tests/test_de_box_order.py")  # fmt: skip
GROUPS = {
    "box": {
        "two_samples_joined": (BOX, 'once = measured.filter(pl.len().over("feature", data.pairs, "group") == 1)', "once = measured", [T_PAIRS]),
        "missing_joined": (BOX, '    measured = data.frame.filter(pl.col("value").is_not_null() & pl.col("value").is_not_nan()', '    measured = data.frame.filter(pl.lit(True)', [T_PAIRS]),
        "groups_swapped": (BOX, "    reference, comparison = data.groups\n    measured", "    comparison, reference = data.groups\n    measured", [T_PAIRS]),
        "files_from_centre": (BOX, "segments = [[(i + OFFSETS[0], a), (i + OFFSETS[1], b)]", "segments = [[(i, a), (i, b)]", [T_PAIRS]),
        "files_over_points": (BOX, "linewidths=LINE_WIDTHS[\"line\"], zorder=2))", "linewidths=LINE_WIDTHS[\"line\"], zorder=4))", [T_PAIRS]),
        "files_in_ink": (BOX, "LineCollection(segments, colors=theme.reference,", "LineCollection(segments, colors=theme.ink,", [T_PAIRS]),
        "page_no_lines": (BOX, "    if not lines:\n        return None\n    x: list", "    if True:\n        return None\n    x: list", [T_PAIRS]),
        "page_boxes_grouped": (BOX, "x=[at[f] + offset for f in rows[\"feature\"]],", "x=[at[f] for f in rows[\"feature\"]],", [T_PAIRS]),
        "wrapped_break_dropped": (BOX, 'return "<br>".join(_lines(feature))', 'return "".join(_lines(feature))', [T_PAIRS, T_OWN]),
        "wrapped_measured_whole": (BOX, "*(text_width(line, tick) for line in _lines(f)),", "text_width(f, tick),", [T_PAIRS, T_OWN]),
        "phone_kept_when_wider": (BOX, "if sum(narrow[0]) < sum(wide[0]):", "if True:", [T_OWN]),
        "m5_pairs_never_set": ("src/ctx/de/run.py", 'return (replace(data, pairs="donor") if result.contrast.paired else data), missing', "return data, missing", [T_DEBOX]),
        "m5_reading_silent": (DEPAGE, 'extra = (figure("box", "reading_paired"),) if paired else ()', "extra = ()", [T_DEBOX]),
        "f1_unpaired_joined": (BOXES, 'pair=pl.when(pl.col("contrast").is_in(joined)).then("donor"))', 'pair="donor")', [T_OWN]),
        "box_caption_silent": (BOXES, '"figures.box.caption_paired" if paired else "figures.box.caption"', '"figures.box.caption"', [T_OWN]),
        "f1_views_viewed": (GCPAGE, "from ctx.report.lookups import lookup", "from ctx.report.views import viewed as lookup", [T_OWN]),
        "export_object_as_text": ("src/ctx/report/export.py", "return caption if isinstance(caption, str) else str(caption.text())", "return caption", [T_OWN]),
    },
}  # fmt: skip
EXPORT, HEAT, QB, QF = ("src/ctx/report/templates/_export.js.j2", "src/ctx/viz/heatmap.py",
                        "src/ctx/quant/blocks.py", "src/ctx/quant/figures.py")  # fmt: skip
T_QB, T_QF = "tests/test_quant_views_browser.py", "tests/test_quant_figures.py"
GROUPS["heat"] = {
    "export_raster": (EXPORT, "        if (i === 0 && ex.cells && trace.type === \"heatmap\" && ex.cells[heat]) {", "        if (false) {", [T_QB]),
    "runs_across_levels": (EXPORT, "while (end < cols && z[r][end] === v) end++;", "while (end < cols && z[r][end] !== null) end++;", [T_QB]),
    "tracks_in_cell_colours": (HEAT, 'layers: dict[str, Any] = {"cells": [colours[a:b] for a, b in itertools.pairwise(starts)]}', 'layers: dict[str, Any] = {"cells": [colours[: sizes[0]] for _ in sizes]}', [T_QF]),
    "page_tints": (HEAT, "for rgba in self._rgba(THEMES[\"white\"])]", "for rgba in self._rgba(THEMES[\"black\"])]", [T_QF]),
    "thinned_caption_files": (HEAT, '            layers["caption"] = plain(self.caption_text())', '            layers["caption"] = self.caption', [T_QF]),
    "order_clustering_kept": (QF, "        order = grouped_order(order, sheet, by)", "        order = order", [T_QF]),
    "order_levels_sorted": (QF, "    first = {v: i for i, v in enumerate(dict.fromkeys(sheet[by].cast(pl.String)))}", "    first = {v: i for i, v in enumerate(sorted(set(sheet[by].cast(pl.String))))}", [T_QF]),
    "strip_grouping_last": (QF, "    wanted = list(dict.fromkeys([by, *STRIPS] if by and by not in STRIPS else STRIPS))", "    wanted = list(dict.fromkeys([*STRIPS, by] if by and by not in STRIPS else STRIPS))", [T_QF]),
    "grouped_caption_build": (QB, "    data = fig.correlation_heatmap(rho, sheet, text.fill(PAGE, key, **values), by=column)", "    data = fig.correlation_heatmap(rho, sheet, text.fill(PAGE, key, **values))", [T_QB]),
}  # fmt: skip
RUNPY, RUNJS, PAGERUN, ORAJS, LIVE = (
    "src/ctx/function/running.py",
    "src/ctx/report/templates/_running.js.j2",
    "src/ctx/report/templates/_page_running.js.j2",
    "src/ctx/report/templates/_ora.js.j2",
    "src/ctx/function/live.py",
)
T_RUN, T_ORA, T_M6B = ("tests/test_running_page.py", "tests/test_ora_page.py",
                       "tests/test_function_views_browser.py")  # fmt: skip
GROUPS["m6"] = {
    "keep_build_first": (RUNPY, "    for wanted in (current, build):", "    for wanted in (build, current):", [T_RUN]),
    "keep_js_build_first": (RUNJS, "var wanted = [current, build];", "var wanted = [build, current];", [T_RUN]),
    "keep_js_minus_one": (RUNJS, "    return 0;\n  },\n  /* The enrichment score", "    return -1;\n  },\n  /* The enrichment score", [T_RUN]),
    "order_shown_kept_untested": (RUNPY, "    return [*([shown] if shown in enriched else []),", "    return [shown, *([] if shown in enriched else []),", [T_RUN]),
    "page_build_ranking": (PAGERUN, "var i = Number(sets.value), s = contrasts[c].sets[i], weights = weightsOf(c);", "var i = Number(sets.value), s = contrasts[c].sets[i], weights = weightsOf(0);", [T_M6B]),
    "page_sets_not_relisted": (PAGERUN, "      if (c !== shownContrast || (e && e.target === picker)) relist(c);\n", "", [T_M6B]),
    "page_export_kept": (PAGERUN, "if (built) cfg.export = moved ? Object.assign({}, built, { caption: s.export_caption }) : built;", "", [T_M6B]),
    "page_set_reset_first": (PAGERUN, "      sets.value = String(at);\n", '      sets.value = "0";\n', [T_M6B]),
    "ora_featured_calls": (ORAJS, "(cfg.contrasts || []).forEach(function (c) { if (c[0] === name) found = c; });", "found = cfg.contrasts[0];", [T_ORA]),
    "ora_featured_background": (ORAJS, "cfg.backgrounds.forEach(function (b) { if (b[0] === name) background = b; });", "background = cfg.backgrounds[0];", [T_ORA]),
    "ora_up_down_swapped": (ORAJS, "return { up: names(found[2]), down: names(found[3]), background: tested };", "return { up: names(found[3]), down: names(found[2]), background: tested };", [T_ORA]),
}  # fmt: skip
FLAG, QCB = "src/ctx/qc/flag_map.py", "src/ctx/qc/blocks.py"
T_QV = "tests/test_qc_views.py"
GROUPS["m2"] = {
    "flag_needs_no_value": (FLAG, '    if value is None or (isinstance(value, float) and math.isnan(value)):\n        return "none"\n    if flagged:', '    if flagged:\n        return "flagged"\n    if value is None or (isinstance(value, float) and math.isnan(value)):\n        return "none"\n    if flagged:', [T_QV]),
    "unjudged_as_passed": (FLAG, '    return "passed" if judged else "unjudged"', '    return "passed"', [T_QV]),
    "flags_ignored": (FLAG, "            cell = level(row[m.column], (row[\"sample\"], m.id) in flagged, judged_here)", "            cell = level(row[m.column], False, judged_here)", [T_QV]),
    "fitted_never": (FLAG, "    if len(rendered(page, path).encode()) <= MAX_BYTES:", "    if True:", [T_QV]),
    "fitted_always": (FLAG, "    if len(rendered(page, path).encode()) <= MAX_BYTES:", "    if False:", [T_QV]),
    "fitted_drops_map": (FLAG, "            items += [replace(item, controls=None),\n", "            items += [\n", [T_QV]),
    "view_export_builds": (QCB, "    return exported(block, options, text.word(PAGE, \"figures.flag_map.view\"))", "    return viewed(block, options, text.word(PAGE, \"figures.flag_map.view\"))", [T_QV]),
}  # fmt: skip


def run(group: str, name: str) -> str:
    """Plant one mutation, run its tests, restore; the verdict as one line."""
    path, old, new, tests = GROUPS[group][name]
    target = WT / path
    body = target.read_bytes()
    text = body.decode()
    if text.count(old) != 1:
        return f"{name}\tNOT PLANTED: {text.count(old)} matches"
    target.write_text(text.replace(old, new))
    try:
        env = {**os.environ, "PATH": f"{BIN}:/usr/bin:/bin:/usr/sbin:/sbin",
               "PYTHONPATH": str(WT / "src"), "PYTHONDONTWRITEBYTECODE": "1",
               "CTX_CHROMIUM": CHROMIUM}  # fmt: skip
        done = subprocess.run([f"{BIN}/python", "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider",
                               *tests], cwd=WT, env=env, capture_output=True, text=True,
                              check=False)  # fmt: skip
    finally:
        target.write_bytes(body)
    last = (done.stdout.strip().splitlines() or [""])[-1]
    return f"{name}\t{'CAUGHT' if done.returncode else 'SURVIVED'}\t{last}"


if __name__ == "__main__":
    chosen = sys.argv[1]
    for name in sys.argv[2:] or GROUPS[chosen]:
        print(run(chosen, name), flush=True)
