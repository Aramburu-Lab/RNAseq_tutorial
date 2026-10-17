# ruff: noqa: E501 - page scripts and mutations are quoted whole, a line each
"""Package wo13-e2b: plant each new rule's mutation, run the tests that name it, and restore.

A mutation is caught when its tests fail; the file is restored byte for byte after each, and
no bytecode is written (lessons: a same-size mutation can outlive its restore in a .pyc).
Run: python mutate.py <group> [name ...]; a line per mutation to stdout.
"""

import os
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
GRIDJS = "src/ctx/report/templates/_grid.js.j2"
CONTROLS = "src/ctx/report/templates/_page_controls.js.j2"
FIT = "src/ctx/report/templates/_page_fit.js.j2"
VIEWGRID = "src/ctx/report/templates/_view_grid.html.j2"
LOOKUPS = "src/ctx/report/lookups.py"
M5 = "src/ctx/de/figures.py"
M5RUN = "src/ctx/de/run.py"
T_LOOK, T_M5F, T_M5B = (
    "tests/test_lookups.py",
    "tests/test_de_figures.py",
    "tests/test_de_views_browser.py",
)
GROUPS = {
    "m5": {
        "py_key_separator": (LOOKUPS, 'return ",".join(str(c) for c in choices)', 'return ";".join(str(c) for c in choices)', [T_LOOK]),
        "js_key_separator": (GRIDJS, 'return String(c); }).join(",");', 'return String(c); }).join(";");', [T_LOOK]),
        "js_missing_key_not_build": (GRIDJS, "return found === undefined ? 0 : found;", "return found;", [T_LOOK]),
        "grid_start_not_build": (LOOKUPS, '"start": at[0][i]}', '"start": 0}', [T_LOOK]),
        "grid_repeat_accepted": (LOOKUPS, "if len(grid) != n or n != total:", "if n != total:", [T_LOOK]),
        "option_keeps_template": (LOOKUPS, 'if k != "template"}', "}", [T_LOOK]),
        "option_export_is_builds": (LOOKUPS, '**({"export": _exported(data, block.stem)} if i else {})', "", [T_LOOK]),
        "counted_untested_genes": (M5, 'return result.results.filter(pl.col("padj").is_not_null())', "return result.results", [T_M5F]),
        "counted_cutoff_ignored": (M5, "pl.col(setting.basis), pl.col(\"log2FoldChange\"), setting.alpha, setting.lfc_min", "pl.col(setting.basis), pl.col(\"log2FoldChange\"), ALPHA, setting.lfc_min", [T_M5F]),
        "counted_basis_ignored": (M5, "pl.col(setting.basis), pl.col(\"log2FoldChange\"), setting.alpha, setting.lfc_min", "pl.col(\"padj\"), pl.col(\"log2FoldChange\"), setting.alpha, setting.lfc_min", [T_M5F]),
        "build_not_first": (M5, "return sorted(grid, key=lambda s: s[0] != BUILD)", "return grid", [T_M5F]),
        "agreement_sides_swapped": (M5, '.when(x).then(pl.lit("x_only"))\n            .when(y).then(pl.lit("y_only"))', '.when(x).then(pl.lit("y_only"))\n            .when(y).then(pl.lit("x_only"))', [T_M5F]),
        "agreement_opposite_as_same": (M5, '.when(x & y & (call_x == call_y)).then(pl.lit("same"))', '.when(x & y).then(pl.lit("same"))', [T_M5F]),
        "deseq_label_from_mwu": (M5, 'return figure("box", f"deseq_{key}", p=f"{p:.2g}", lfc=lfc)', "return Caption.label(self, feature, key)", [T_M5B]),
        "deseq_padj_wrong_gene": (M5, "padj = {s: found.get(g, (None, None))[0] for s, g in symbols}", "padj = {s: found.get(symbols[0][1], (None, None))[0] for s, g in symbols}", [T_M5B]),
        "grid_lead_last_on_tie": ("src/ctx/de/page.py", "most = max(labels, key=lambda n: counts[n])", "most = max(reversed(list(labels)), key=lambda n: counts[n])", [T_M5B]),
        "view_export_builds": (CONTROLS, "if (built.export) cfg.export = (i && o.export) || built.export;", "", [T_M5B]),
        "view_no_refit": (CONTROLS, '      div.dispatchEvent(new Event("ctx-view"));\n', "", [T_M5B]),
        "refit_keeps_build_phone": (FIT, "      narrow = m.narrow;\n", "", [T_M5B]),
        "grid_select_start_first": (VIEWGRID, '{% if loop.index0 == axis.start %} selected{% endif %}', '{% if loop.first %} selected{% endif %}', [T_M5B]),
        "grid_sync_on_any_change": (GRIDJS, 'if (e.target.hasAttribute && e.target.hasAttribute("data-axis")) sync();', "sync();", ["tests/test_view_browser.py"]),
        "box_no_export": ("src/ctx/viz/box.py", '    def export_layers(self) -> dict[str, Any]:\n        """Nothing: the page draws every box, point and label the figure files draw."""\n        return {}\n', "", [T_M5B]),
    },
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
        env = {**os.environ, "PATH": f"{BIN}:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin",
               "PYTHONPATH": str(WT / "src"), "PYTHONDONTWRITEBYTECODE": "1"}  # fmt: skip
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
