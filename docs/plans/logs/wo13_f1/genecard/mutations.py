"""wo13-f1 M9: plant each mutation, run the tests that should catch it, restore the file.

Usage: wo13_f1_genecard_py.sh python docs/plans/logs/wo13_f1/genecard/mutations.py [name ...]
Each mutation is one exact replacement; the file is restored even when the run fails.
"""

import subprocess
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[5]
T = "tests/test_figures_m9.py"
EVERY = f"{T}::test_every_figure_has_a_spec_and_redraws_byte_for_byte"
SHORT = f"{T}::test_a_page_short_of_cards_redraws_the_whole_list"
DRAW, KINDS = "src/ctx/figures/draw_m9.py", "src/ctx/figures/kinds_m9.py"
PAGE = "src/ctx/genecard/page.py"
MUTATIONS = [
    ("figure tables not written", PAGE,
     "    figure_tables.write(writer, ctx.chosen, ctx.evidence, shown.chosen, ctx.family,\n"
     "                        ctx.options.hpa_release)  # fmt: skip\n", "", EVERY),
    ("evidence view of the kept cards", PAGE,
     "figure_tables.write(writer, ctx.chosen, ctx.evidence,",
     "figure_tables.write(writer, ctx.chosen, shown.evidence,", SHORT),
    ("evidence view drawn over the first card", DRAW,
     'listed = _numbered(_read(inputs, "evidence_listed"))',
     'listed = _numbered(_read(inputs, "evidence_listed")).filter(pl.col("symbol") == '
     'pl.col("symbol").first())', SHORT),
    ("page laid out before the panels", PAGE,
     "built = _module_page(ctx, again,", "built = _module_page(ctx, drawn[:kept],", EVERY),
    ("box family of the kept cards", PAGE,
     "shown.chosen, ctx.family,",
     "shown.chosen, {g: t for g, t in ctx.family.items() if g in "
     "{c.gene_id for c in shown.chosen}},", SHORT),
    ("dropped samples read as none", DRAW,
     'card.exons_dropped = int(dropped["exons_dropped"][0])', "card.exons_dropped = 0", EVERY),
    ("MANE exons left out", DRAW,
     'card.gtf, card.rest = _mane(mane, "gtf"), _mane(mane, "rest")',
     "card.gtf, card.rest = None, None", EVERY),
    ("every contrast unpaired", DRAW,
     "{c: paired.get(c, False) for c in frame", "{c: False for c in frame", EVERY),
    ("another card's rows", DRAW,
     'rows = _numbered(frame).filter(pl.col("number") == number)',
     'rows = _numbered(frame).filter(pl.col("number") != number)', EVERY),
    ("HPA release from the constant", DRAW,
     'release = _read(inputs, "genecard_options")["hpa_release"][0]',
     "release = hpa.RELEASE",
     f"{T}::test_the_immune_panel_names_the_release_the_run_named"),
    ("card pattern not checked", KINDS,
     "found = [] if CARD.fullmatch(p.gene) else [", "found = [] if p.gene else [",
     f"{T}::test_a_card_is_named_as_m9_names_it"),
    ("words check left out", KINDS,
     "if kind in OWN and not words_unique(spec):", "if kind in OWN and not True:",
     f"{T}::test_contrasts_sharing_their_words_leave_the_own_panels_without_a_spec"),
]  # fmt: skip


def run(names: list[str]) -> int:
    """Plant, test, restore; print one line per mutation and return the number survived."""
    survived = 0
    for name, path, old, new, tests in MUTATIONS:
        if names and name not in names:
            continue
        file = W / path
        original = file.read_text()
        assert original.count(old) == 1, (name, original.count(old))
        file.write_text(original.replace(old, new))
        try:
            done = subprocess.run(["pytest", "-q", "-x", "-p", "no:cacheprovider", *tests.split()],
                                  cwd=W, capture_output=True, text=True, check=False)  # fmt: skip
        finally:
            file.write_text(original)
        caught = done.returncode != 0
        survived += not caught
        lines = done.stdout.strip().splitlines()
        last = lines[-1] if lines else ""
        why = next((x.strip() for x in lines if x.startswith("E ")), "")[:160]
        print(f"{'CAUGHT ' if caught else 'SURVIVED'} {name}: {last} | {why}", flush=True)
    return survived


if __name__ == "__main__":
    sys.exit(run(sys.argv[1:]))
