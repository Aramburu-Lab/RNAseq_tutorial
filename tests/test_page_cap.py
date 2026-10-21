"""Y.1: the page cap is a property of the finished page (``ctx.report.cap``).

A page's optional parts leave past the cap as its page as written decides, Reproduce panels
and all: M2's flag view and M9's cards after the first. Each fixture page is built once in full,
then again with the cap lowered to one byte under it. Each test names the mutation it catches.
"""

import subprocess
import sys

import pytest
from genecard_root import SLUG
from test_figures_m9 import FOLDER, GENES, _built
from test_qc_run import ENV
from test_qc_views import page  # noqa: F401 - the M2 fixture page, with its flag view

from ctx.report import cap, text
from ctx.report.builder import unpacked

# ctx qc as test_qc_views's page ran it, with the cap lowered in the process that writes it.
LOWERED = "import sys\nfrom ctx.report import cap\ncap.MAX_BYTES = int(sys.argv[1])\n" \
          "from ctx.cli import app\napp(sys.argv[2:], standalone_mode=False)\n"  # fmt: skip


def _note(limit: int) -> str:
    """M2's note in the flag view's place, under a cap of ``limit`` bytes."""
    return text.fill("qc", "figures.flag_map.left_out", cap=limit)


def test_m2_drops_its_flag_view_under_a_cap_one_byte_below_its_page(page, tmp_path):  # noqa: F811
    """M2's fixture page with its flag view, rebuilt under a cap one byte below its size: the
    view leaves with its note and the page is written under the cap. Catches the view found by
    a stem no real band map has (``_band_map`` at the end, before Y.1), which left it on the
    page and stopped the run."""
    from test_qc_run import RNASEQ, SPEC  # noqa: PLC0415

    root, out, _, html = page
    full = html.read_text()
    assert "every flag of the run" in full and "flag view of this map is left out" not in full
    limit = html.stat().st_size - 1
    done = subprocess.run(
        [sys.executable, "-c", LOWERED, str(limit), "--root", str(root), "qc", "--dataset",
         "desynth", "--dataset-file", str(SPEC), "--rnaseq-dir", str(RNASEQ)],
        cwd=root, env=ENV, capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    (again,) = out.glob("*_overview_v1.html")
    shown = again.read_text()
    assert again.stat().st_size <= limit
    assert _note(limit) in shown and "every flag of the run" not in shown


def test_m9_decides_its_cards_on_the_page_with_its_panels(tmp_path, monkeypatch):
    """M9's three-card fixture page, rebuilt under a cap one byte below it: the cards fitted
    without their panels all fit, the page with them does not, so the last card leaves, with
    no file of its own, and the page is written under the cap. Catches the fit made before the
    panels existed, which kept every card and stopped the run at write."""
    full = tmp_path / "a" / "root"
    _built(full)
    (page_path,) = (full / "outputs" / FOLDER[0]).glob("*_overview_v1.html")
    limit = page_path.stat().st_size - 1
    monkeypatch.setattr(cap, "MAX_BYTES", limit)
    decided = []
    measured = cap.kept

    def spy(pages, parts, least=0):
        """``cap.kept``, each answer recorded."""
        decided.append(measured(pages, parts, least))
        return decided[-1]

    monkeypatch.setattr(cap, "kept", spy)
    lowered = tmp_path / "b" / "root"
    _built(lowered)
    folder = lowered / "outputs" / FOLDER[0]
    (shown,) = folder.glob("*_overview_v1.html")
    assert shown.stat().st_size <= limit
    assert decided[0] == len(GENES) and decided[-1] < len(GENES), decided
    left = GENES[decided[-1] :]
    assert f"Left out, in list order: {', '.join(left)}." in unpacked(shown.read_text())
    named = {p.name for p in folder.iterdir()}
    assert not [n for n in named if any(f"_{g.lower()}_" in n for g in left)], sorted(named)
    assert [p.name for p in folder.iterdir() if f"genecard_{SLUG}_box_" in p.name]


def test_one_cap_for_every_reader(monkeypatch, tmp_path):
    """A cap lowered in ``ctx.report.cap`` holds for ``builder.write`` and for a fit alike, one
    byte either side. Catches a module that binds its own copy of the number, as M2's and M9's
    did, which a lowered cap leaves at 15 MB."""
    import plotly.graph_objects as go  # noqa: PLC0415

    from ctx.report.builder import FigureBlock, FigurePage, write  # noqa: PLC0415

    block = FigureBlock(go.Figure(), "No test.", [], "a.csv")
    size = write(FigurePage("Cap", block), tmp_path / "a.html").stat().st_size
    monkeypatch.setattr(cap, "MAX_BYTES", size)
    assert write(FigurePage("Cap", block), tmp_path / "b.html").stat().st_size == size
    assert cap.kept(lambda k: ("x" * (size - 1 + k), ()), 1) == 1
    monkeypatch.setattr(cap, "MAX_BYTES", size - 1)
    with pytest.raises(ValueError, match=f"over the {size - 1}-byte cap"):
        write(FigurePage("Cap", block), tmp_path / "c.html")
    assert cap.kept(lambda k: ("x" * (size - 1 + k), ()), 1) == 0


def test_parts_without_sizes_leave_one_at_a_time(monkeypatch):
    """Three optional parts of 100 bytes whose sizes the page does not give, under a cap that
    holds two: one leaves, and the page is measured again before another does. Catches the
    excess left uncovered by an unmeasured part, which took every part out at once."""
    monkeypatch.setattr(cap, "MAX_BYTES", 1250)
    assert cap.kept(lambda k: ("x" * (1000 + 100 * k), ()), 3) == 2  # noqa: PLR2004
    assert cap.kept(lambda k: ("x" * (1000 + 100 * k), ()), 3, least=3) == 3  # noqa: PLR2004
