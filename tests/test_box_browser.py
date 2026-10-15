"""A box figure at a phone's width in Chromium: no sideways scroll, no label over another (C3.4).

The kidney page's stage box scrolled inside its card at 375 px (2026-10-14): its upright
labels needed 480 px. Skipped where no Chromium is found; set ``CTX_CHROMIUM``.
"""

from datetime import date

import pytest
from test_box_width import kidney_like
from test_export_browser import CHROMIUM, _browse

from ctx.io.naming import stem
from ctx.report.builder import (
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    write,
)
from ctx.viz.box import box_data
from ctx.viz.save import figure_block, write_figure

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
# Every tick label's box and every shown annotation's, and the plot's box against its content.
MEASURE = """new Promise(r => setTimeout(r, 1500)).then(() => {
  const box = document.querySelector('.fig-plot');
  const rect = e => {
    const b = e.getBoundingClientRect(); return [b.left, b.top, b.right, b.bottom];
  };
  return {client: box.clientWidth, scroll: box.scrollWidth,
          ticks: Array.from(box.querySelectorAll('.xtick text')).map(rect),
          notes: Array.from(box.querySelectorAll('.annotation'))
            .filter(a => a.getBoundingClientRect().width > 0).length};
})"""


def _apart(boxes: list[list[float]]) -> bool:
    """No two boxes overlap."""
    return all(a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1]
               for i, a in enumerate(boxes) for b in boxes[i + 1 :])  # fmt: skip


@pytest.fixture(scope="module")
def ran(tmp_path_factory):
    """The module page holding the kidney-like box, and the box opened alone, at 375 and 1280."""
    out = tmp_path_factory.mktemp("box")
    columns = {
        "value": ("log2 CPM", "x"),
        "feature": ("none", "Gene."),
        "group": ("none", "Group."),
        "sample": ("none", "Sample."),
    }
    data = box_data(kidney_like(), columns=columns, feature_label="Gene",
                    value_label="Expression", groups=("a", "b"))  # fmt: skip
    figure_stem = stem("clinical", "box_demo", 1, day=date(2026, 10, 14))
    text = FigureText("Four genes by group", reading=("Each box is one gene in one group.",))
    write_figure(data, out, figure_stem, title="Box", data_releases={}, text=text)
    page = ModulePage(
        title="M10 Demo: box", purpose="Purpose.", input="Input.", method="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Figures", (figure_block(data, figure_stem, text),)),),
        interpretation="Interpretation.", questions=(SelfCheck("A question.", "An answer."),),
        pitfalls=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    day = date(2026, 10, 14)
    module = write(page, out / f"{stem('clinical', 'box_overview', 1, day=day)}.html")
    pages = {"page": module, "alone": out / f"{figure_stem}.html"}
    tasks = [{"name": f"{k}_{w}", "url": p.as_uri(), "width": w, "height": 900, "act": MEASURE}
             for k, p in pages.items() for w in (375, 1280)]  # fmt: skip
    return _browse(tasks, tmp_path_factory.mktemp("downloads"))


@pytest.mark.parametrize("where", ["page", "alone"])
def test_a_phone_draws_the_box_whole_with_its_labels_apart(ran, where):
    """At 375 px the box fits its card, its turned labels apart, the labels over the boxes
    hidden; at 1,280 px the upright labels and the four p and r labels show.

    Mutations caught: no narrow layout (the box scrolls, 495 px in a 317 px card); the least
    width left at the upright layout's (it scrolls again)."""
    phone, desk = ran[f"{where}_375"], ran[f"{where}_1280"]
    assert phone["errors"] == [] and desk["errors"] == []
    seen = phone["acted"]
    assert seen["scroll"] <= seen["client"], seen
    assert len(seen["ticks"]) == 4 and _apart(seen["ticks"])  # noqa: PLR2004 - four genes
    assert seen["notes"] == 0
    assert desk["acted"]["notes"] == 4 and _apart(desk["acted"]["ticks"])  # noqa: PLR2004
