"""M3 on nf-core's alignment in Chromium: the core edition's lite page and the full page draw
every figure, raise no console error and scroll no sideways, at 1,280 and 375 px (D-111, D-116).

Skipped where no Chromium is found; set ``CTX_CHROMIUM``.
"""

import pytest
from test_export_browser import CHROMIUM, _browse
from test_lite_browser import LOOK
from test_process_nfcore import FIGURES, core, full  # noqa: F401 - the session builds

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
SIZES = {"wide": 1280, "phone": 375}


@pytest.fixture(scope="module")
def looked(full, core, tmp_path_factory):  # noqa: F811 - the fixtures imported above
    """Each edition's page at each width, as LOOK reads it."""
    tasks = []
    for edition, (root, _) in (("full", full), ("core", core)):
        (page,) = (root / "outputs" / "03_process").glob("*_overview_v1.html")
        tasks += [{"name": f"{edition}_{size}", "url": page.as_uri(), "width": width,
                   "height": 900, "scheme": "light", "act": LOOK}
                  for size, width in SIZES.items()]  # fmt: skip
    return _browse(tasks, tmp_path_factory.mktemp("downloads"))


@pytest.mark.parametrize("edition", ["full", "core"])
@pytest.mark.parametrize("size", list(SIZES))
def test_every_figure_is_drawn_without_sideways_scroll(looked, edition, size):
    """Each of the eleven figures is a drawn Plotly figure; no console error. Mutation: a
    contract the page cannot draw, which leaves its plot empty."""
    done = looked[f"{edition}_{size}"]
    seen = done["acted"]
    assert done["errors"] == []
    assert seen["figures"] == seen["drawn"] == len(FIGURES)
    assert seen["scroll"] <= seen["page"]
    if edition == "core":
        assert seen["sections"] == ["results"]
