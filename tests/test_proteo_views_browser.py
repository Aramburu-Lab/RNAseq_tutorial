"""M7's density views (P7.10) in Chromium on the made-up roots (tests/proteo_root.py), at
1,280 px light and 375 px dark: F4 at each fold change a call may need, F7 by class or by
kind of call. Each view's numbers are computed here with ``classify`` from the figure's own
CSV, never read from the page's JSON. Each test names the mutations it catches.
"""

import json
import re
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from proteo_root import across, cross
from test_export_browser import CHROMIUM
from views_browser import walk

from ctx.de.stats import ALPHA, classify
from ctx.proteo import densities, folds
from ctx.report import text
from ctx.report.builder import unpacked

pytestmark = pytest.mark.skipif(CHROMIUM is None, reason="no Chromium; set CTX_CHROMIUM")
SIZES = (("wide_light", 1280, "light"), ("phone_dark", 375, "dark"))


@pytest.fixture(scope="module")
def m7(tmp_path_factory):
    """Both roots' M7 folders and every view walked at both sizes."""
    seen = {}
    for name, make in (("across", across), ("cross", cross)):
        root = make(tmp_path_factory.mktemp(f"m7_{name}"))
        out = root / "outputs" / "07_proteo"
        (page,) = out.glob("*_overview_v1.html")
        seen[name] = (out, walk(page, tmp_path_factory.mktemp(f"walk_{name}"), SIZES))
    return seen


def _counts(frame: pl.DataFrame, lfc_min: float) -> dict[str, int]:
    """The lead's numbers, each layer called here by classify, apart from the module's code."""
    rna = frame.select(classify(pl.col("rna_padj"), pl.col("rna_lfc"), ALPHA, lfc_min))
    protein = frame.select(classify(pl.col("protein_padj"), pl.col("protein_lfc"), ALPHA, lfc_min))
    r, p = rna.to_series().to_list(), protein.to_series().to_list()
    out = {f"rna_{d}": r.count(d) for d in ("up", "down")}
    out |= {f"protein_{d}": p.count(d) for d in ("up", "down")}
    both = [(a, b) for a, b in zip(r, p, strict=True) if a != "ns" and b != "ns"]
    return out | {"both": len(both), "opposite": sum(a != b for a, b in both)}


def _kind_words() -> dict[str, str]:
    return {k: text.word("proteo", f"figures.contrast.kinds.{k}") for k in folds.KINDS}


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_each_fold_change_regroups_the_genes_as_classify_does(m7, size):
    """Mutations caught: every view drawn at the build's fold change; the protein layer called
    at the RNA's fold change only; the build not the first view."""
    out, seen = m7["across"]
    assert seen[size]["errors"] == []
    acted = seen[size]["acted"]
    assert acted["scroll"] <= acted["box"]
    views = [v for v in acted["views"] if v["label"] == "Calls at log2 FC at least"]
    assert views
    words = _kind_words()
    frames = [pl.read_csv(p) for p in out.glob("*_proteo_*_contrast_*_v1.csv")
              if not p.stem.endswith("_calls_v1")]  # fmt: skip
    for view in views:
        assert len(view["options"]) == len(densities.FOLDS)
        lfcs = [float(o["line"].split("at least ")[1].split(",")[0]) if o["line"] else 1.0
                for o in view["options"]]  # fmt: skip
        assert lfcs[0] == 1.0 and sorted(lfcs) == list(densities.FOLDS)
        leads = [[text.fill("proteo", "figures.contrast.lead", **_counts(f, lfc)) for lfc in lfcs]
                 for f in frames]  # fmt: skip
        assert [o["lead"] for o in view["options"]] in leads
        frame = frames[leads.index([o["lead"] for o in view["options"]])]
        for lfc_min, option in zip(lfcs, view["options"], strict=True):
            n = densities.counted(densities.kinds(frame, ALPHA, lfc_min))
            drawn = {t["name"].removesuffix(" 50 %").removesuffix(" 90 %")
                     for t in option["traces"]}  # fmt: skip
            assert {words[k] for k in folds.KINDS if n[k]} <= drawn | {words["both"]}


@pytest.mark.parametrize("size", ["wide_light", "phone_dark"])
def test_f7_groups_by_class_or_by_kind_of_call(m7, size):
    """Mutations caught: the kind view drawn by class; the class marks dropped by kind."""
    _, seen = m7["cross"]
    assert seen[size]["errors"] == []
    acted = seen[size]["acted"]
    (view,) = [v for v in acted["views"] if v["label"] == "A density per"]
    words = set(_kind_words().values())
    by_class, by_kind = ({t["name"] for t in o["traces"]} for o in view["options"])
    assert not by_class & words and by_kind & words
    assert view["options"][1]["moved"] and not view["options"][0]["moved"]
    assert np.all([o["wide"] <= acted["box"] for o in view["options"]])


def test_the_grids_travel_in_16_bits_with_their_levels_on_that_scale(m7):
    """Mutations caught: a level left in density units over a 16-bit grid; the grid's
    maximum not its last step."""
    grids = []
    for out, _ in m7.values():
        (page,) = out.glob("*_overview_v1.html")
        for found in re.finditer(
            r'id="ctx-controls-\d+">(.*?)</script>', unpacked(page.read_text()), re.S
        ):
            controls = json.loads(found.group(1))
            if controls.get("kind") == "view":
                grids += [t["z"] for o in controls["options"] for t in o["figure"]["data"]
                          if t.get("type") == "contour"]  # fmt: skip
    assert grids and all(g["dtype"] == "u2" for g in grids)
    trace = {
        "type": "contour",
        "z": [[0.0, 1.0], [2.0, 4.0]],
        "contours": {"start": 2.0, "end": 2.0, "size": 1},
    }
    packed = densities._typed(trace)  # noqa: SLF001 - the packing rule itself
    z = densities._grid(packed["z"])  # noqa: SLF001
    assert z.max() == densities.QUANTA and packed["contours"]["start"] == densities.QUANTA / 2
    assert Path(page).stat().st_size > 0
