"""Opened alone at twice the size, a linked view across the card, and the export's files.

Rytis, 2026-10-01: "Open alone" drew a figure as small as the page did (B); the linked view's
parallel sets were squeezed beside the reading guide (C); and the PDF and SVG buttons gave the
build-time files whatever a reader changed (D). The browser half of D is
``tests/test_export_browser.py``; the page script's pure half is in ``test_page_script.py``.
"""

import base64
import dataclasses
import hashlib
import io
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from fontTools.ttLib import TTFont
from synth import counts, demo_band_map, demo_linked, demo_pca, demo_volcano
from test_viz import filtered_bars

from ctx.report import export, text
from ctx.report.builder import (
    FigurePage,
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    render,
)
from ctx.viz.save import figure_block
from ctx.viz.style import (
    ALONE_SCALE,
    AVERAGE_CHAR_EM,
    PAGE_PLOT_PX,
    PHONE_PLOT_PX,
    SIZES,
    THEMES,
    WIDTHS_MM,
    _is_length,
    on_page,
    plotly_templates,
    zoom,
)

ROOT = Path(__file__).resolve().parents[1]
READING = FigureText(title="Demo", reading=("Each mark is one row of the demo.",))
# Rytis, 2026-10-01: "Open alone" draws a figure 2x bigger than on the page. Stated here apart
# from style.ALONE_SCALE, so a test cannot take the constant it checks as its own oracle.
TWICE = 2
DEMOS = {
    "pca": lambda: demo_pca(counts()),
    "band_map": demo_band_map,
    "linked": demo_linked,
    "volcano": lambda: demo_volcano(counts(), width="double"),
    "bars": filtered_bars,
}


def _block(name: str):
    text = FigureText(title="Demo") if name == "volcano" else READING
    return figure_block(DEMOS[name](), "2026-10-01_qc_demo_v1", text)


def _widths(html: str) -> list[str]:
    """The CSS width of each plot's box, in the order the page draws them."""
    return re.findall(r'<div style="height:[^;"]*; width:([^"]*);">\s*<div id="ctx-', html)


def _drawn(html: str) -> list[tuple[str, dict, dict]]:
    """(div id, data, layout) of every Plotly.newPlot call in a page."""
    decoder, out = json.JSONDecoder(), []
    for match in re.finditer(r"Plotly\.newPlot\(\s*", html):
        values, at = [], match.end()
        for _ in range(3):
            value, at = decoder.raw_decode(html, at)
            values.append(value)
            at = re.compile(r"\s*,\s*").match(html, at).end()
        out.append(tuple(values))
    return out


def _named(parent: str, key: str) -> bool:
    """What Rytis named, stated here apart from style.zoom: a font's size, a marker's, a
    line's width, and the layout's width and height."""
    if key == "size":
        return parent == "marker" or parent.endswith("font")
    return (key, parent) in {("width", "line"), ("width", "layout"), ("height", "layout")}


def _compare(page, alone, factor, parent="", path=""):
    """Every length of ``alone`` is ``factor`` times the page's; everything else is equal.

    A length is what Rytis named (``_named``), or another length on screen ``style.zoom``
    scales, such as a margin or a tick's length.
    """
    if isinstance(page, dict):
        assert set(alone) >= set(page), path
        for key, value in page.items():
            if (_named(parent, key) or _is_length(parent, key)) and value is not None:
                scaled = [v * factor for v in value] if isinstance(value, list) else value * factor
                assert alone[key] == pytest.approx(scaled), f"{path}.{key}"
            else:
                _compare(value, alone[key], factor, key, f"{path}.{key}")
    elif isinstance(page, list):
        assert len(page) == len(alone), path
        for i, (a, b) in enumerate(zip(page, alone, strict=True)):
            _compare(a, b, factor, parent, f"{path}[{i}]")
    else:
        assert alone == page, path


def test_zoom_scales_the_lengths_on_screen_and_nothing_else():
    """Mutations: a "size" scaled under any parent (contours and bins are data); a "width"
    scaled on a bar (a data width); the margins left alone."""
    figure = {
        "data": [{"type": "bar", "width": 0.8, "marker": {"size": [3, 4], "line": {"width": 1}},
                  "contours": {"size": 1.0, "labelfont": {"size": 9}}, "xbins": {"size": 2},
                  "node": {"pad": 18, "thickness": 16}, "opacity": 0.5}],
        "layout": {"height": 400, "width": 600, "margin": {"l": 8, "t": 48, "pad": 2},
                   "title": {"x": 0.05, "font": {"size": 20}}, "xaxis": {"ticklen": 4,
                   "domain": [0, 0.5], "title": {"standoff": 8}},
                   "annotationdefaults": {"arrowhead": 0, "arrowwidth": 1}},
    }  # fmt: skip
    zoomed = zoom(figure, 2.0)
    trace, layout = zoomed["data"][0], zoomed["layout"]
    assert (trace["width"], trace["opacity"], trace["contours"]["size"]) == (0.8, 0.5, 1.0)
    assert trace["xbins"]["size"] == 2
    assert trace["marker"] == {"size": [6, 8], "line": {"width": 2}}
    assert trace["contours"]["labelfont"]["size"] == 18
    assert trace["node"] == {"pad": 36, "thickness": 32}
    assert (layout["height"], layout["width"], layout["margin"]) == (
        800, 1200, {"l": 16, "t": 96, "pad": 4})  # fmt: skip
    assert (layout["title"]["x"], layout["title"]["font"]["size"]) == (0.05, 40)
    assert layout["xaxis"] == {"ticklen": 8, "domain": [0, 0.5], "title": {"standoff": 16}}
    assert layout["annotationdefaults"] == {"arrowhead": 0, "arrowwidth": 2}
    assert figure["layout"]["height"] == 400, "zoom changed its input"


@pytest.mark.parametrize("name", list(DEMOS))
def test_alone_is_the_page_figure_at_twice_its_size(name):
    """Every font, marker, line, margin and the height twice the page's; the width twice the
    page's at full layout, or the box's where that is less (Rytis 2026-10-02), and never under
    twice a phone's plot or the least width the figure names at twice its type; data and
    positions the same; the theme toggle's templates at the same scale. Mutations:
    ALONE_SCALE not applied to the templates; the width left to the window alone; the width
    fixed at twice the page's, wider than a laptop's card; a marker size left out of the
    lengths; ALONE_SCALE 1.5 (which passed every test that took it as its own oracle, review
    2026-10-01)."""
    assert ALONE_SCALE == TWICE
    block = _block(name)
    alone = render(FigurePage("Demo", block))
    drawn = _drawn(alone)
    pages = [block.figure, *([block.below] if block.below is not None else [])]
    parts = ["beside", "across"]
    assert len(drawn) == len(pages) == len(_widths(alone))
    for (_, data, layout), page, part, width in zip(
        drawn, pages, parts, _widths(alone), strict=False
    ):
        page_json = json.loads(page.to_json())
        assert "width" not in page_json["layout"] and "width" not in layout
        cap = f"min(100%, {TWICE * PAGE_PLOT_PX[part]:.0f}px)"
        needs = (page_json["layout"].get("meta") or {}).get("least_width") or []
        least = max([TWICE * PHONE_PLOT_PX, *(TWICE * scaled + fixed for scaled, fixed in needs)])
        assert width == f"max({cap}, {least:.0f}px)"
        _compare(page_json["data"], data, TWICE, "data")
        _compare(page_json["layout"], layout, TWICE, "layout")
    (templates,) = re.findall(r'id="ctx-templates">(.*?)</script>', alone, re.S)
    _compare(json.loads(json.dumps(plotly_templates())), json.loads(templates), TWICE)


def test_the_page_widths_are_the_stylesheets():
    """The width a module page gives a plot, from the tokens its stylesheet reads: 1,072 px of
    content less a card's padding and border each side, less the reading guide and the gap
    beside it; Chromium measured 686 and 1,022 px at a 1,440 px window (2026-10-01).
    Mutation: the card's border left out, 2 px too wide."""
    assert PAGE_PLOT_PX == {"across": 1072 - 2 * 25, "beside": 1072 - 2 * 25 - 320 - 16}
    css = render(FigurePage("Demo", _block("pca")))
    assert ".fig{margin:var(--space-5) 0 var(--space-6);padding:var(--fig-pad);" in css
    assert "gap:var(--fig-gap)" in css and "--fig-pad: 24px; --fig-gap: 16px;" in css


def _module(*blocks) -> ModulePage:
    return ModulePage(
        title="M2 Demo: synthetic", question="Purpose.", concept="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Figures", tuple(blocks)),), practice="Interpretation.",
        think=(SelfCheck("A question.", "An answer."),), mistakes=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip


def test_linked_view_puts_the_sets_across_the_card_and_the_scatter_beside_the_guide():
    """The scatter in the plot column, the reading guide beside it, the parallel sets after
    both and before the caption, with room for their outer labels; each plot holds its own
    trace, as the brush addresses them. Mutations: the parallel sets drawn inside the plot
    column again; their margins left at the template's, which cut the outer labels off."""
    block = _block("linked")
    assert [t.type for t in block.figure.data] == ["scatter"]
    assert [t.type for t in block.below.data] == ["parcats"]
    assert (block.controls["scatter"], block.controls["parcats"]) == (0, 0)
    html = render(_module(block))
    plot = html.index('<div class="fig-plot">')
    reading, below = html.index('<aside class="reading"'), html.index('<div class="fig-below">')
    caption = html.index("<figcaption>")
    assert plot < html.index('id="ctx-figure-0"') < reading < below < caption
    assert below < html.index('id="ctx-below-0"') < caption
    # plotly sets the first dimension's levels left of the plot and the last's right of it.
    data, tick = demo_linked(), on_page(SIZES["tick"]) * AVERAGE_CHAR_EM
    ends = (data.dimensions[0], data.dimensions[-1])
    first, last = (max(map(len, data.frame[c].unique())) for c in ends)
    assert (
        block.below.layout.margin.l >= first * tick and block.below.layout.margin.r >= last * tick
    )


def test_vendored_files_are_their_manifests():
    """Each file's bytes are the SHA-256 its manifest gives, with its licences beside it; the
    fonts are Liberation Sans 2.1.5 under the OFL, as released, before a page renames its
    subset. Mutations: a manifest hash edited by one digit; a licence file the manifest names
    removed."""
    entries = export.manifest()
    assert set(entries) == {*export.SCRIPTS, *export.FONTS.values()}
    for name, entry in entries.items():
        data = (export.VENDOR / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == entry["sha256"], name
        assert entry["licence_files"] and all(
            (export.VENDOR / f).is_file() for f in entry["licence_files"])  # fmt: skip
    for style, name in export.FONTS.items():
        names = TTFont(export.VENDOR / name)["name"]
        assert names.getDebugName(6) == {"normal": "LiberationSans", "bold": "LiberationSans-Bold"}[
            style]  # fmt: skip
        assert names.getDebugName(5) == "Version 2.1.5"
        assert "SIL Open Font License, Version 1.1" in names.getDebugName(13)
    assert {entries[s]["licence"] for s in export.SCRIPTS} == {"MIT"}


def test_a_changed_vendored_file_is_refused(tmp_path, monkeypatch):
    """Mutation: the check skipped, so an edited library would ship on every page."""
    for path in export.VENDOR.iterdir():
        (tmp_path / path.name).write_bytes(path.read_bytes())
    script = tmp_path / export.SCRIPTS[1]
    script.write_bytes(script.read_bytes() + b";")
    monkeypatch.setattr(export, "VENDOR", tmp_path)
    export.vendored.cache_clear()
    try:
        with pytest.raises(export.VendorError, match="SHA-256"):
            export.vendored(export.SCRIPTS[1])
    finally:
        export.vendored.cache_clear()


def test_the_libraries_ride_once_and_only_on_a_page_with_controls():
    """Two figures with controls carry the libraries once; a page with none carries neither
    them nor the fonts. Mutation: the libraries inlined on every page."""
    banner = "jsPDF - PDF Document creation from JavaScript"
    with_controls = render(_module(_block("linked"), _block("volcano"), _block("pca")))
    without = render(_module(_block("pca"), _block("band_map")))
    assert with_controls.count(banner) == 1 and with_controls.count('id="ctx-export"') == 1
    assert "sourceMappingURL" not in with_controls
    assert banner not in without and 'id="ctx-export"' not in without
    assert 'data-export="pdf"' not in without


def _embedded_fonts(html: str) -> dict[str, bytes]:
    """The fonts a page carries for its PDF of a view, by style."""
    (config,) = re.findall(r'id="ctx-export">(.*?)</script>', html, re.S)
    return {style: base64.b64decode(data) for style, data in json.loads(config)["fonts"].items()}


def _characters_set(html: str) -> set[str]:
    """Every character a page's PDF of a view can set, read off the page itself: every string
    its plots draw and its controls carry, and the package's page text."""
    found = {c for _, data, layout in _drawn(html) for s in export.strings([data, layout])
             for c in s}  # fmt: skip
    for controls in re.findall(r'id="ctx-controls-\d+">(.*?)</script>', html, re.S):
        found |= {c for s in export.strings(json.loads(controls)) for c in s}
    content = ROOT / "src" / "ctx" / "report"
    for path in [*(content / "content").glob("*.toml"), content / "glossary.yaml"]:
        found |= set(path.read_text(encoding="utf-8"))
    return found - set("\n\t\r")


def test_the_fonts_keep_every_character_a_page_sets():
    """Rytis, 2026-10-01: the embedded Liberation Sans is cut to the characters the pages set.
    Each font a page carries keeps every character its controlled figures set, a cell type in
    Cyrillic among them, every character of the content and the glossary, and the stated margin
    (Latin, Greek, punctuation); the OFL text and the note beside it, once; and it is a fraction
    of the whole font. Its names are ``test_the_subset_is_renamed_and_keeps_its_notices``'s.
    Mutations: the page's own characters left out of the subset, so "Жук" had no glyph; the
    package's text left out of the rule, which lost the minus sign and the
    greater-than-or-equal sign; the margin dropped; the fonts' notice left off; the whole fonts
    inlined; numpy arrays of text skipped."""
    data = demo_linked()
    frame = data.frame.with_columns(cell_type=pl.col("cell_type").replace({"LDN": "Жук"}))
    block = figure_block(dataclasses.replace(data, frame=frame), "2026-10-01_qc_demo_v1", READING)
    html = render(_module(block, _block("volcano")))
    wanted = _characters_set(html)
    assert {"Ж", "у", "к", "\u2212", "\u2265", "ρ"} <= wanted
    margin = {chr(c) for low, high in export.SUBSET_MARGIN for c in range(low, high + 1)}
    for style, data in _embedded_fonts(html).items():
        font, whole = TTFont(io.BytesIO(data)), TTFont(export.VENDOR / export.FONTS[style])
        held, every = set(map(chr, font.getBestCmap())), set(map(chr, whole.getBestCmap()))
        assert (wanted | margin) & every <= held, sorted((wanted | margin) & every - held)[:10]
        assert len(data) < (export.VENDOR / export.FONTS[style]).stat().st_size / 4
    assert html.count("with Reserved Font Name Liberation") == 1
    assert html.count(export.manifest()[export.FONTS["normal"]]["note"]) == 1
    # A figure may hold its text in a numpy array, as plotly keeps what it was given.
    held = {"x": np.array(["Жук", "b"], dtype=object), "y": np.array([1.5]), "z": np.array(["ü"])}
    assert export.characters(held) == {"Ж", "у", "к", "b", "ü"}


# The fonts' Reserved Font Names, from their licence's first lines: "with Reserved Font Arimo,
# Tinos and Cousine" (Google) and "with Reserved Font Name Liberation" (Red Hat).
RESERVED = ("Liberation", "Arimo", "Tinos", "Cousine")
# Rytis, 2026-10-02: the subset takes its own names, stated here apart from export's constants,
# by name ID: family (1), unique (3), full (4) and PostScript (6).
RENAMED = {
    "normal": {1: "Ctx Sans", 3: "Ctx Sans Regular; Version 2.1.5", 4: "Ctx Sans Regular",
               6: "CtxSans-Regular"},
    "bold": {1: "Ctx Sans", 3: "Ctx Sans Bold; Version 2.1.5", 4: "Ctx Sans Bold",
             6: "CtxSans-Bold"},
}  # fmt: skip
# The name IDs that present a font's name to users (OpenType name table): family, unique, full
# and PostScript; typographic family and subfamily, compatible full, WWS family and subfamily,
# and the variations PostScript prefix.
PRIMARY = {1, 3, 4, 6, 16, 17, 18, 21, 22, 25}
# Copyright, licence and licence URL, which the OFL asks every copy to carry.
NOTICES = {0, 13, 14}
OFL_TITLE = "SIL OPEN FONT LICENSE Version 1.1 - 26 February 2007"


def _records(table) -> dict[tuple[int, int, int, int], str]:
    return {(r.nameID, r.platformID, r.platEncID, r.langID): r.toUnicode() for r in table.names}


def test_the_subset_is_renamed_and_keeps_its_notices():
    """Rytis, 2026-10-02: the page's subset is a Modified Version, which may not keep a Reserved
    Font Name (OFL condition 3), so it takes its own. Read off the page's own fonts, on every
    platform's records: no Reserved Font Name in a name presented to users; the family, unique,
    full and PostScript names Ctx Sans's, each style's own; the copyright, the licence and its
    URL, and every other record, the original's; the OFL text once on the page, with the note
    that the face is a Modified Version of Liberation Sans. Mutations: rename() never called;
    only the Windows records renamed, so the Mac family still read Liberation Sans; the name
    records cut to the subsetter's default, which dropped the licence, its URL and the Mac
    records; the fonts' note left off the page; the full name set to the family alone, which
    gave regular and bold the one full name "Ctx Sans" (review, 2026-10-02); the unique name
    left without its version."""
    licence = (export.VENDOR / "LICENSE.liberation-fonts.txt").read_text(encoding="utf-8")
    heads = [line for line in licence.splitlines() if "with Reserved Font" in line]
    assert all(any(name in line for line in heads) for name in RESERVED)
    html = render(_module(_block("linked")))
    for style, data in _embedded_fonts(html).items():
        records = _records(TTFont(io.BytesIO(data))["name"])
        original = _records(TTFont(export.VENDOR / export.FONTS[style])["name"])
        assert set(records) == set(original) and {key[0] for key in records} >= NOTICES
        for key, value in records.items():
            if key[0] in PRIMARY:
                assert not any(name in value for name in RESERVED), (key, value)
            else:
                assert value == original[key], key
        for name_id, name in RENAMED[style].items():
            assert {v for k, v in records.items() if k[0] == name_id} == {name}, name_id
    assert html.count(OFL_TITLE) == 1 and html.count("with Reserved Font Name Liberation") == 1
    (note,) = {export.manifest()[name]["note"] for name in export.FONTS.values()}
    assert html.count(note) == 1
    assert "Ctx Sans" in note and "Modified Version of Liberation Sans" in note


BUILD = """
import hashlib, json
from ctx.report import export
config = json.loads(str(export.page_assets(frozenset("\u0416\u0443\u043a"))["config"]))
print(json.dumps({s: hashlib.sha256(f.encode()).hexdigest() for s, f in config["fonts"].items()}))
"""


def test_two_builds_give_the_same_fonts():
    """A page's fonts are the same bytes from two builds, each in a fresh interpreter with its
    own hash seed, so a page rebuilt from the same characters is the same page. Mutation: the
    unique name stamped with the build's time."""
    source = str(Path(export.__file__).parents[2])
    env = {**os.environ, "PYTHONPATH": os.pathsep.join([source, os.environ.get("PYTHONPATH", "")])}
    run = [sys.executable, "-c", BUILD]
    runs = [
        subprocess.run(
            run, capture_output=True, text=True, check=True, env={**env, "PYTHONHASHSEED": seed}
        ).stdout
        for seed in ("1", "2")
    ]
    assert runs[0] == runs[1] and len(json.loads(runs[0])) == len(export.FONTS)


# The copyright lines of every notice a vendored script carries: its own and its bundle's.
COPYRIGHTS = {
    "jspdf.umd.min.js": ("(c) 2010-2025 James Hall", "Copyright (c) 2015 Michaël Zasso",
                         "Copyright (c) 2026 Arjun Barrett",
                         "Copyright (C) 2014-2017 by Vitaly Puzrin and Andrei Tuputcyn"),
    "svg2pdf.umd.min.js": ("Copyright (c) 2015-2026 yWorks GmbH",
                           "Copyright (c) 2013-2015 by Vitaly Puzrin",
                           "Copyright Mathias Bynens", "Copyright (c) 2015 Taro Hanamura",
                           "Copyright (c) 2016 Keegan Street and others",
                           "Copyright (c) Microsoft Corporation."),
}  # fmt: skip


def test_every_vendored_notice_rides_on_a_page_with_controls():
    """MIT asks for each copyright and permission notice in every copy, and a page is one: each
    script's own and those of the packages its bundle holds, before its code, once. The
    svg2pdf.js bundle printed none (review, 2026-10-01). Mutation: page_assets() inlining a
    script without notice()."""
    html = render(_module(_block("linked"), _block("volcano")))
    for script, lines in COPYRIGHTS.items():
        entry = export.manifest()[script]
        for line in lines:
            assert line in "".join((export.VENDOR / f).read_text(encoding="utf-8")
                                   for f in entry["licence_files"]), line  # fmt: skip
            assert html.count(line) >= 1, line
        head = f"/*! {entry['name']} {entry['version']}, MIT, from {entry['source']}"
        assert html.count(head) == 1
        code = export.vendored(script).decode()[:200]
        assert html.index(head) < html.index(code)
    assert html.count("Permission is hereby granted") >= 11


def test_dash_patterns_are_one_definition_for_both_renderers():
    """matplotlib's rc takes its dashes from style.DASH_PATTERNS, and the export maps plotly's
    names to the same lengths. Mutations: apply_matplotlib not setting the patterns (they
    equal matplotlib's defaults, so the rc is spoiled first); the export's settings without
    them, which left plotly's own 9 on 9 off for "dash" and 3 on 3 off for "dot"."""
    import matplotlib as mpl  # noqa: PLC0415

    from ctx.viz.style import DASH_PATTERNS, apply_matplotlib, plotly_dashes  # noqa: PLC0415

    with mpl.rc_context():
        for name in DASH_PATTERNS:
            mpl.rcParams[f"lines.{name}_pattern"] = [9.0, 9.0]
        apply_matplotlib(THEMES["white"])
        for name, pattern in DASH_PATTERNS.items():
            assert list(mpl.rcParams[f"lines.{name}_pattern"]) == list(pattern)
    assert plotly_dashes() == {"dash": (3.7, 1.6), "dot": (1.0, 1.65),
                               "dashdot": (6.4, 1.6, 1.0, 1.6)}  # fmt: skip
    config = json.loads(str(export.page_assets()["config"]))
    assert config["dashes"] == {k: list(v) for k, v in plotly_dashes().items()}


def test_a_thinned_volcano_ships_the_points_the_page_leaves_out(monkeypatch):
    """The page draws every 3rd of the points no setting can call; its PDF and SVG of a view
    get the rest, so with the page's they are every point the figure files draw, and they
    take the figure files' caption, not the page's thinning sentence. Mutations: the drawn
    points shipped instead of the omitted ones; caption_text() passed to the export."""
    import base64  # noqa: PLC0415

    import numpy as np  # noqa: PLC0415

    from ctx.viz import scatter  # noqa: PLC0415

    monkeypatch.setattr(scatter, "THIN_MAX", 300)
    data = demo_volcano(counts(), width="double")
    block = figure_block(data, "2026-10-01_qc_demo_v1", FigureText(title="Demo"))
    ex = block.controls["export"]
    assert "every 3rd" in block.caption and ex["caption"] == data.caption
    ns = block.controls["traces"]["ns"]
    assert list(ex["thinned"]) == [str(ns)]
    xy = np.frombuffer(base64.b64decode(ex["thinned"][str(ns)]), "<f4").reshape(-1, 2)
    page = data.on_page().filter(pl.col("call") == "ns").select("lfc", "neg_log10_p").rows()
    every = data.frame.filter(pl.col("call") == "ns").select("lfc", "neg_log10_p").rows()
    assert len(xy) + len(page) == len(every) and len(xy) == data.frame.height - len(block.rows)
    drawn = sorted(np.float32(v) for pair in page for v in pair)
    whole = sorted(np.float32(v) for pair in every for v in pair)
    assert sorted([*drawn, *xy.ravel()]) == whole


def test_the_export_carries_the_labels_the_figure_file_prints():
    """Each label the white figure file prints, at its point, with the offset and alignment
    static() placed it at. Mutation: no labels in the export settings, so a moved volcano's
    files named no gene (review, 2026-10-01)."""
    from matplotlib.text import Annotation  # noqa: PLC0415

    data = demo_volcano(counts(), width="double")
    labels = _block("volcano").controls["export"]["labels"]
    fig = data.static(THEMES["white"])
    printed = [(n.get_text(), n.xyann, n.get_ha()) for n in fig.axes[0].texts
               if isinstance(n, Annotation)]  # fmt: skip
    assert len(printed) == 6
    assert [(p["text"], (p["dx"], p["dy"]), p["anchor"]) for p in labels] == printed
    named = data.frame.filter(pl.col("label") != "")
    for label in labels:
        row = named.filter(pl.col("label") == label["text"])
        assert (label["x"], label["y"]) == pytest.approx((row["lfc"][0], row["neg_log10_p"][0]))


def test_the_linked_export_keeps_room_for_the_outer_levels():
    """A painted view's PDF draws the scatter's layout; its parallel sets' outer levels need
    their own margins, at the figure files' 9 pt. Mutation: the sides left out, so a long
    level ran 35 pt off the page (review, 2026-10-01)."""
    from ctx.viz.metrics import text_width  # noqa: PLC0415

    data = demo_linked()
    sides = _block("linked").controls["export"]["sides"]
    first = max(text_width(v, SIZES["tick"]) for v in data.frame[data.dimensions[0]].unique())
    last = max(text_width(v, SIZES["tick"]) for v in data.frame[data.dimensions[-1]].unique())
    assert sides["l"] > first and sides["r"] > last
    assert sides["l"] - first == pytest.approx(sides["r"] - last)


@pytest.mark.parametrize("name", ["volcano", "bars", "linked"])
def test_export_settings_are_the_figure_files_size(name):
    """A changed view's PDF and SVG take the figure files' width, 85 or 170 mm, and for a
    scatter or bars their height too. Mutation: a size taken from the page's plot."""
    data = DEMOS[name]()
    settings = _block(name).controls["export"]
    assert settings["width_pt"] == pytest.approx(WIDTHS_MM[data.width] / 25.4 * 72)
    if name != "linked":
        size = data.static(THEMES["white"]).get_size_inches() * 72
        assert (settings["width_pt"], settings["height_pt"]) == pytest.approx(tuple(size))
    assert settings["zoom"] == 1 and settings["file"] == "2026-10-01_qc_demo_v1"
    assert "[[" not in settings["caption"]
    assert settings["note"] == text.word("ui", "export.note") and settings["note"]
