"""A page's own PDF and SVG of a changed view (D-033, amended 2026-10-01).

Once a reader has moved a figure's page controls from their defaults, its PDF and SVG buttons
draw the view as it stands, in the browser, at the sizes and fonts of the figure's own files
(``templates/_export.js.j2``); at the defaults they link the matplotlib files. For that a page
inlines jsPDF, svg2pdf.js and Liberation Sans from ``vendor/``, once, and only when one of its
figures has page controls. ``vendor/MANIFEST.toml`` names each file's source, version, licence
and SHA-256, and a file whose bytes differ is refused. The fonts ride subset to the characters
the page can set (``subset_font``, Rytis 2026-10-01): whole, they were 1,100,224 of the
1,630,308 bytes the export added to a page; cut, 195,608 of 730,540 (measured 2026-10-02).
The subset is a Modified Version under the OFL and is renamed Ctx Sans (Rytis 2026-10-02).
"""

import base64
import hashlib
import io
import json
import re
import tomllib
from collections.abc import Iterable, Iterator, Mapping
from functools import cache
from pathlib import Path
from typing import Any

from markupsafe import Markup

from ctx.report import text
from ctx.report.marks import plain
from ctx.viz.style import (
    CAPTION_BOTTOM_PT,
    CAPTION_INSET,
    FONT_CHAIN,
    LINE_SPACING,
    LINE_WIDTHS,
    PAGE_PX_PER_PT,
    SIZES,
    THEMES,
    journal_template,
    plotly_dashes,
)

VENDOR = Path(__file__).parent / "vendor"
SCRIPTS = ("jspdf.umd.min.js", "svg2pdf.umd.min.js")
# Liberation Sans, Arial's metric twin, is the face the PDF embeds, first in the chain, from
# these vendored files. Cut to a page's characters it is a Modified Version, which the OFL
# (condition 3) forbids to keep the Reserved Font Name Liberation, so the subset is renamed
# (Rytis 2026-10-02): FAMILY, the family and full names; FACE, the family the drawing's text is
# set in for jsPDF; NAMES, each style's PostScript name, which the PDF carries so regular and
# bold differ.
FONTS = {"normal": "LiberationSans-Regular.ttf", "bold": "LiberationSans-Bold.ttf"}
FAMILY = "Ctx Sans"
FACE = "CtxSans"
NAMES = {"normal": "CtxSans-Regular", "bold": "CtxSans-Bold"}
# The bundles end by naming a source map the page does not carry; the inlined copy drops it.
SOURCE_MAP = re.compile(r"\n//# sourceMappingURL=\S+\s*$")
# The characters a page's embedded fonts keep (Rytis 2026-10-01): every character of the ctx
# package's text, which holds content/*.toml, the glossary and the code that still names axes,
# levels and captions (ctx.viz, ctx.de, ctx.qc); every character a page's figures with controls
# carry, their data's names among them (``characters``); and this margin of Unicode blocks, for
# a name or a sign the package has not set yet: printable ASCII, Latin-1 and Latin Extended-A,
# Greek, general punctuation. A character the full font lacks has no glyph to keep.
SUBSET_MARGIN = ((0x20, 0x7E), (0xA0, 0x17F), (0x370, 0x3FF), (0x2000, 0x206F))
PACKAGE = Path(__file__).parents[1]
PACKAGE_TEXT = (".py", ".toml", ".yaml", ".j2")


class VendorError(RuntimeError):
    """A vendored file is missing, unlisted, or not the bytes its manifest names."""


@cache
def manifest() -> dict[str, dict[str, Any]]:
    """``vendor/MANIFEST.toml``, by file name."""
    entries = tomllib.loads((VENDOR / "MANIFEST.toml").read_text())["file"]
    return {entry["file"]: entry for entry in entries}


@cache
def vendored(name: str) -> bytes:
    """A vendored file's bytes, checked against its manifest entry.

    Raises:
        VendorError: If the manifest does not list the file or its SHA-256 differs.
    """
    entry = manifest().get(name)
    if entry is None:
        raise VendorError(f"vendor/{name} is not in vendor/MANIFEST.toml")
    data = (VENDOR / name).read_bytes()
    if (found := hashlib.sha256(data).hexdigest()) != entry["sha256"]:
        raise VendorError(f"vendor/{name}: SHA-256 {found}, the manifest says {entry['sha256']}")
    return data


def notice(*names: str) -> str:
    """Vendored files' licence notices as a comment that minifiers and readers keep.

    MIT asks for its copyright and permission notice in every copy, and a page is one. The
    svg2pdf.js bundle carries no banner, so every script takes its manifest's licence files:
    its own, and those of what its bundle holds. Files that share a licence print it once.
    """
    entries = [manifest()[name] for name in names]
    files = dict.fromkeys(f for entry in entries for f in entry["licence_files"])
    texts = [(VENDOR / f).read_text(encoding="utf-8").strip() for f in files]
    heads = [f"{e['name']} {e['version']}, {e['licence']}, from {e['source']}" for e in entries]
    notes = [e["note"] for e in entries if "note" in e]
    body = "\n\n".join([*heads, *dict.fromkeys(notes), *texts])
    return "/*! " + body.replace("*/", "* /").replace("</", "<\\/") + " */"


def strings(node: Any) -> Iterator[str]:
    """Every string in a plotly figure's JSON or a page's settings: titles, labels and text."""
    if isinstance(node, str):
        yield node
    elif isinstance(node, Mapping):
        for value in node.values():
            yield from strings(value)
    elif isinstance(node, list | tuple):
        for value in node:
            yield from strings(value)
    elif getattr(getattr(node, "dtype", None), "kind", "") in ("O", "U", "S"):
        # A numpy array of text, by its dtype: numpy is not imported before the run's profile.
        yield from strings(node.tolist())


def characters(*nodes: Any) -> frozenset[str]:
    """Every character of every string in ``nodes``: what a page's PDF of a view can set."""
    return frozenset("".join(s for node in nodes for s in strings(node)))


@cache
def package_characters() -> frozenset[str]:
    """``SUBSET_MARGIN`` and every character of the ctx package's text files."""
    found = {chr(c) for low, high in SUBSET_MARGIN for c in range(low, high + 1)}
    for path in sorted(PACKAGE.rglob("*")):
        if path.suffix in PACKAGE_TEXT:
            found.update(path.read_text(encoding="utf-8"))
    return frozenset(found)


@cache
def subset_font(style: str, extra: frozenset[str] = frozenset()) -> bytes:
    """A style's vendored font cut to ``package_characters()`` and ``extra``, and renamed.

    fontTools' subsetter cuts it; ``rename`` gives it its own names. Every other name record
    stays, the copyright, the licence and its URL among them; the licence text with its
    Reserved Font Name rides in the fonts' notice (``page_assets``). Glyph outlines and hinting
    are the font's own; a missing character draws the notdef box, never nothing. No timestamp
    is renewed, so the same characters give the same bytes.
    """
    from fontTools import subset  # type: ignore[import-untyped]  # noqa: PLC0415 - heavy

    options = subset.Options(name_IDs=["*"], name_languages=["*"], name_legacy=True,
                             notdef_outline=True, recalc_timestamp=False)  # fmt: skip
    # FontForge's own timestamp table, which the subsetter cannot cut and would warn about.
    options.drop_tables += ["FFTM"]
    font = subset.load_font(io.BytesIO(vendored(FONTS[style])), options)
    held = font.getBestCmap()
    cutter = subset.Subsetter(options)
    cutter.populate(unicodes=sorted(c for c in map(ord, package_characters() | extra) if c in held))
    cutter.subset(font)
    rename(font, style)
    out = io.BytesIO()
    subset.save_font(font, out, options)
    return out.getvalue()


def rename(font: Any, style: str) -> None:
    """Give a subset font its own family, full, unique and PostScript names, on every platform.

    The OFL's Reserved Font Name restriction covers the primary font name as presented to users
    (condition 3); the copyright, trademark, version and licence records keep their text.
    """
    table = font["name"]
    subfamily, version = table.getDebugName(2), table.getDebugName(5)
    names = {1: FAMILY, 3: f"{FAMILY} {subfamily}; {version}", 4: f"{FAMILY} {subfamily}",
             6: NAMES[style]}  # fmt: skip
    for record in table.names:
        if record.nameID in names:
            record.string = names[record.nameID]


@cache
def page_assets(extra: frozenset[str] = frozenset()) -> dict[str, Markup]:
    """What a page with page controls inlines: the two libraries, and the export's settings.

    Args:
        extra: Characters the page's figures set beyond ``package_characters()``.

    Returns:
        ``script``, the libraries as one script, each after its licence notices, then the
        fonts' notice; ``config``, JSON of the fonts subset to the page's characters and
        renamed, in base64, the face's and each style's PostScript name, the journal template,
        the dash patterns and the regular face's widths of those characters, which
        ``ctxExport`` reads, safe inside a ``<script>`` element.
    """
    from ctx.viz.metrics import table  # noqa: PLC0415 - metrics reads this module's fonts

    script = "\n".join(notice(name) + "\n" + SOURCE_MAP.sub("\n", vendored(name).decode())
                       for name in SCRIPTS)  # fmt: skip
    script += "\n" + notice(*FONTS.values())
    fonts = {style: base64.b64encode(subset_font(style, extra)).decode() for style in FONTS}
    # A caption wraps by the figure files' own measure, so the view's breaks where theirs do.
    config = {"fonts": fonts, "names": NAMES, "face": FACE, "template": journal_template(),
              "dashes": plotly_dashes(),
              "metrics": table(package_characters() | extra)}  # fmt: skip
    text = json.dumps(config, allow_nan=False).replace("<", "\\u003c")
    return {"script": Markup(script), "config": Markup(text)}


def assets(prepared: Iterable[Mapping[str, Any]]) -> dict[str, Markup] | None:
    """The libraries and settings, on a page one of whose figures has page controls.

    The fonts keep every character those figures carry (``figures.prepare``). None for a page
    with none, which then carries neither.
    """
    exporting = [f for f in prepared if f.get("export")]
    if not exporting:
        return None
    carried = frozenset().union(*(f["characters"] for f in exporting))
    return page_assets(carried - package_characters())


def file_caption(caption: Any) -> str:
    """The caption a contract's own files print, from its ``caption`` attribute.

    A text as given (a scatter's, without the page's thinning line); a caption object, as a
    box's ``Caption`` or M9's ``OwnCaption``, as its ``text()`` on one line. The one place a
    caption object becomes text for a view's PDF and SVG: ``views.viewed`` passed M9b's
    ``OwnCaption`` on as text, and ``ctx genecard`` stopped on Step-1 (V2).
    """
    return caption if isinstance(caption, str) else str(caption.text())


def settings(
    inches: tuple[float, float], caption: Any, stem: str, layers: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """What a figure's controls carry for its PDF and SVG of a changed view.

    Args:
        inches: The size the drawing takes, as the figure's own files take theirs.
        caption: The contract's caption, text or object (``file_caption``): the one its own
            files print, not the page's; a glossary mark is its display text off the page.
        stem: The figure's file stem; the view is saved as ``<stem>_view.pdf`` and ``.svg``.
        layers: What the figure files draw that the page's plots lack, from the contract's
            ``export_layers()``: points the page thins, labels it leaves out, room it gives.
            A contract whose view's file is the page's drawing, as a scatter drawn as
            densities (D-118), gives its page caption here as ``caption``, which replaces it.

    Returns:
        The size in pt, the caption, the stem and the caption's type and insets; ``marker``
        turns a mark's page px back into pt, at the zoom a plot is drawn at, which its layout's
        ``meta.zoom`` gives; the labels' type, ink and leader width; ``note``, the sentence
        that closes the file's caption and says it is the page's drawing.
    """
    width, height = inches
    return {
        "width_pt": width * 72, "height_pt": height * 72, "caption": plain(file_caption(caption)),
        "note": text.word("ui", "export.note"),
        "file": stem, "size": SIZES["caption"], "spacing": LINE_SPACING, "inset": CAPTION_INSET,
        "bottom": CAPTION_BOTTOM_PT, "ink": THEMES["white"].ink, "marker": 1 / PAGE_PX_PER_PT,
        "family": ", ".join([*(f"'{f}'" for f in FONT_CHAIN), "sans-serif"]),
        "label_size": SIZES["annotation"], "leader": LINE_WIDTHS["line"], **(layers or {}),
    }  # fmt: skip
