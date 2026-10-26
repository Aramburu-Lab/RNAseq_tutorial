#!/usr/bin/env python3
# /// codebase-meta
# exempt: ships inside the student kit; checks a copy of the kit with the standard library alone
# ///
"""Check a course kit: every file against MANIFEST.json, every link, and every name on Windows.

Standard library only, so it runs wherever Python does, with nothing installed; the kit
builder imports the same functions, so a kit is checked by one definition when it is built
and when a copy of it is checked.

Usage:
    python check_kit.py            # the kit this file sits in
    python check_kit.py KIT_FOLDER
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

# Every attribute through which a browser follows or fetches a URL.
ATTRIBUTES = ("href", "src", "srcset", "data", "poster", "xlink:href")
_CSS = re.compile(r"""url\(\s*['"]?([^'")\s]+)|@import\s+['"]([^'"]+)""")
# Windows' MAX_PATH is 260 with the terminating NUL. Explorer's "Extract all" unpacks into
# Downloads\<zip name>\, and the zip holds <kit>/ again; a 20-character user name is assumed.
WINDOWS_LIMIT = 259
WINDOWS_PREFIX = "C:\\Users\\" + "x" * 20 + "\\Downloads\\{kit}\\"
_FORBIDDEN = re.compile(r'[<>:"|?*\\\x00-\x1f]')
_RESERVED = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
             *(f"LPT{i}" for i in range(1, 10))}  # fmt: skip
# The kit's record of the files left out of it by rule, each with its rebuild commands.
LEFT_OUT = "LEFT_OUT.json"
# A packed part of a page (ctx.report.builder.PACKED): its HTML gzipped and base64'd in a
# script, drawn when the reader reaches it. Its links are links of the page.
PACKED = re.compile(
    r'<script type="application/octet-stream" id="ctx-packed-([0-9]+)">([^<]*)</script>'
)
# The modules that read a gene list from a link (ctx.report.entry.CONSUMERS), and the entry
# page each one's folder must hold: a link from another page names it, never a dated page.
CONSUMERS = ("function", "genecard", "clinical")
ENTRY = "index.html"
_MODULE_FOLDER = re.compile(r"[0-9]{2}_([a-z_]+)")
_DATED_PAGE = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}_.+[.]html")


class _Links(HTMLParser):
    """Every URL a page's tags or style carry, as (tag, attribute, value)."""

    def __init__(self) -> None:
        super().__init__()
        self.found: list[tuple[str, str, str]] = []
        self._style = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        """Keep link attributes, each srcset candidate, and URLs inside style attributes."""
        self._style = tag == "style"
        for key, value in attrs:
            if not value:
                continue
            if key == "srcset":
                # Candidates end at a comma and a space; a data: URI holds a bare comma.
                self.found += [
                    (tag, key, c.split()[0]) for c in re.split(r",\s+", value) if c.strip()
                ]
            elif key in ATTRIBUTES:
                self.found.append((tag, key, value))
            elif key == "style":
                self.found += [(tag, key, a or b) for a, b in _CSS.findall(value)]

    def handle_endtag(self, tag: str) -> None:
        """Leave a style element."""
        if tag == "style":
            self._style = False

    def handle_data(self, data: str) -> None:
        """URLs in a style element's CSS; a script's text is not read."""
        if self._style:
            self.found += [("style", "url", a or b) for a, b in _CSS.findall(data)]


def unpacked(html: str) -> str:
    """A page with each packed part's HTML in place of its script (twin of builder.unpacked)."""
    return PACKED.sub(lambda m: gzip.decompress(base64.b64decode(m[2])).decode(), html)


def links(page: Path) -> list[tuple[str, str, str]]:
    """The (tag, attribute, value) of every URL on one page, its packed parts' included."""
    parser = _Links()
    parser.feed(unpacked(page.read_text(encoding="utf-8", errors="replace")))
    return parser.found


def entry_problems(root: Path) -> list[str]:
    """Each folder of a gene-list-reading module that holds dated pages and no entry page.

    ``root`` is a ctx root (``outputs/NN_<module>``) or a kit (``<folder>/outputs/...``).
    """
    problems = []
    for folder in sorted({*root.glob("outputs/*"), *root.glob("*/outputs/*")}):
        found = _MODULE_FOLDER.fullmatch(folder.name)
        if not folder.is_dir() or not found or found[1] not in CONSUMERS:
            continue
        dated = any(_DATED_PAGE.fullmatch(p.name) for p in folder.glob("*.html"))
        if dated and not (folder / ENTRY).is_file():
            where = folder.relative_to(root).as_posix()
            problems.append(f"{where}: no {ENTRY}, which a link to this module names")
    return problems


def fetches(tag: str, attribute: str, value: str) -> bool:
    """Whether a URL leaves the kit: anything off the machine but a clicked web link.

    A ``file:`` URL names the builder's own disk, so even a clicked one fails elsewhere.
    """
    url = urlsplit(value)
    if url.scheme == "data" or not (url.scheme or url.netloc):
        return False
    clicked = tag == "a" and attribute == "href" and url.scheme in ("http", "https", "mailto")
    return not clicked


def target(page: Path, value: str) -> Path | None:
    """The file a relative URL names, resolved; None for a URL that names no local file."""
    url = urlsplit(value)
    if url.scheme or url.netloc or not url.path:
        return None
    return (page.parent / unquote(url.path)).resolve()


def declared(root: Path) -> set[Path]:
    """The files the kit declares left out (``LEFT_OUT.json``), resolved; a link to one resolves.

    GDC count matrices stay out of a kit (Rytis 2026-10-09); a page that links one keeps its
    link, and the declaration names the file, its bytes and the commands that rebuild it.
    """
    record = root / LEFT_OUT
    if not record.is_file():
        return set()
    entries = json.loads(record.read_text(encoding="utf-8"))["files"]
    return {(root / e["path"]).resolve() for e in entries}


def link_problems(root: Path, pages: list[Path]) -> list[str]:
    """Every URL that fetches from off the machine, or names a file missing or outside root.

    A file the kit declares left out (``declared``) counts as present.
    """
    root = root.resolve()
    absent = declared(root)
    problems = []
    for page in pages:
        where = page.resolve().relative_to(root)
        for tag, attribute, value in links(page):
            found = target(page, value)
            if fetches(tag, attribute, value):
                problems.append(f"{where}: loads {value[:80]}")
            elif (
                found is not None
                and found not in absent
                and not (found.is_file() and found.is_relative_to(root))
            ):
                problems.append(f"{where}: {value} does not resolve")
    return problems


def windows_problems(names: list[str], kit: str) -> list[str]:
    """Names a Windows copy of the kit cannot hold, each relative to the kit's folder."""
    prefix = len(WINDOWS_PREFIX.format(kit=kit)) + len(kit) + 1
    problems, seen = [], {}
    for name in names:
        parts = name.split("/")
        if prefix + len(name) > WINDOWS_LIMIT:
            problems.append(f"{name}: {prefix + len(name)} characters after extraction")
        if any(_FORBIDDEN.search(p) or p.endswith((".", " ")) for p in parts):
            problems.append(f"{name}: a character or ending Windows refuses")
        if any(p.split(".")[0].upper() in _RESERVED for p in parts):
            problems.append(f"{name}: a name Windows reserves")
        if (other := seen.setdefault(name.lower(), name)) != name:
            problems.append(f"{name}: differs from {other} only by case")
    return problems


def sha256(path: Path) -> str:
    """A file's SHA-256, read in blocks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_problems(kit: Path) -> list[str]:
    """Files missing, changed or unlisted against the kit's MANIFEST.json."""
    record = json.loads((kit / "MANIFEST.json").read_text(encoding="utf-8"))
    listed = {e["path"]: e for e in record["entries"]}
    present = {p.relative_to(kit).as_posix() for p in kit.rglob("*") if p.is_file()}
    problems = [f"{name}: listed, not present" for name in sorted(listed.keys() - present)]
    problems += [f"{name}: present, not listed" for name in sorted(present - listed.keys())
                 if name != "MANIFEST.json"]  # fmt: skip
    for name in sorted(listed.keys() & present):
        entry, path = listed[name], kit / name
        if path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
            problems.append(f"{name}: not the file the course built")
    return problems


def check(kit: Path) -> list[str]:
    """Every problem with one kit folder: its manifest, its links and its names."""
    names = sorted(p.relative_to(kit).as_posix() for p in kit.rglob("*") if p.is_file())
    return [
        *manifest_problems(kit),
        *link_problems(kit, sorted(kit.rglob("*.html"))),
        *entry_problems(kit),
        *windows_problems(names, kit.name),
    ]


def main() -> None:
    """Check the kit named on the command line, or the one this file sits in."""
    kit = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    problems = check(kit)
    print("\n".join(problems) or f"{kit.name}: every file, link and name checks out")
    sys.exit(1 if problems else 0)


# A dataset's hub and its module pages, by file name.
_HUB = re.compile(r"\d{4}-\d{2}-\d{2}_report_(?P<dataset>.+)_hub_v\d+\.html")
_OVERVIEW = re.compile(r"\d{4}-\d{2}-\d{2}_[a-z_]+?_(?P<dataset>.+)_overview_v\d+\.html")


def unlinked_pages(kept: list[Path]) -> list[str]:
    """Each module page its dataset's hub does not link: the hub was written before it.

    ``ctx run`` writes the hub last and ``ctx genecard`` runs after it, so Lung-tissue-LUAD's
    first hub linked no M9 page while the kit shipped them (2026-10-26).
    """
    hubs = {m["dataset"]: p for p in kept if (m := _HUB.fullmatch(p.name))}
    problems = []
    for dataset, hub in sorted(hubs.items()):
        linked = {target(hub, v) for _, _, v in links(hub)}
        for p in kept:
            m = _OVERVIEW.fullmatch(p.name)
            if m and m["dataset"] == dataset and p not in linked:
                problems.append(f"{p}: {hub.name} does not link it; write the hub again with "
                                f"`ctx report --dataset {dataset}`")  # fmt: skip
    return problems


if __name__ == "__main__":
    main()
