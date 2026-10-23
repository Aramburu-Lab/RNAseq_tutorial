"""The course home page of a student kit: every module in course order, each dataset's page.

Built from the dataset folders a kit already holds, so every link it makes points inside the
kit and works wherever the kit is unpacked (Phase 2K). Text is in ``content/course.toml``.
"""

import html
import os
import re
from collections.abc import Sequence
from pathlib import Path

from ctx import editions
from ctx.io.naming import MODULES
from ctx.report import text
from ctx.report.builder import CourseLink, CoursePage, CourseRow, KeyNumber, shown, write
from ctx.report.hub import NOT_MODULES, PAGE, TITLES

# The page a student double-clicks. A script (.command, .cmd) runs through the student's own
# shell, whose startup can eat keys: oh-my-zsh's update prompt took the "/" of Terminal's
# typed path on 2026-09-30. An HTML file opens in the browser on every system with no shell.
INDEX = "START.html"
_TITLE = re.compile(r"<title>(.*?)</title>", re.S)


def label(page: Path) -> str:
    """The dataset a page names after the colon of its title; M0 names none."""
    found = _TITLE.search(page.read_text(encoding="utf-8", errors="replace"))
    title = html.unescape(found.group(1)) if found else ""
    return title.partition(": ")[2] or text.word("course", "machine")


def newest(folder: Path, suffix: str) -> list[Path]:
    """The newest page per dataset whose name ends ``_<suffix>_v<N>.html``, by name.

    Args:
        folder: A module folder, such as ``outputs/05_de``.
        suffix: ``overview``, ``hub``, or ``environment`` for M0, which names no dataset.

    Returns:
        One page per dataset, sorted by file name.
    """
    best: dict[str, tuple[str, int, Path]] = {}
    for p in folder.glob("*.html"):
        m = PAGE.fullmatch(p.name)
        if not m or not (m["what"] == suffix or m["what"].endswith(f"_{suffix}")):
            continue
        key = (m.group(1), int(m["version"]), p)
        if m["what"] not in best or key[:2] > best[m["what"]][:2]:
            best[m["what"]] = key
    return sorted((v[2] for v in best.values()), key=lambda p: p.name)


def _links(
    kit: Path, roots: Sequence[tuple[Path, bool]], folder: str, suffix: str
) -> list[CourseLink]:
    """One link per dataset page in one module folder across every root."""
    return [
        CourseLink(label(p), Path(os.path.relpath(p, kit)).as_posix(), real)
        for root, real in roots
        # M0 checks the computer, the same in every root: made-up roots would only repeat it.
        if real or suffix != "environment"
        for p in newest(root / "outputs" / folder, suffix)
    ]


def page(kit: Path, roots: Sequence[tuple[Path, bool]]) -> CoursePage:
    """The course home page of a kit.

    Args:
        kit: The kit's folder; the page is written at its top.
        roots: Each dataset root inside the kit, with whether its data are real.

    Returns:
        The page, one row per module of the selected edition (``ctx.editions``), in course
        order.
    """
    rows, data, edition = [], [], editions.active()
    for module, number in MODULES:
        if module in NOT_MODULES or not edition.has_module(module):
            continue
        suffix = "environment" if module == "setup" else "overview"
        tag, _, name = TITLES[module].partition(" ")
        found = _links(kit, roots, f"{number:02d}_{module}", suffix)
        links = [link for link in found if shown(kit / link.href)]
        if found and not links:
            continue  # every page of the module is a lite page with nothing to show
        if module == "setup":
            # One computer built every root's pages, so M0 shows its check once: the newest,
            # by the date its file name carries (Rytis 2026-10-20). max keeps the first of a
            # tie, the first real root in the kit's order (K.17).
            links = [max(links, key=lambda link: _dated(Path(link.href)))] if links else []
        row = CourseRow(tag, name, tuple(links))
        rows.append(row)
        # M0 checks the computer and holds no data, so the real-data count leaves it out.
        if module != "setup":
            data.append(row)
    report = f"{dict(MODULES)['report']:02d}_report"
    hubs = _links(kit, roots, report, "hub")
    counts = (
        ("modules", sum(bool(r.links) for r in rows), len(rows)),
        ("real", sum(any(link.real for link in r.links) for r in data), len(data)),
    )
    keys = tuple(
        KeyNumber(text.word("course", f"keys.{key}"),
                  text.fill("course", "keys.value", n=n, total=total))
        for key, n, total in counts
    )  # fmt: skip
    # The newest glossary any root holds; each root's pages link their own copy (D-048).
    book, guide = (_newest_across(kit, roots, report, what) for what in ("glossary", "guide"))
    return CoursePage(text.word("course", "title"), tuple(rows), tuple(hubs), keys, book, guide)


def _newest_across(kit: Path, roots: Sequence[tuple[Path, bool]], folder: str, what: str) -> str:
    """The newest page of one kind any root holds, linked from the kit's top, or empty.

    The glossary (D-048) and the custom-data guide: each root's pages link their own copy.
    """
    pages = [p for root, _ in roots for p in newest(root / "outputs" / folder, what)]
    found = max(pages, key=_dated, default=None)
    return Path(os.path.relpath(found, kit)).as_posix() if found else ""


def _dated(page: Path) -> tuple[str, int]:
    """A dated page's day and version, the order of its copies."""
    m = PAGE.fullmatch(page.name)
    return (m.group(1), int(m["version"])) if m else ("", 0)


def write_home(kit: Path, roots: Sequence[tuple[Path, bool]]) -> Path:
    """Write the course home page at the top of a kit.

    Args:
        kit: The kit's folder.
        roots: Each dataset root inside the kit, with whether its data are real.

    Returns:
        The written ``START.html``.
    """
    return write(page(kit, roots), kit / INDEX)
