"""The course home page of a student kit: every module in course order, each dataset's page.

Built from the dataset folders a kit already holds, so every link it makes points inside the
kit and works wherever the kit is unpacked (Phase 2K). Text is in ``content/course.toml``.
"""

import html
import os
import re
from collections.abc import Sequence
from pathlib import Path

from ctx.io.naming import MODULES
from ctx.report import text
from ctx.report.builder import CourseLink, CoursePage, CourseRow, KeyNumber, write
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
        The page, one row per module in course order.
    """
    rows, data = [], []
    for module, number in MODULES:
        if module in NOT_MODULES:
            continue
        suffix = "environment" if module == "setup" else "overview"
        tag, _, name = TITLES[module].partition(" ")
        row = CourseRow(tag, name, tuple(_links(kit, roots, f"{number:02d}_{module}", suffix)))
        rows.append(row)
        # M0 checks the computer and holds no data, so the real-data count leaves it out.
        if module != "setup":
            data.append(row)
    hubs = _links(kit, roots, f"{dict(MODULES)['report']:02d}_report", "hub")
    counts = (
        ("modules", sum(bool(r.links) for r in rows), len(rows)),
        ("real", sum(any(link.real for link in r.links) for r in data), len(data)),
    )
    keys = tuple(
        KeyNumber(text.word("course", f"keys.{key}"),
                  text.fill("course", "keys.value", n=n, total=total))
        for key, n, total in counts
    )  # fmt: skip
    return CoursePage(text.word("course", "title"), tuple(rows), tuple(hubs), keys)


def write_home(kit: Path, roots: Sequence[tuple[Path, bool]]) -> Path:
    """Write the course home page at the top of a kit.

    Args:
        kit: The kit's folder.
        roots: Each dataset root inside the kit, with whether its data are real.

    Returns:
        The written ``START.html``.
    """
    return write(page(kit, roots), kit / INDEX)
