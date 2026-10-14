"""A module folder's entry page, ``outputs/NN_<module>/index.html`` (GS.2, D-049).

Pages are dated and M5 is built before M6, so M5 cannot name M6's file. A module that reads a
gene list from a link writes this page on every build, so a link to the folder's entry page
always reaches the module's newest page, its gene list carried in the fragment. The kit's link
walk follows it (``scripts/python/build_kit.py``) and ``check_kit.py`` refuses a consuming
module's folder without one.
"""

import os
from pathlib import Path

from ctx.io.naming import ENTRY
from ctx.report import text
from ctx.report.builder import EntryPage, write
from ctx.report.course import label, newest
from ctx.report.hub import TITLES

# The modules that read a gene list from a link: M6's live ORA, M9's cards, M10's survival.
CONSUMERS = ("function", "genecard", "clinical")


def write_entry(folder: Path, module: str) -> Path:
    """Write ``index.html`` in a module folder, linking each dataset's newest overview page.

    Args:
        folder: The module's output folder, as ``outputs/11_genecard``.
        module: The module's token, one of ``CONSUMERS``.

    Returns:
        The entry page.

    Raises:
        ValueError: If the module reads no gene list, or the folder holds no overview page.
    """
    if module not in CONSUMERS:
        raise ValueError(f"{module} reads no gene list from a link; only {CONSUMERS} do")
    pages = newest(folder, "overview")
    if not pages:
        raise ValueError(f"{folder}: no overview page for the entry page to link")
    links = tuple((label(p), Path(os.path.relpath(p, folder)).as_posix()) for p in pages)
    title = text.fill("ui", "entry.title", module=TITLES[module])
    return write(EntryPage(title, links), folder / ENTRY)
