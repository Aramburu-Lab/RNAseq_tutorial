"""No run overwrites another's output (TODO.md ED.5b, Rytis 2026-10-22): writers write here.

A same-day rerun of a module writes its files under the same names. Each file is written to
``.staging/<name>`` beside its target first. Where the target exists with other bytes, it moves
into ``superseded_<run stamp>/`` beside it before the new file takes its place; with the same
bytes the staged copy is dropped, so a rerun that changes nothing moves nothing. Readers glob a
module folder's top level and the naming lint skips hidden folders, so neither folder is read.
"""

import contextlib
import filecmp
from collections.abc import Callable
from datetime import datetime
from pathlib import Path

# One stamp per process, so every file one run moves aside lands in the same folder.
RUN_STAMP = datetime.now().strftime("%Y-%m-%dT%H%M%S")
STAGING = ".staging"
# The files this process has written: a run rewriting its own file replaces it, since only an
# earlier run's output is kept (M9 writes a page twice in one run).
_OURS: set[Path] = set()


def written(path: Path, write: Callable[[Path], object]) -> Path:
    """Write ``path`` through ``write(staged_path)``, keeping any older, different copy.

    Args:
        path: The file's final place.
        write: Writes the file to the path it is given; its return value is ignored.

    Returns:
        ``path``.
    """
    staged = path.parent / STAGING / path.name
    staged.parent.mkdir(parents=True, exist_ok=True)
    try:
        write(staged)
        _settle(path, staged)
    except BaseException:
        staged.unlink(missing_ok=True)  # a failed write leaves nothing new behind
        raise
    finally:
        # The staging folder goes once empty: readers list a module folder's entries as files.
        with contextlib.suppress(OSError):
            staged.parent.rmdir()
    return path


def _settle(path: Path, staged: Path) -> None:
    """Put ``staged`` at ``path``: an earlier run's different copy moves aside, a same one stays."""
    ours = path.resolve() in _OURS
    _OURS.add(path.resolve())
    if path.exists() and not ours:
        if filecmp.cmp(path, staged, shallow=False):
            staged.unlink()
            return
        aside = path.parent / f"superseded_{RUN_STAMP}"
        aside.mkdir(exist_ok=True)
        target, n = aside / path.name, 1
        while target.exists():  # one run writing one name twice keeps both earlier copies
            n += 1
            target = aside / f"{path.stem}.{n}{path.suffix}"
        path.rename(target)
    staged.replace(path)
