"""Content-keyed folders: one per distinct set of inputs, reused only while intact (D-037).

An index, a trimmed run, an alignment or a count is built into ``<prefix>_<key12>.part`` and
renamed to ``<prefix>_<key12>`` only when its build returned, with ``ctx_manifest.json`` listing
the inputs and every file's size. A changed input is a new key and a new folder beside the old
one, so a student who reruns with another counting option reuses the alignment, and no folder
an earlier page describes is overwritten. A folder that no longer matches its manifest, and a
``.part`` a crashed build left, are moved aside, never deleted and never reused.
"""

import hashlib
import json
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

MANIFEST = "ctx_manifest.json"
KEY_CHARS = 12


def key(inputs: Mapping[str, Any]) -> str:
    """SHA-256 of the inputs as canonical JSON: sorted keys, no spaces.

    Raises:
        TypeError: If an input is not plain JSON; a path must be given by its SHA-256.
        ValueError: If an input is NaN or infinite.
    """
    text = json.dumps(inputs, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode()).hexdigest()


def keyed(
    parent: Path, prefix: str, inputs: Mapping[str, Any], build: Callable[[Path], None]
) -> tuple[Path, bool]:
    """The folder for these inputs: reused when intact, else built.

    Args:
        parent: Where the keyed folders sit.
        prefix: The folder's name before the key, e.g. ``"star"``.
        inputs: Everything the content depends on, as plain JSON: tool and version, the
            arguments, and each input file's SHA-256.
        build: Writes the content into the folder it is given.

    Returns:
        The folder, absolute, and whether it was built now (False: reused).
    """
    digest = key(inputs)
    # Absolute: a builder runs its tool inside the folder, where a relative path misses.
    folder = parent.resolve() / f"{prefix}_{digest[:KEY_CHARS]}"
    if _intact(folder, digest):
        return folder, False
    part = folder.with_name(folder.name + ".part")
    for stale in (folder, part):
        if stale.exists():
            _set_aside(stale)
    part.mkdir(parents=True)
    build(part)
    files = {str(p.relative_to(part)): p.stat().st_size
             for p in sorted(part.rglob("*")) if p.is_file()}  # fmt: skip
    (part / MANIFEST).write_text(
        json.dumps({"key": digest, "inputs": inputs, "files": files}, indent=1, sort_keys=True)
        + "\n"
    )
    part.rename(folder)
    return folder, True


def _intact(folder: Path, digest: str) -> bool:
    """Whether the folder holds a finished build for this key.

    Its manifest must name the key, and every file it lists must be there at its listed size.
    """
    manifest = folder / MANIFEST
    if not manifest.is_file():
        return False
    try:
        recorded = json.loads(manifest.read_text())
    except json.JSONDecodeError:
        return False
    if recorded.get("key") != digest:
        return False
    return all(
        (folder / name).is_file() and (folder / name).stat().st_size == size
        for name, size in recorded.get("files", {}).items()
    )


def _set_aside(path: Path) -> Path:
    """Rename a folder that cannot be reused to the first free ``<name>.stale<N>``."""
    n = 1
    while (aside := path.with_name(f"{path.name}.stale{n}")).exists():
        n += 1
    path.rename(aside)
    return aside
