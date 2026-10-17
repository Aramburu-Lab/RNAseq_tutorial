#!/usr/bin/env python3
# /// codebase-meta
# exempt: build_kit's helper, not a pipeline step; names the dataset file a rebuild command reads
# ///
"""The dataset file a kit's rebuild command names: the codebase's committed bytes (D-075, K.12).

A gather README records its ``--dataset-file`` path and SHA-256. Bytes some commit holds are
the codebase's, whichever clone ran them, so the command names ``conf/datasets/<file>``. Bytes
the kit commit's copy does not hold come from the newest commit holding them, written beside
the command by ``git show <commit>:conf/datasets/<file> > <file>`` and named as ``<file>``:
tcga_kirc.toml changed at 5c524fe after ~/ctx_kidney's 2026-10-14 gather read it.
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path
from typing import NamedTuple

from ctx.runtime import shown_path

REPO = Path(__file__).resolve().parents[2]
DATASETS = "conf/datasets"
# Starts a commit line in the log, git's %x00, so a commit and its raw lines read apart; argv
# cannot carry the byte itself.
COMMIT = "\x00"


class Version(NamedTuple):
    """One committed version of a dataset file.

    Attributes:
        path: Its path in the repository.
        commit: The newest commit holding these bytes, or None when the kit commit's copy at
            ``path`` holds them.
    """

    path: str
    commit: str | None


def _git(repo: Path, *args: str, given: bytes | None = None) -> bytes:
    """What git prints for one command in ``repo``."""
    return subprocess.run(["git", "-C", str(repo), *args], input=given, capture_output=True,
                          check=True).stdout  # fmt: skip


def codebase_datasets(repo: Path = REPO) -> dict[str, Version]:
    """SHA-256 of every version of every dataset file ``repo`` committed, to that version.

    Bytes a commit holds are the repository's, from whichever clone or checkout ran them; a
    reader's own file, even one saved in ``conf/datasets/`` of a clone, has bytes no commit
    holds. A path test would take that file for the codebase's, and needs the tree on disk.
    """
    # Every parent of a merge too: a version can first appear in a merge (tcga_kirc.toml at
    # 789f8ad, which ~/ctx_kidney's gather read), and default log follows one parent. Newest
    # first, so the first commit to show a version is the newest holding it.
    raw = _git(repo, "log", "--full-history", "-m", "--format=%x00%H", "--raw",
               "--no-abbrev", "HEAD", "--", f"{DATASETS}/*.toml")  # fmt: skip
    blobs: dict[str, tuple[str, str]] = {}
    commit = ""
    for line in raw.decode().splitlines():
        meta, *paths = line.split("\t")
        if line.startswith(COMMIT):
            commit = line[1:]
        elif line.startswith(":") and meta.split()[3].strip("0"):
            blobs.setdefault(meta.split()[3], (paths[-1], commit))
    head = {}
    for line in _git(repo, "ls-tree", "HEAD", "--", f"{DATASETS}/").decode().splitlines():
        meta, path = line.split("\t")
        head[path] = meta.split()[2]
    batch = _git(repo, "cat-file", "--batch", given="".join(f"{b}\n" for b in blobs).encode())
    found, at = {}, 0
    for blob, (path, commit) in blobs.items():  # cat-file answers in the order asked
        end = batch.index(b"\n", at)
        size = int(batch[at:end].split()[2])
        sha = hashlib.sha256(batch[end + 1 : end + 1 + size]).hexdigest()
        found[sha] = Version(path, None if head.get(path) == blob else commit)
        at = end + size + 2
    return found


def shown_command(
    command: str, recorded: dict[str, str] | None = None,
    codebase: dict[str, Version] | None = None,
) -> tuple[str, str | None, str | None]:  # fmt: skip
    """A command as a README recorded it, with ``ctx`` for its entry point and ~ for home.

    A ``--dataset-file`` whose recorded SHA-256 (``recorded``, file name to hash) is a version
    the codebase committed (``codebase``) is shown by its repository path, or, where the kit
    commit's copy differs, by its file name after the command that writes those bytes there.

    Returns:
        The command shown; the dataset file it names that is the builder's own, or None; and
        the earlier version of the codebase's it names, for ``restore``, or None.
    """
    words = command.split()
    shown, own, earlier = [words[0]], None, None
    if words[0].endswith(("__main__.py", "/ctx")):
        shown = ["ctx"]
    for flag, word in zip(words, words[1:], strict=False):
        mapped = None
        if flag == "--dataset-file":
            version = (codebase or {}).get((recorded or {}).get(Path(word).name, ""))
            own = None if version else shown_path(word)
            if version is not None:
                mapped = Path(version.path).name if version.commit else version.path
                earlier = version if version.commit else None
        shown.append(mapped or (shown_path(word) if word.startswith("/") else word))
    return " ".join(shown), own, earlier


def restore(version: Version) -> str:
    """The command that writes an earlier version's bytes to its file name, from a clone."""
    return f"git show {version.commit}:{version.path} > {Path(version.path).name}"
