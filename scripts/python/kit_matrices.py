#!/usr/bin/env python3
# /// codebase-meta
# exempt: build_kit's helper, not a pipeline step; keeps raw GDC count matrices out of a kit
# ///
"""A raw GDC count matrix a kit leaves out, and the commands that rebuild it (build_kit.py).

A raw GDC count matrix never ships (Rytis 2026-10-09; kit v8's were 94 to 205 MB). The rule
reads the data, not a list of names: a gather table of the counts kind whose README names a
GDC folder is left out, with its README kept; the kit carries that folder's GDC manifest and
fetch script, and ``LEFT_OUT.json``, README.md and MANIFEST.json name each file with its bytes
and the commands that rebuild it, the gather command read from the table's own README. That
row is provenance, the command as run, and stays verbatim in the shipped README; the command
shown to a reader names a ``--dataset-file`` the codebase committed as ``conf/datasets/<file>``
(``kit_datasets.codebase_datasets``), and README.md names the builder's ctx root and any
dataset file of the builder's own, for the reader to replace (K.9). Where the kit commit's
copy is not the bytes the run read, a ``git show`` from the newest commit holding them comes
first among the commands and LEFT_OUT.json names that commit (K.12). A page
that links a left-out matrix keeps its link: the link check counts a declared file as present
(``check_kit.declared``). Every derived table ships, and a count matrix not from GDC ships.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any

import check_kit
from kit_datasets import Version, restore, shown_command

from ctx.io import naming
from ctx.runtime import shown_path

# A gather table of this kind, from a GDC folder, is a raw count matrix (counts, miRNA counts).
MATRIX_KIND = "_counts"
# The rows ctx's table READMEs carry (ctx.io.readme_companion): the folders read, the command.
# A sheet joining several projects names each folder after its project (gdc_run, D-076).
GDC_FOLDER_ROW = re.compile(r"^\| (?:[\w-]+ )?GDC (?:\w+ )?folder \| (.+?) \|$", re.MULTILINE)
COMMAND_ROW = re.compile(r"^\| Command \| `(.+)` \|$", re.MULTILINE)
# The dataset file's SHA-256 as gather read it, from the README's Data releases rows.
DATASET_SHA_ROW = re.compile(r"^\| sha256:([^|/ ]+\.toml) \| ([0-9a-f]{64}) \|$", re.MULTILINE)
FILE_ROW = re.compile(r"^\| `([^`]+)` \| (\d+) \| `([0-9a-f]{64})` \|$", re.MULTILINE)
# How a GDC folder's own fetch script is run again, into the folder it fetched.
FETCH = "OUT={folder} bash {script}"


def table_readme(path: Path) -> Path:
    """The README ctx writes beside a table, shared by its CSV and Parquet copies."""
    return path.with_name(path.name.split(".")[0] + "_README.md")


def _after(words: list[str], flag: str) -> str | None:
    """The value following ``flag`` in a split command, or None."""
    return next((v for f, v in zip(words, words[1:], strict=False) if f == flag), None)


def gdc_matrix(path: Path, codebase: dict[str, Version] | None = None) -> dict[str, Any] | None:
    """What a raw GDC count matrix's README records, or None for any other file.

    A gather table of the counts kind whose README names a GDC folder; its bytes and SHA-256
    from the README's own file table, the folders and the gather command from its rows, the
    command as ``shown_command`` shows it against ``codebase`` (``codebase_datasets``).
    """
    note = table_readme(path)
    kind = (naming.description(path.name) or "").endswith(MATRIX_KIND)
    if path.parent.name != "01_gather" or not kind or path == note or not note.is_file():
        return None
    written = note.read_text(encoding="utf-8")
    folders = [Path(f) for f in GDC_FOLDER_ROW.findall(written)]
    if not folders:
        return None
    command = COMMAND_ROW.search(written)
    recorded = dict(DATASET_SHA_ROW.findall(written))
    shown, own, earlier = (shown_command(command[1], recorded, codebase) if command
                           else (None, None, None))  # fmt: skip
    root = _after(command[1].split(), "--root") if command else None
    sizes = {name: (int(n), sha) for name, n, sha in FILE_ROW.findall(written)}
    return {"folders": folders, "command": shown, "own_dataset_file": own, "earlier": earlier,
            "root": shown_path(root) if root else None,
            "bytes": sizes.get(path.name, (path.stat().st_size, ""))[0],
            "sha256": sizes.get(path.name, (0, ""))[1] or check_kit.sha256(path)}  # fmt: skip


def gdc_sources(folder: Path) -> list[Path]:
    """A GDC folder's manifests and fetch script, which the kit carries for a left-out matrix."""
    found = sorted(folder.glob("gdc_manifest_*.txt"))
    script = sorted(folder.glob("_fetch_*.sh"))
    return [*found, *script] if found and script else []


def carry_matrices(out: Path, folder: str, base: Path,
                   matrices: dict[Path, dict[str, Any]]) -> list[dict[str, Any]]:  # fmt: skip
    """Copy each left-out matrix's GDC manifest and fetch script; its LEFT_OUT.json entries."""
    entries = []
    for p, m in sorted(matrices.items()):
        carried, fetches, earlier = [], [], m["earlier"]
        for source in m["folders"]:
            place = out / folder / "gdc" / source.name
            place.mkdir(parents=True, exist_ok=True)
            for f in gdc_sources(source):
                carried.append(Path(shutil.copy2(f, place / f.name)).relative_to(out).as_posix())
                if f.suffix == ".sh":
                    fetches.append(FETCH.format(folder=shown_path(source),
                                                script=carried[-1]))  # fmt: skip
        entries.append({
            "path": (Path(folder) / "outputs" / p.relative_to(base)).as_posix(),
            "bytes": m["bytes"], "sha256": m["sha256"], "carried": carried,
            "rebuild": [*fetches, *([restore(earlier)] if earlier else []), m["command"]],
            "builder_root": m["root"], "builder_dataset_file": m["own_dataset_file"],
            "dataset_file_commit": earlier.commit if earlier else None,
            "dataset_file_path": earlier.path if earlier else None,
        })  # fmt: skip
    return entries
