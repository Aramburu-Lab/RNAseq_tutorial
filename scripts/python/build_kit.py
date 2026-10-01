#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side; packs the generated pages and tables into the student kit (Phase 2K)
# ///
r"""Build the student kit: each root's current outputs and the course home page, START.html.

The kit is what a student downloads: the pages, figures and tables the course generated, and
nothing to install. ``START.html``, the course home page, links every dataset's pages. The
checks are ``check_kit.py``'s, which ships in the kit with ``check_kit.ps1`` for Windows:
every link resolves inside it, no page loads anything from the web, and every name fits
Windows. The kit is built in
``<out>.partial`` and renamed only once every check passes. ``MANIFEST.json`` holds every
file's SHA-256; ``--zip`` packs one folder.

Usage:
    pixi run python scripts/python/build_kit.py --out ~/ctx_kit/2026-09-30_course_kit_v2 \
        --real step1=~/ctx_tier1 --made-up preview=~/ctx_preview/synthetic --zip
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_kit  # noqa: E402 - the kit's own checker, beside this file

from ctx.report import text  # noqa: E402
from ctx.report.course import label, newest, write_home  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
# The same checks for a Windows reader without Python, reading the kit's own MANIFEST.json.
CHECK_PS1 = REPO / "scripts" / "powershell" / "check_kit.ps1"
# Reads, alignments and their indexes: a kit never carries them (Rytis 2026-09-30).
REFUSED = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".bam", ".bai", ".sam", ".cram", ".crai",
           ".sra")  # fmt: skip
# Above the largest output measured on 2026-09-30, Step-1's counts CSV at 39.5 MB.
CAP_BYTES = 50_000_000
DATED = re.compile(r"\d{4}-\d{2}-\d{2}_")
# Scripts for readers who look for START in their own system's form; each opens START.html.
# A double-clicked .command is typed into the user's shell, which runs it from the home
# folder, so each script finds the page beside itself; %~dp0 is the .cmd file's own folder.
START = {
    "START.command": ('#!/bin/sh\ncd "$(dirname "$0")" && open START.html\n', 0o755),
    "START.sh": ('#!/bin/sh\ncd "$(dirname "$0")" && xdg-open START.html\n', 0o755),
    "START.cmd": ('@echo off\r\nstart "" "%~dp0START.html"\r\n', 0o644),
}
FOLDER = re.compile(r"[a-z0-9_]+")
# Code a kit is built from; a change here, even an untracked file, marks the commit dirty.
CODE = ("src", "scripts", "conf", "pixi.toml", "pixi.lock", "pyproject.toml")


def current(outputs: Path) -> set[Path]:
    """The newest dated copy of every output in each module folder; older days are dropped."""
    best: dict[tuple[Path, str], Path] = {}
    for p in outputs.glob("*/*"):
        if p.is_file() and DATED.match(p.name):
            key = (p.parent, p.name[11:])
            if key not in best or p.name > best[key].name:
                best[key] = p
    return set(best.values())


def reached(pages: set[Path], inside: Path) -> tuple[set[Path], list[str]]:
    """Every file the pages link to, followed through linked pages, and every link problem.

    A URL that fetches from off the machine is a problem here too, before anything is copied.
    """
    seen, todo = set(), sorted(pages)
    while todo:
        page = todo.pop()
        if page in seen:
            continue
        seen.add(page)
        if page.suffix != ".html":
            continue
        for _, _, value in check_kit.links(page):
            found = check_kit.target(page, value)
            if found is not None and found.is_file() and found.is_relative_to(inside):
                todo.append(found)
    html = sorted(p for p in seen if p.suffix == ".html")
    return seen, check_kit.link_problems(inside, html)


def select(root: Path) -> tuple[list[Path], list[str]]:
    """A root's files for the kit, and every reason one of them cannot ship."""
    # Resolved, so a symlinked outputs/ holds the same paths the links resolve to.
    outputs = (root / "outputs").resolve()
    chosen = current(outputs)
    found, problems = reached({p for p in chosen if p.suffix == ".html"}, outputs)
    for p in sorted(chosen | found):
        if p.name.lower().endswith(REFUSED):
            problems.append(f"{p}: reads or alignments never ship")
        elif p.stat().st_size > CAP_BYTES:
            problems.append(f"{p}: {p.stat().st_size:,} bytes, over the {CAP_BYTES:,} cap")
    return sorted(chosen | found), problems


def _commit() -> str:
    """The codebase commit that built the kit, marked dirty when its code had any change.

    Untracked code counts: a builder not yet committed cannot be rebuilt from the commit.
    Untracked notes elsewhere in the tree do not build the kit, so they do not count.
    """

    def git(*args: str) -> str:
        return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True,
                              check=True).stdout.strip()  # fmt: skip

    changed = git("status", "--porcelain", "--untracked-files=all", "--", *CODE)
    return git("rev-parse", "--short", "HEAD") + ("-dirty" if changed else "")


def manifest(kit: Path, roots: list[tuple[str, bool]], commit: str) -> Path:
    """``MANIFEST.json``: every other file of the kit with its size and SHA-256."""
    out = kit / "MANIFEST.json"
    files = [
        {"path": p.relative_to(kit).as_posix(), "bytes": p.stat().st_size,
         "sha256": check_kit.sha256(p)}
        for p in sorted(kit.rglob("*")) if p.is_file()
    ]  # fmt: skip
    record = {
        "built": date.today().isoformat(),
        "commit": commit,
        "roots": [{"folder": f, "real": real} for f, real in roots],
        "files": len(files),
        "bytes": sum(f["bytes"] for f in files),
        "entries": files,
    }
    out.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    return out


def readme(kit: Path, roots: list[tuple[str, bool]], shipped: list[Path], commit: str) -> Path:
    """The kit's README.md, its text from ``content/course.toml``."""
    lines = []
    for folder, real in roots:
        pages = [p for d in sorted((kit / folder / "outputs").iterdir()) if d.is_dir()
                 for p in newest(d, "overview")]  # fmt: skip
        names = sorted({label(p) for p in pages})
        kind = text.word("course", "real" if real else "made_up").lower()
        lines.append(text.fill("course", "readme.dataset", folder=folder,
                               names=", ".join(names), kind=kind))  # fmt: skip
    size = f"{sum(p.stat().st_size for p in shipped) / 1e6:,.0f} MB"
    body = text.fill("course", "readme.text", files=f"{len(shipped):,}", size=size,
                     datasets="\n".join(lines), built=date.today().isoformat(),
                     commit=commit)  # fmt: skip
    out = kit / "README.md"
    out.write_text(body, encoding="utf-8")
    return out


def pack(kit: Path) -> Path:
    """One zip holding the kit's folder; file modes are kept, so the scripts stay runnable."""
    out = kit.with_suffix(".zip")
    with zipfile.ZipFile(out, "x", zipfile.ZIP_DEFLATED) as archive:
        for p in sorted(kit.rglob("*")):
            archive.write(p, Path(kit.name) / p.relative_to(kit))
    return out


def _roots(args: argparse.Namespace) -> list[tuple[str, Path, bool]]:
    """Each ``folder=path`` argument as (folder, root, real); a bad one stops the build."""
    roots = []
    for spec, real in [(s, True) for s in args.real] + [(s, False) for s in args.made_up]:
        folder, _, path = spec.partition("=")
        root = Path(path).expanduser()
        if not FOLDER.fullmatch(folder) or not (root / "outputs").is_dir():
            sys.exit(f"{spec}: want <folder>=<ctx root>, folder in a-z 0-9 _, root with outputs/")
        roots.append((folder, root, real))
    return roots


def build(out: Path, roots: list[tuple[str, Path, bool]]) -> list[Path]:
    """Copy each root's files and write the home page; return what was copied."""
    chosen, problems = {}, []
    for folder, root, _ in roots:
        chosen[folder], found = select(root)
        problems += found
    if problems:
        sys.exit("kit refused, nothing written:\n" + "\n".join(problems))
    out.mkdir(parents=True)
    shipped = []
    for folder, root, _ in roots:
        base = (root / "outputs").resolve()
        for p in chosen[folder]:
            dest = out / folder / "outputs" / p.relative_to(base)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shipped.append(Path(shutil.copy2(p, dest)))
    write_home(out, [(out / folder, real) for folder, _, real in roots])
    for name, (body, mode) in START.items():
        (out / name).write_bytes(body.encode())
        (out / name).chmod(mode)
    shutil.copy2(check_kit.__file__, out / "check_kit.py")
    shutil.copy2(CHECK_PS1, out / "check_kit.ps1")
    return shipped


def main() -> None:
    """Build the kit, check it, write its README and manifest, and pack it if asked."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, required=True, help="the kit folder; must not exist")
    parser.add_argument("--real", action="append", default=[], metavar="FOLDER=ROOT")
    parser.add_argument("--made-up", action="append", default=[], metavar="FOLDER=ROOT")
    parser.add_argument("--zip", action="store_true", help="also write <out>.zip")
    args = parser.parse_args()
    out = args.out.expanduser()
    partial = out.with_name(out.name + ".partial")
    if taken := [p for p in (out, partial, out.with_suffix(".zip")) if p.exists()]:
        sys.exit(f"{', '.join(map(str, taken))} exists; a kit is never written over")
    roots = _roots(args)
    shipped = build(partial, roots)
    names = [p.relative_to(partial).as_posix() for p in partial.rglob("*") if p.is_file()]
    pages = sorted(partial.rglob("*.html"))
    if problems := check_kit.link_problems(partial, pages) + check_kit.windows_problems(
        names, out.name
    ):
        sys.exit(f"kit checks fail; left in {partial} for inspection:\n" + "\n".join(problems))
    commit, folders = _commit(), [(f, real) for f, _, real in roots]
    readme(partial, folders, shipped, commit)
    record = json.loads(manifest(partial, folders, commit).read_text())
    partial.rename(out)
    print(f"{out}: {record['files']:,} files, {record['bytes'] / 1e6:,.1f} MB, commit {commit}")
    if args.zip:
        packed = pack(out)
        print(f"{packed}: {packed.stat().st_size / 1e6:,.1f} MB")


if __name__ == "__main__":
    main()
