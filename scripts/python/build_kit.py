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
        --real step1=~/ctx_tier1 --made-up preview=~/ctx_preview/synthetic \
        --leave-out tier1_teach --zip

``--leave-out`` keeps a dataset's outputs out of the kit: the teaching subset is the optional
laptop exercise, built in the Step-1 root and never shipped (D-050, Rytis 2026-10-07).

Every page is rebuilt from one commit first, so the pages and ``MANIFEST.json`` name the same
commit. The steps for the Step-1 root, kept here since a session's scratch copy was lost twice
(kits v5 and v6); ``build_preview.py`` refuses an existing root, so move the last one aside:
    R=~/ctx_tier1; C=conf/datasets; T=$R/data/raw/tier1_teach
    pixi run ctx --root $R setup
    pixi run ctx --root $R gather --dataset tier1 --fetchngs-dir $R/data/raw/tier1_meta \
        --dataset-file $C/tier1.toml --rnaseq-dir $R/data/raw/tier1_rnaseq
    pixi run ctx --root $R gather --dataset tier1_teach --fetchngs-dir $T \
        --dataset-file $C/tier1_teach.toml
    pixi run ctx --root $R qc --dataset tier1 --dataset-file $C/tier1.toml \
        --rnaseq-dir $R/data/raw/tier1_rnaseq
    pixi run ctx --root $R qc --dataset tier1_teach --dataset-file $C/tier1_teach.toml --reads \
        --fetchngs-dir $T
    pixi run ctx --root $R process --dataset tier1_teach --dataset-file $C/tier1_teach.toml \
        --genome $R/data/raw/reference/GRCh38.primary_assembly.chr22.fa.gz --aligner hisat2 \
        --compare-dataset tier1 --fetchngs-dir $T
    pixi run ctx --root $R de --dataset tier1 --dataset-file $C/tier1.toml
    pixi run ctx --root $R function --dataset tier1 --dataset-file $C/tier1.toml
    pixi run ctx --root $R report --dataset tier1 --dataset-file $C/tier1.toml
    pixi run python scripts/python/build_preview.py --root ~/ctx_preview/synthetic
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

from ctx.io import naming  # noqa: E402
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


def slugs_in(root: Path) -> set[str]:
    """The datasets a root holds: one per gather page, since every dataset enters through M1."""
    pages = newest((root / "outputs" / "01_gather").resolve(), "overview")
    return {(naming.description(p.name) or "").removesuffix("_overview") for p in pages}


def belongs(name: str, slug: str, slugs: set[str]) -> bool:
    """Whether a dated output is the dataset ``slug``'s.

    The owner is the longest known slug its description starts with, so leaving out
    ``tier1`` keeps ``tier1_teach``'s files, whose descriptions start with ``tier1_`` too.
    """
    desc = naming.description(name) or ""
    owners = [s for s in slugs | {slug} if desc == s or desc.startswith(f"{s}_")]
    return max(owners, key=len, default=None) == slug


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


def select(
    root: Path, left_out: tuple[str, ...] = (), slugs: frozenset[str] = frozenset()
) -> tuple[list[Path], list[str]]:
    """A root's files for the kit, less the left-out datasets', and every reason one cannot ship.

    A shipped page that links to a left-out dataset's file is a problem: following the link
    would put the dataset back into the kit.
    """
    # Resolved, so a symlinked outputs/ holds the same paths the links resolve to.
    outputs = (root / "outputs").resolve()
    out = {p for p in current(outputs) if any(belongs(p.name, s, set(slugs)) for s in left_out)}
    chosen = current(outputs) - out
    found, problems = reached({p for p in chosen if p.suffix == ".html"}, outputs)
    problems += [f"{p}: a shipped page links to it, but its dataset is left out"
                 for p in sorted(found & out)]  # fmt: skip
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


def manifest(
    kit: Path, roots: list[tuple[str, bool]], commit: str, left_out: tuple[str, ...] = ()
) -> Path:
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
        "left_out": sorted(left_out),
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


def build(
    out: Path, roots: list[tuple[str, Path, bool]], left_out: tuple[str, ...] = ()
) -> list[Path]:
    """Copy each root's files, less the left-out datasets', and write the home page.

    Returns:
        What was copied.
    """
    slugs = frozenset().union(*(slugs_in(root) for _, root, _ in roots))
    if unknown := sorted(set(left_out) - slugs):
        sys.exit(f"--leave-out {', '.join(unknown)}: no root holds that dataset; nothing written")
    chosen, problems = {}, []
    for folder, root, _ in roots:
        chosen[folder], found = select(root, left_out, slugs)
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
    # Its load proves START.html is beside the kit; see course.html.j2.
    (out / "kit.js").write_text("window.ctxKit = true;\n", encoding="utf-8")
    shutil.copy2(check_kit.__file__, out / "check_kit.py")
    shutil.copy2(CHECK_PS1, out / "check_kit.ps1")
    return shipped


def main() -> None:
    """Build the kit, check it, write its README and manifest, and pack it if asked."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, required=True, help="the kit folder; must not exist")
    parser.add_argument("--real", action="append", default=[], metavar="FOLDER=ROOT")
    parser.add_argument("--made-up", action="append", default=[], metavar="FOLDER=ROOT")
    parser.add_argument("--leave-out", action="append", default=[], metavar="DATASET",
                        help="a dataset slug whose outputs stay out of the kit")  # fmt: skip
    parser.add_argument("--zip", action="store_true", help="also write <out>.zip")
    args = parser.parse_args()
    out = args.out.expanduser()
    partial = out.with_name(out.name + ".partial")
    if taken := [p for p in (out, partial, out.with_suffix(".zip")) if p.exists()]:
        sys.exit(f"{', '.join(map(str, taken))} exists; a kit is never written over")
    roots = _roots(args)
    left_out = tuple(args.leave_out)
    shipped = build(partial, roots, left_out)
    names = [p.relative_to(partial).as_posix() for p in partial.rglob("*") if p.is_file()]
    pages = sorted(partial.rglob("*.html"))
    if problems := check_kit.link_problems(partial, pages) + check_kit.windows_problems(
        names, out.name
    ):
        sys.exit(f"kit checks fail; left in {partial} for inspection:\n" + "\n".join(problems))
    commit, folders = _commit(), [(f, real) for f, _, real in roots]
    readme(partial, folders, shipped, commit)
    record = json.loads(manifest(partial, folders, commit, left_out).read_text())
    partial.rename(out)
    print(f"{out}: {record['files']:,} files, {record['bytes'] / 1e6:,.1f} MB, commit {commit}")
    if left_out:
        print(f"left out: {', '.join(sorted(left_out))}")
    if args.zip:
        packed = pack(out)
        print(f"{packed}: {packed.stat().st_size / 1e6:,.1f} MB")


if __name__ == "__main__":
    main()
