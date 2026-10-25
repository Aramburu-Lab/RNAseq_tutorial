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
        --real step1=~/ctx_Non-cancerous-blood --made-up preview=~/ctx_preview/synthetic --zip

``--leave-out`` keeps a dataset's outputs out of the kit, as the teaching subset was kept out
until it was removed (D-050, D-115).

``--release`` builds the kit students get (K.6, D-050): it refuses a ``--made-up`` root, and a
page carrying the preview's made-up note inside a ``--real`` root, before anything is written;
it then counts the home page's Made-up links, which must be 0, and states the release in
README.md and MANIFEST.json. A development kit keeps the previews, marked Made-up.

``--note FOLDER=TEXT`` writes one root's note into README.md, on its dataset line, and into
MANIFEST.json (K.10), so no README is edited after a build: for example
``--note "kidney=M9 on the TCGA dataset shows the [clinical] genes only."``. A note passes
the prose lint; one naming no ``--real`` or ``--made-up`` folder stops the build before
anything is written.

``--edition NAME`` (``CTX_EDITION``, else full; ``ctx.editions``) ships only that edition's
module folders, refuses a root holding a course dataset outside it, and is recorded in
README.md and MANIFEST.json; START.html lists only its modules. Every ctx page records the
edition it was built under (``<html data-edition>``, ED.2b); a page of another edition stops
the build, named with its edition, and a page recording none counts as full's.

A raw GDC count matrix never ships, with what rebuilds it: ``kit_matrices.py``.

A source file never ships whole (M7 plan, decision 3: the CPTAC freeze rests on a data use
agreement silent on redistribution, Rieckmann's table is the publisher's): a file whose SHA-256
equals any ``[[resources]]`` entry of ``conf/datasets/*.toml`` stops the build, whatever it is
named. Derived tables and figure data ship.

Every page is rebuilt from one commit first, so the pages and ``MANIFEST.json`` name the same
commit. The steps for the Step-1 root, kept here since a session's scratch copy was lost twice
(kits v5 and v6); ``build_preview.py`` refuses an existing root, so move the last one aside:
    N=Non-cancerous-blood; R=~/ctx_$N; C=conf/datasets
    pixi run ctx --root $R setup
    pixi run ctx --root $R gather --dataset $N --fetchngs-dir $R/data/raw/${N}_meta \
        --dataset-file $C/$N.toml --rnaseq-dir $R/data/raw/${N}_rnaseq
    pixi run ctx --root $R qc --dataset $N --dataset-file $C/$N.toml \
        --rnaseq-dir $R/data/raw/${N}_rnaseq
    pixi run ctx --root $R de --dataset $N --dataset-file $C/$N.toml
    pixi run ctx --root $R function --dataset $N --dataset-file $C/$N.toml
    pixi run ctx --root $R report --dataset $N --dataset-file $C/$N.toml
    pixi run python scripts/python/build_preview.py --root ~/ctx_preview/synthetic
"""

from __future__ import annotations

import argparse
import html
import json
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_kit  # noqa: E402 - the kit's own checker, beside this file
from kit_datasets import Version, codebase_datasets  # noqa: E402
from kit_matrices import carry_matrices, gdc_matrix, gdc_sources  # noqa: E402

from ctx import editions  # noqa: E402
from ctx.dataset import load  # noqa: E402
from ctx.io import naming  # noqa: E402
from ctx.report import prose, text  # noqa: E402
from ctx.report.builder import shown  # noqa: E402
from ctx.report.course import label, newest, write_home  # noqa: E402
from ctx.report.hub import PAGE  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
# The same checks for a Windows reader without Python, reading the kit's own MANIFEST.json.
CHECK_PS1 = REPO / "scripts" / "powershell" / "check_kit.ps1"
# Reads, alignments and their indexes: a kit never carries them (Rytis 2026-09-30).
REFUSED = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".bam", ".bai", ".sam", ".cram", ".crai",
           ".sra")  # fmt: skip
# Above the largest output measured, Step-1's M4 counts CSV at 59,259,687 B written at Salmon's
# 3 decimals (Rytis 2026-10-09, raised from 50 MB).
CAP_BYTES = 64_000_000
# The files a lite kit left out for their size, read back into README.md.
LARGE = "LEFT_OUT_LARGE.json"
DATED = re.compile(r"\d{4}-\d{2}-\d{2}_")
# A kit folder may carry a dataset name, so it is an output description (D-114).
FOLDER = re.compile(naming.DESCRIPTION)
# Code a kit is built from; a change here, even an untracked file, marks the commit dirty.
CODE = ("src", "scripts", "conf", "pixi.toml", "pixi.lock", "pyproject.toml")
# How the home page marks a made-up link (course.html.j2's pill), counted in a release kit.
MADE_UP_LINK = '<span class="pill" data-kind="made-up">'
# Report pages no module run writes again: `ctx guide` writes the guide, and the glossary is
# written with the first page of each day that links it (D-048).
NOT_MODULE = ("guide", "glossary")
# A page's Reproduce command (module.html.j2).
PAGE_COMMAND = re.compile(r'<code id="ctx-command">(.*?)</code>', re.S)


def current(outputs: Path, slugs: frozenset[str] = frozenset()) -> set[Path]:
    """Each module folder's newest run of each dataset, and the newest copy of any other output.

    A run writes all its files on one day, so a dataset's newest day in a folder is its run;
    keeping the newest copy per file name instead shipped 114 files of an older M9 run whose
    genes the newer run did not draw (kit v9, 2026-10-16). A file no known dataset owns keeps
    the per-name rule.
    """
    newest_day: dict[tuple[Path, str], str] = {}
    files = [p for p in outputs.glob("*/*") if p.is_file() and DATED.match(p.name)]
    for p in files:
        key = _run_key(p, slugs)
        newest_day[key] = max(newest_day.get(key, ""), p.name[:10])
    return {p for p in files if p.name[:10] == newest_day[_run_key(p, slugs)]}


def _run_key(path: Path, slugs: frozenset[str]) -> tuple[Path, str]:
    """A file's module folder and its dataset, or its undated name when no dataset owns it."""
    desc = naming.description(path.name) or ""
    owners = [s for s in slugs if desc == s or desc.startswith(f"{s}_")]
    return path.parent, max(owners, key=len) if owners else "name:" + path.name[11:]


def slugs_in(root: Path) -> set[str]:
    """The datasets a root holds: one per gather page, since every dataset enters through M1."""
    pages = newest((root / "outputs" / "01_gather").resolve(), "overview")
    return {(naming.description(p.name) or "").removesuffix("_overview") for p in pages}


def belongs(name: str, slug: str, slugs: set[str]) -> bool:
    """Whether a dated output is the dataset ``slug``'s.

    The owner is the longest known slug its description starts with, then ``_``, so leaving
    out a dataset keeps the files of another whose slug starts with its slug and ``_``.
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


def declared_sources(folder: Path = REPO / "conf" / "datasets") -> dict[str, str]:
    """SHA-256 of every resource a dataset file declares, to the file and resource naming it."""
    return {r.sha256: f"{path.name} {r.name}" for path in sorted(folder.glob("*.toml"))
            for r in load(path).resources}  # fmt: skip


def select(
    root: Path, left_out: tuple[str, ...] = (), slugs: frozenset[str] = frozenset(),
    sources: dict[str, str] | None = None, codebase: dict[str, Version] | None = None,
) -> tuple[list[Path], dict[Path, dict[str, Any]], list[str], list[Path]]:  # fmt: skip
    """A root's files, its GDC matrices apart, what stops a ship, and lite's oversized files.

    A shipped page that links to a left-out dataset's file is a problem: following the link
    would put the dataset back into the kit. A matrix left out by rule needs its rebuild
    command and its GDC folder's manifest and fetch script. A file equal to a declared source
    (``declared_sources``, the default for ``sources``) is a problem by its bytes, not its name.
    """
    sources = declared_sources() if sources is None else sources
    # Resolved, so a symlinked outputs/ holds the same paths the links resolve to.
    outputs = (root / "outputs").resolve()
    newest = current(outputs, frozenset(slugs))
    out = {p for p in newest if any(belongs(p.name, s, set(slugs)) for s in left_out)}
    # The module folders the selected edition holds: no other folder's file ships.
    edition = editions.active()
    held = {d.name for d in outputs.iterdir() if edition.holds(d.name)}
    chosen = {p for p in newest - out if p.relative_to(outputs).parts[0] in held}
    # A lite page with nothing to show stays home, as the hub and START.html leave it out.
    chosen = {p for p in chosen if p.suffix != ".html" or shown(p)}
    # Undated, so current() never picks it: each module folder's entry page (GS.2) ships with
    # the folder's pages, and its links are walked like theirs.
    folders = {c.parent for c in chosen}
    chosen |= {p for p in outputs.glob(f"*/{check_kit.ENTRY}") if p.parent in folders}
    found, problems = reached({p for p in chosen if p.suffix == ".html"}, outputs)
    # An entry page is needed only where its folder ships: "outputs/<folder>: ..." names it.
    problems += [f"{root}: {e}" for e in check_kit.entry_problems(outputs.parent)
                 if Path(e.partition(":")[0]).name in held]  # fmt: skip
    problems += [f"{p}: a shipped page links to it, but its dataset is left out"
                 for p in sorted(found & out)]  # fmt: skip
    problems += [f"{p}: a shipped page links to it, outside the {edition.name} edition"
                 for p in sorted(found) if p.relative_to(outputs).parts[0] not in held]  # fmt: skip
    matrices = {p: m for p in chosen | found if (m := gdc_matrix(p, codebase)) is not None}
    for p, m in sorted(matrices.items()):
        problems += [f"{p}: left out as a GDC matrix, but its README records no command"
                     ] if m["command"] is None else []  # fmt: skip
        problems += [f"{p}: its GDC folder {f} holds no manifest or fetch script"
                     for f in m["folders"] if not gdc_sources(f)]  # fmt: skip
    kept = (chosen | found) - matrices.keys()
    # A lite edition's pages offer no Reproduce panel, so a table over the cap, a redraw input
    # such as M3's every-gene comparison (633 MB on Breast-tissue-BRCA, 2026-10-25), stays home
    # and README.md names it; a full kit, whose panels redraw from it, still refuses it.
    large = sorted(
        p
        for p in kept
        if editions.active().page == "lite"
        and p.stat().st_size > CAP_BYTES
        and not p.name.lower().endswith(REFUSED)
    )
    kept -= set(large)
    problems += stale_pages(sorted(kept))
    problems += other_edition(sorted(kept), outputs, edition)
    for p in sorted(kept):
        if p.name.lower().endswith(REFUSED):
            problems.append(f"{p}: reads or alignments never ship")
        elif p.stat().st_size > CAP_BYTES:
            problems.append(f"{p}: {p.stat().st_size:,} bytes, over the {CAP_BYTES:,} cap")
        elif (named := sources.get(check_kit.sha256(p))) is not None:
            problems.append(f"{p}: the source file {named}; derived tables ship, sources never")
    return sorted(kept), matrices, problems, large


def other_edition(kept: list[Path], outputs: Path, edition: editions.Edition) -> list[str]:
    """Each ctx page built under another edition than the kit's, naming that edition.

    A ctx page is a dated page or an entry page in a module folder; MultiQC's report records
    none. A page built before pages recorded one (ED.2b) was built under full: a full kit takes
    it, any other refuses it, since a full hub lists M7 and M10 as not run.
    """
    pages = [p for p in kept if p.suffix == ".html" and p.parent.parent == outputs
             and (naming.check(p.name) or p.name == check_kit.ENTRY)]  # fmt: skip
    fix = f"rebuild it under {edition.name} (ctx --edition {edition.name})"
    problems = []
    for p in pages:
        built = editions.recorded(p.read_text(encoding="utf-8", errors="replace"))
        if built is None and edition.name != editions.DEFAULT:
            problems.append(f"{p}: records no edition, built before pages recorded one; {fix}")
        elif built is not None and built != edition.name:
            problems.append(f"{p}: built under the {built} edition, not {edition.name}; {fix}")
    return problems


def _not_module(path: Path) -> bool:
    """Whether a file is a guide or glossary page, which no module run writes again."""
    m = PAGE.fullmatch(path.name)
    return bool(m) and m["module"] == "report" and m["what"] in NOT_MODULE


def stale_pages(kept: list[Path]) -> list[str]:
    """Each guide or glossary page no module page links, older than the root's newest module page.

    Kit v10 shipped ~/ctx_tcga_hnsc's 2026-10-18 guide, and the glossary it links, beside that
    root's 2026-10-20 pages (K.15). Per root, because a module run again on a later day leaves
    the root's other pages linking their own day's glossary, which a module page then links.
    """
    pages = [p for p in kept if p.suffix == ".html" and PAGE.fullmatch(p.name)]
    modules = [p for p in pages if not _not_module(p)]
    latest = max((p.name for p in modules), key=lambda name: name[:10], default="")
    old = [p for p in pages if _not_module(p) and p.name[:10] < latest[:10]]
    if not old:
        return []
    links = {p: {check_kit.target(p, v) for _, _, v in check_kit.links(p)} for p in pages}
    problems = []
    for p in old:
        if any(p in links[m] for m in modules):
            continue
        command = PAGE_COMMAND.search(p.read_text(encoding="utf-8", errors="replace"))
        linkers = sorted(q.name for q in pages if p in links[q])
        fix = (f"write it again with `{html.unescape(command[1])}`" if command else
               f"{', '.join(linkers) or 'no page'} links it; write that page again")  # fmt: skip
        problems.append(f"{p}: from {p.name[:10]}, older than {latest}, and no module page "
                        f"links it; {fix}")  # fmt: skip
    return problems


def _git(*args: str) -> bytes:
    """What git prints for one command in this repository."""
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True,
                          check=True).stdout  # fmt: skip


def _commit() -> str:
    """The codebase commit that built the kit, marked dirty when its code had any change.

    Untracked code counts: a builder not yet committed cannot be rebuilt from the commit.
    Untracked notes elsewhere in the tree do not build the kit, so they do not count.
    """
    changed = _git("status", "--porcelain", "--untracked-files=all", "--", *CODE).strip()
    return _git("rev-parse", "--short", "HEAD").decode().strip() + ("-dirty" if changed else "")


def made_up_pages(files: list[Path]) -> list[str]:
    """Every page that says it is drawn on made-up data: it carries the preview's note."""
    note = text.word("ui", "preview.note")
    pages = [p for p in files if p.suffix == ".html"]
    return [f"{p}: a made-up page; a release kit holds none" for p in pages
            if note in p.read_text(encoding="utf-8", errors="replace")]  # fmt: skip


def release_problems(roots: list[tuple[str, Path, bool]], chosen: dict[str, list[Path]]
                     ) -> list[str]:  # fmt: skip
    """Why these roots cannot make a release kit: a made-up root, or a made-up page."""
    problems = [f"--made-up {folder}={root}: a release kit holds no made-up root"
                for folder, root, real in roots if not real]  # fmt: skip
    return problems + made_up_pages([p for files in chosen.values() for p in files])


def left_out_files(kit: Path) -> list[dict[str, Any]]:
    """The kit's LEFT_OUT.json entries; none where nothing was left out by rule."""
    record = kit / check_kit.LEFT_OUT
    return json.loads(record.read_text(encoding="utf-8"))["files"] if record.is_file() else []


def made_up_links(kit: Path) -> int:
    """Made-up links on the kit's home page."""
    return (kit / "START.html").read_text(encoding="utf-8").count(MADE_UP_LINK)


def manifest(
    kit: Path,
    roots: list[tuple[str, bool]],
    commit: str,
    left_out: tuple[str, ...] = (),
    release: bool = False,
    notes: dict[str, str] | None = None,
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
        "roots": [
            {"folder": f, "real": real, "note": (notes or {}).get(f)} for f, real in roots
        ],  # fmt: skip
        "left_out": sorted(left_out),
        "release": release,
        "edition": editions.active().name,
        "made_up_links": made_up_links(kit),
        "left_out_files": left_out_files(kit),
        "files": len(files),
        "bytes": sum(f["bytes"] for f in files),
        "entries": files,
    }
    out.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    return out


def left_out_text(entry: dict[str, Any]) -> str:
    """One left-out file in README.md: its commands, and which of their paths are the builder's."""
    said = [text.fill("course", "readme.left_out_file", path=entry["path"],
                      bytes=f"{entry['bytes']:,}",
                      commands="\n".join(f"    {c}" for c in entry["rebuild"]))]  # fmt: skip
    if entry.get("builder_root"):
        said.append(text.fill("course", "readme.left_out_root", root=entry["builder_root"]))
    if entry.get("builder_dataset_file"):
        said.append(text.fill("course", "readme.left_out_own",
                              path=entry["builder_dataset_file"]))  # fmt: skip
    if entry.get("dataset_file_commit"):
        path = entry["dataset_file_path"]
        said.append(text.fill("course", "readme.left_out_earlier", name=Path(path).name,
                              path=path, commit=entry["dataset_file_commit"]))  # fmt: skip
    return "\n".join(said)


def readme(kit: Path, roots: list[tuple[str, bool]], shipped: list[Path], commit: str,
           release: bool = False, notes: dict[str, str] | None = None) -> Path:  # fmt: skip
    """The kit's README.md, its text from ``content/course.toml``; a root's note on its line."""
    lines = []
    for folder, real in roots:
        pages = [p for d in sorted((kit / folder / "outputs").iterdir()) if d.is_dir()
                 for p in newest(d, "overview")]  # fmt: skip
        names = sorted({label(p) for p in pages})
        kind = text.word("course", "real" if real else "made_up").lower()
        lines.append(text.fill("course", "readme.dataset", folder=folder, names=", ".join(names),
                               kind=kind, note=(notes or {}).get(folder, "")).rstrip())  # fmt: skip
    size = f"{sum(p.stat().st_size for p in shipped) / 1e6:,.0f} MB"
    body = text.fill("course", "readme.text", files=f"{len(shipped):,}", size=size,
                     datasets="\n".join(lines), built=date.today().isoformat(),
                     commit=commit, edition=edition_text())  # fmt: skip
    if release:
        body += "\n" + text.word("course", "readme.release") + "\n"
    if entries := left_out_files(kit):
        listed = "\n".join(left_out_text(e) for e in entries)
        body += "\n" + text.fill("course", "readme.left_out", files=listed) + "\n"
    if (kit / LARGE).is_file():
        big = json.loads((kit / LARGE).read_text(encoding="utf-8"))["files"]
        listed = "\n".join(text.fill("course", "readme.too_large_file", path=e["path"],
                                      bytes=f"{e['bytes']:,}") for e in big)  # fmt: skip
        body += "\n" + text.fill("course", "readme.too_large", cap=f"{CAP_BYTES / 1e6:g}",
                                  files=listed) + "\n"  # fmt: skip
    out = kit / "README.md"
    out.write_text(body, encoding="utf-8")
    return out


def edition_text() -> str:
    """The README's line on the selected edition: its name, label and modules."""
    edition = editions.active()
    return text.fill("course", "readme.edition", label=edition.label, name=edition.name,
                     modules=", ".join(edition.modules))  # fmt: skip


def pack(kit: Path) -> Path:
    """One zip holding the kit's folder; file modes are kept, so the scripts stay runnable."""
    out = kit.with_suffix(".zip")
    with zipfile.ZipFile(out, "x", zipfile.ZIP_DEFLATED) as archive:
        for p in sorted(kit.rglob("*")):
            archive.write(p, Path(kit.name) / p.relative_to(kit))
    return out


def _notes(args: argparse.Namespace, roots: list[tuple[str, Path, bool]]) -> dict[str, str]:
    """Each ``--note folder=text`` as {folder: text}; a bad one stops the build, writing nothing.

    A note goes into the README beside the course's own text, so it passes the same prose lint.
    """
    notes: dict[str, str] = {}
    folders, problems = {f for f, _, _ in roots}, []
    for spec in args.note:
        folder, _, note = spec.partition("=")
        where = f"--note {folder}"
        if folder not in folders:
            problems.append(f"{where}: no --real or --made-up folder of that name")
        if folder in notes or not note.strip():
            problems.append(f"{where}: one note per folder, and not empty")
        problems += [*prose.check_text(where, note, None), *prose.asks(where, note)]
        notes[folder] = note.strip()
    if problems:
        sys.exit("kit refused, nothing written:\n" + "\n".join(problems))
    return notes


def _roots(args: argparse.Namespace) -> list[tuple[str, Path, bool]]:
    """Each ``folder=path`` argument as (folder, root, real); a bad one stops the build."""
    roots = []
    for spec, real in [(s, True) for s in args.real] + [(s, False) for s in args.made_up]:
        folder, _, path = spec.partition("=")
        root = Path(path).expanduser()
        if not FOLDER.fullmatch(folder) or not (root / "outputs").is_dir():
            sys.exit(f"{spec}: want <folder>=<ctx root>, folder of letters, digits, - and _, "
                     "root with outputs/")  # fmt: skip
        roots.append((folder, root, real))
    return roots


def build(
    out: Path,
    roots: list[tuple[str, Path, bool]],
    left_out: tuple[str, ...] = (),
    release: bool = False,
) -> list[Path]:
    """Copy each root's files, less the left-out datasets', and write the home page.

    A release build refuses every made-up root and page first, and writes nothing then.

    Returns:
        What was copied.
    """
    slugs = frozenset().union(*(slugs_in(root) for _, root, _ in roots))
    if unknown := sorted(set(left_out) - slugs):
        sys.exit(f"--leave-out {', '.join(unknown)}: no root holds that dataset; nothing written")
    chosen, matrices, problems, sources = {}, {}, [], declared_sources()
    large: list[dict[str, Any]] = []
    problems += [f"{folder}={root}: {refused}" for folder, root, _ in roots
                 for s in sorted(slugs_in(root) - set(left_out))
                 if (refused := editions.active().refusal(dataset=s))]  # fmt: skip
    codebase = codebase_datasets()
    for folder, root, _ in roots:
        chosen[folder], matrices[folder], found, big = select(root, left_out, slugs, sources,
                                                              codebase)  # fmt: skip
        problems += found
        base = (root / "outputs").resolve()
        large += [{"path": f"{folder}/outputs/{p.relative_to(base).as_posix()}",
                   "bytes": p.stat().st_size} for p in big]  # fmt: skip
    if release:
        problems += release_problems(roots, chosen)
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
    entries = [e for folder, root, _ in roots
               for e in carry_matrices(out, folder, (root / "outputs").resolve(),
                                       matrices[folder])]  # fmt: skip
    if entries:
        (out / check_kit.LEFT_OUT).write_text(json.dumps({"files": entries}, indent=1) + "\n",
                                              encoding="utf-8")  # fmt: skip
    if large:
        (out / LARGE).write_text(json.dumps({"files": large}, indent=1) + "\n", encoding="utf-8")
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
    parser.add_argument("--note", action="append", default=[], metavar="FOLDER=TEXT",
                        help="a sentence on one root, for README.md and MANIFEST.json")  # fmt: skip
    parser.add_argument("--zip", action="store_true", help="also write <out>.zip")
    parser.add_argument("--release", action="store_true",
                        help="the students' kit: no made-up root or page")  # fmt: skip
    parser.add_argument("--edition", help="conf/editions/<name>.toml; CTX_EDITION, else full")
    args = parser.parse_args()
    try:
        editions.choose(args.edition)
        editions.active()  # CTX_EDITION too, checked before anything is written
    except editions.EditionError as error:
        sys.exit(f"kit refused, nothing written: {error}")
    out = args.out.expanduser()
    partial = out.with_name(out.name + ".partial")
    if taken := [p for p in (out, partial, out.with_suffix(".zip")) if p.exists()]:
        sys.exit(f"{', '.join(map(str, taken))} exists; a kit is never written over")
    roots = _roots(args)
    notes = _notes(args, roots)
    left_out = tuple(args.leave_out)
    shipped = build(partial, roots, left_out, args.release)
    names = [p.relative_to(partial).as_posix() for p in partial.rglob("*") if p.is_file()]
    pages = sorted(partial.rglob("*.html"))
    problems = check_kit.link_problems(partial, pages) + check_kit.windows_problems(names, out.name)
    if args.release and (links := made_up_links(partial)):
        problems.append(f"START.html: {links} Made-up links; a release kit has 0")
    if problems:
        sys.exit(f"kit checks fail; left in {partial} for inspection:\n" + "\n".join(problems))
    commit, folders = _commit(), [(f, real) for f, _, real in roots]
    readme(partial, folders, shipped, commit, args.release, notes)
    record = json.loads(manifest(partial, folders, commit, left_out, args.release,
                                 notes).read_text())  # fmt: skip
    partial.rename(out)
    print(f"{out}: {record['files']:,} files, {record['bytes'] / 1e6:,.1f} MB, commit {commit}")
    if left_out:
        print(f"left out: {', '.join(sorted(left_out))}")
    if args.release:
        print(f"release kit: {record['made_up_links']} Made-up links on START.html")
    if args.zip:
        packed = pack(out)
        print(f"{packed}: {packed.stat().st_size / 1e6:,.1f} MB")


if __name__ == "__main__":
    main()
