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

``--release`` builds the kit students get (K.6, D-050): it refuses a ``--made-up`` root, and a
page carrying the preview's made-up note inside a ``--real`` root, before anything is written;
it then counts the home page's Made-up links, which must be 0, and states the release in
README.md and MANIFEST.json. A development kit keeps the previews, marked Made-up.

``--note FOLDER=TEXT`` writes one root's note into README.md, on its dataset line, and into
MANIFEST.json (K.10), so no README is edited after a build: for example
``--note "kidney=M9 on the TCGA dataset shows the [clinical] genes only."``. A note passes
the prose lint; one naming no ``--real`` or ``--made-up`` folder stops the build before
anything is written.

A raw GDC count matrix never ships (Rytis 2026-10-09; kit v8's were 94 to 205 MB). The rule
reads the data, not a list of names: a gather table of the counts kind whose README names a
GDC folder is left out, with its README kept; the kit carries that folder's GDC manifest and
fetch script, and ``LEFT_OUT.json``, README.md and MANIFEST.json name each file with its bytes
and the commands that rebuild it, the gather command read from the table's own README. That
row is provenance, the command as run, and stays verbatim in the shipped README; the command
shown to a reader names a ``--dataset-file`` the codebase committed as ``conf/datasets/<file>``
(``codebase_datasets``), and README.md names the builder's ctx root and any dataset file of
the builder's own, for the reader to replace (K.9). A page
that links a left-out matrix keeps its link: the link check counts a declared file as present
(``check_kit.declared``). Every derived table ships, and a count matrix not from GDC ships.

A source file never ships whole (M7 plan, decision 3: the CPTAC freeze rests on a data use
agreement silent on redistribution, Rieckmann's table is the publisher's): a file whose SHA-256
equals any ``[[resources]]`` entry of ``conf/datasets/*.toml`` stops the build, whatever it is
named. Derived tables and figure data ship.

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
import hashlib
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

from ctx.dataset import load  # noqa: E402
from ctx.io import naming  # noqa: E402
from ctx.report import prose, text  # noqa: E402
from ctx.report.course import label, newest, write_home  # noqa: E402
from ctx.runtime import shown_path  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
# The same checks for a Windows reader without Python, reading the kit's own MANIFEST.json.
CHECK_PS1 = REPO / "scripts" / "powershell" / "check_kit.ps1"
# Reads, alignments and their indexes: a kit never carries them (Rytis 2026-09-30).
REFUSED = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".bam", ".bai", ".sam", ".cram", ".crai",
           ".sra")  # fmt: skip
# Above the largest output measured, Step-1's M4 counts CSV at 59,259,687 B written at Salmon's
# 3 decimals (Rytis 2026-10-09, raised from 50 MB).
CAP_BYTES = 64_000_000
DATED = re.compile(r"\d{4}-\d{2}-\d{2}_")
FOLDER = re.compile(r"[a-z0-9_]+")
# Code a kit is built from; a change here, even an untracked file, marks the commit dirty.
CODE = ("src", "scripts", "conf", "pixi.toml", "pixi.lock", "pyproject.toml")
# How the home page marks a made-up link (course.html.j2's pill), counted in a release kit.
MADE_UP_LINK = '<span class="pill" data-kind="made-up">'
# A gather table of this kind, from a GDC folder, is a raw count matrix (counts, miRNA counts).
MATRIX_KIND = "_counts"
# The rows ctx's table READMEs carry (ctx.io.readme_companion): the folders read, the command.
# A sheet joining several projects names each folder after its project (gdc_run, D-076).
GDC_FOLDER_ROW = re.compile(r"^\| (?:[\w-]+ )?GDC (?:\w+ )?folder \| (.+?) \|$", re.MULTILINE)
COMMAND_ROW = re.compile(r"^\| Command \| `(.+)` \|$", re.MULTILINE)
# The dataset file's SHA-256 as gather read it, from the README's Data releases rows.
DATASET_SHA_ROW = re.compile(r"^\| sha256:([^|/ ]+\.toml) \| ([0-9a-f]{64}) \|$", re.MULTILINE)
DATASETS = "conf/datasets"
FILE_ROW = re.compile(r"^\| `([^`]+)` \| (\d+) \| `([0-9a-f]{64})` \|$", re.MULTILINE)
# How a GDC folder's own fetch script is run again, into the folder it fetched.
FETCH = "OUT={folder} bash {script}"


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


def table_readme(path: Path) -> Path:
    """The README ctx writes beside a table, shared by its CSV and Parquet copies."""
    return path.with_name(path.name.split(".")[0] + "_README.md")


def _after(words: list[str], flag: str) -> str | None:
    """The value following ``flag`` in a split command, or None."""
    return next((v for f, v in zip(words, words[1:], strict=False) if f == flag), None)


def shown_command(
    command: str, recorded: dict[str, str] | None = None, codebase: dict[str, str] | None = None
) -> tuple[str, str | None]:
    """A command as a README recorded it, with ``ctx`` for its entry point and ~ for home.

    A ``--dataset-file`` whose recorded SHA-256 (``recorded``, file name to hash) is a version
    the codebase committed (``codebase``, hash to path) is shown by that repository path.

    Returns:
        The command shown, and the dataset file it names that is the builder's own, or None.
    """
    words = command.split()
    shown, own = ["ctx" if words[0].endswith(("__main__.py", "/ctx")) else words[0]], None
    for flag, word in zip(words, words[1:], strict=False):
        mapped = None
        if flag == "--dataset-file":
            mapped = (codebase or {}).get((recorded or {}).get(Path(word).name, ""))
            own = None if mapped else shown_path(word)
        shown.append(mapped or (shown_path(word) if word.startswith("/") else word))
    return " ".join(shown), own


def codebase_datasets() -> dict[str, str]:
    """SHA-256 of every version of every dataset file the codebase committed, to its path.

    Bytes a commit holds are the repository's, from whichever clone or checkout ran them; a
    reader's own file, even one saved in ``conf/datasets/`` of a clone, has bytes no commit
    holds. A path test would take that file for the codebase's, and needs the tree on disk.
    """
    # Every parent of a merge too: a version can first appear in a merge (tcga_kirc.toml at
    # 789f8ad, which ~/ctx_kidney's gather read), and default log follows one parent.
    raw = _git("log", "--full-history", "-m", "--format=", "--raw", "--no-abbrev", "HEAD", "--",
               f"{DATASETS}/*.toml")  # fmt: skip
    blobs: dict[str, str] = {}
    for line in raw.decode().splitlines():
        meta, *paths = line.split("\t")
        if line.startswith(":") and meta.split()[3].strip("0"):
            blobs.setdefault(meta.split()[3], paths[-1])
    asked = "".join(f"{b}\n" for b in blobs).encode()
    batch = subprocess.run(["git", "-C", str(REPO), "cat-file", "--batch"], input=asked,
                           capture_output=True, check=True).stdout  # fmt: skip
    found, at = {}, 0
    for path in blobs.values():  # cat-file answers in the order asked
        end = batch.index(b"\n", at)
        size = int(batch[at:end].split()[2])
        found[hashlib.sha256(batch[end + 1 : end + 1 + size]).hexdigest()] = path
        at = end + size + 2
    return found


def gdc_matrix(path: Path, codebase: dict[str, str] | None = None) -> dict[str, Any] | None:
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
    shown, own = shown_command(command[1], recorded, codebase) if command else (None, None)
    root = _after(command[1].split(), "--root") if command else None
    sizes = {name: (int(n), sha) for name, n, sha in FILE_ROW.findall(written)}
    return {"folders": folders, "command": shown, "own_dataset_file": own,
            "root": shown_path(root) if root else None,
            "bytes": sizes.get(path.name, (path.stat().st_size, ""))[0],
            "sha256": sizes.get(path.name, (0, ""))[1] or check_kit.sha256(path)}  # fmt: skip


def declared_sources(folder: Path = REPO / "conf" / "datasets") -> dict[str, str]:
    """SHA-256 of every resource a dataset file declares, to the file and resource naming it."""
    return {r.sha256: f"{path.name} {r.name}" for path in sorted(folder.glob("*.toml"))
            for r in load(path).resources}  # fmt: skip


def select(
    root: Path, left_out: tuple[str, ...] = (), slugs: frozenset[str] = frozenset(),
    sources: dict[str, str] | None = None, codebase: dict[str, str] | None = None,
) -> tuple[list[Path], dict[Path, dict[str, Any]], list[str]]:  # fmt: skip
    """A root's files for the kit, its raw GDC count matrices apart, and what stops a ship.

    A shipped page that links to a left-out dataset's file is a problem: following the link
    would put the dataset back into the kit. A matrix left out by rule needs its rebuild
    command and its GDC folder's manifest and fetch script. A file equal to a declared source
    (``declared_sources``, the default for ``sources``) is a problem by its bytes, not its name.
    """
    sources = declared_sources() if sources is None else sources
    # Resolved, so a symlinked outputs/ holds the same paths the links resolve to.
    outputs = (root / "outputs").resolve()
    out = {p for p in current(outputs) if any(belongs(p.name, s, set(slugs)) for s in left_out)}
    chosen = current(outputs) - out
    # Undated, so current() never picks it: each module folder's entry page (GS.2) ships with
    # the folder's pages, and its links are walked like theirs.
    folders = {c.parent for c in chosen}
    chosen |= {p for p in outputs.glob(f"*/{check_kit.ENTRY}") if p.parent in folders}
    found, problems = reached({p for p in chosen if p.suffix == ".html"}, outputs)
    problems += [f"{root}: {e}" for e in check_kit.entry_problems(outputs.parent)]
    problems += [f"{p}: a shipped page links to it, but its dataset is left out"
                 for p in sorted(found & out)]  # fmt: skip
    matrices = {p: m for p in chosen | found if (m := gdc_matrix(p, codebase)) is not None}
    for p, m in sorted(matrices.items()):
        problems += [f"{p}: left out as a GDC matrix, but its README records no command"
                     ] if m["command"] is None else []  # fmt: skip
        problems += [f"{p}: its GDC folder {f} holds no manifest or fetch script"
                     for f in m["folders"] if not gdc_sources(f)]  # fmt: skip
    kept = (chosen | found) - matrices.keys()
    for p in sorted(kept):
        if p.name.lower().endswith(REFUSED):
            problems.append(f"{p}: reads or alignments never ship")
        elif p.stat().st_size > CAP_BYTES:
            problems.append(f"{p}: {p.stat().st_size:,} bytes, over the {CAP_BYTES:,} cap")
        elif (named := sources.get(check_kit.sha256(p))) is not None:
            problems.append(f"{p}: the source file {named}; derived tables ship, sources never")
    return sorted(kept), matrices, problems


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
        carried, fetches = [], []
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
            "rebuild": [*fetches, m["command"]], "builder_root": m["root"],
            "builder_dataset_file": m["own_dataset_file"],
        })  # fmt: skip
    return entries


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
                     commit=commit)  # fmt: skip
    if release:
        body += "\n" + text.word("course", "readme.release") + "\n"
    if entries := left_out_files(kit):
        listed = "\n".join(left_out_text(e) for e in entries)
        body += "\n" + text.fill("course", "readme.left_out", files=listed) + "\n"
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
            sys.exit(f"{spec}: want <folder>=<ctx root>, folder in a-z 0-9 _, root with outputs/")
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
    codebase = codebase_datasets()
    for folder, root, _ in roots:
        chosen[folder], matrices[folder], found = select(root, left_out, slugs, sources, codebase)
        problems += found
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
    args = parser.parse_args()
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
