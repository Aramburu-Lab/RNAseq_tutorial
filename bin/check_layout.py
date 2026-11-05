#!/usr/bin/env python3
# /// codebase-meta
# exempt: Template tooling (D-019) that checks the layout on disk and writes nothing.
# ///
"""Fail when the working tree has drifted from the layout AGENTS.md promises (D-019).

Checks what is on disk, untracked files included: sanctioned top-level entries, `outputs/`
holding only `NN_<step>` and `_superseded`, no older version beside a newer one. Paths in
`.labtemplate.yml:layout_exempt` are skipped. Usage: [--stop-hook] [project-root]; as a
Claude Code Stop hook it exits 2 to keep the session working, never twice in a row.
"""

import json
import re
import sys
from fnmatch import fnmatch
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # default; a path argument overrides it

DIRS = {
    "bin",
    "conf",
    "containers",
    "deliverables",
    "knowledge",
    "logs",
    "notebooks",
    "outputs",
    "scripts",
    "settings",
    "tests",
}
# nf-core's scaffold plus what a run writes; analysis projects get none of these.
PIPELINE_DIRS = {"assets", "docs", "modules", "subworkflows", "workflows", "results", "work"}
FILES = (
    "*.md",
    "LICENSE*",
    "CITATION.cff",
    "pixi.toml",
    "pixi.lock",
    "ruff.toml",
    "pyproject.toml",
    "pytest.ini",
    "conftest.py",
    "explain_codebase.html",
    "codebase_manifest.json",
    "hpc_sync.conf",
    "hpc_sync.conf.example",
)
PIPELINE_FILES = (
    "main.nf",
    "nextflow.config",
    "nextflow_schema.json",
    "modules.json",
    "nf-test.config",
    "tower.yml",
    "ro-crate-metadata.json",
)

STEP = re.compile(r"^\d+_[A-Za-z0-9][A-Za-z0-9_.-]*$")
# A version marker: _v3, a trailing _2026-10-09, or the leading date deliverables carry.
VERSION = re.compile(r"(^\d{4}-\d{2}-\d{2}_)|(_v\d+$)|(_\d{4}-\d{2}-\d{2}$)")
SUPERSEDED = "_superseded"


def read_config(root: Path) -> tuple[str, list[str]]:
    """Return the flavour and the layout_exempt list, without needing PyYAML."""
    answers = root / ".copier-answers.yml"
    flavour = "analysis"
    if answers.exists():
        hit = re.search(r"^flavour:\s*(\S+)", answers.read_text(), re.M)
        flavour = hit.group(1) if hit else flavour
    exempt, inside = [], False
    labtemplate = root / ".labtemplate.yml"
    for line in labtemplate.read_text().splitlines() if labtemplate.exists() else []:
        if line.startswith("layout_exempt:"):
            inside = True
            exempt += re.findall(r"[^\s\[\],'\"]+", line.split(":", 1)[1].split("#")[0])
        elif inside and re.match(r"^\s+-\s", line):
            exempt.append(line.split("-", 1)[1].split("#")[0].strip().strip("'\""))
        elif inside and line.strip() and not line.lstrip().startswith("#"):
            inside = False
    return flavour, [e.rstrip("/") for e in exempt]


def stem(base: str) -> str:
    """A name without its extension, with every version marker removed."""
    previous = None
    while previous != base:
        previous, base = base, VERSION.sub("", base)
    return base


def natural(name: str) -> list:
    """Sort key that puts v9 before v10."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", name)]


def versioned_twins(folder: Path, rel: str) -> list[str]:
    """Every version but the newest of each thing in one folder.

    Grouped by stem AND extension: out_x_v2.tsv beside out_x_v2.pdf is one version in
    two formats, not two versions. A group counts only if a member carries a marker,
    so out_x.tsv beside out_x.pdf is never flagged.
    """
    groups: dict[tuple[str, str], list[str]] = {}
    for entry in folder.iterdir():
        if entry.name.startswith(".") or entry.name == SUPERSEDED:
            continue
        base, _, ext = entry.name.partition(".")
        groups.setdefault((stem(base), ext), []).append(entry.name)
    found = []
    for names in groups.values():
        if len(names) > 1 and any(VERSION.search(n.partition(".")[0]) for n in names):
            found += [f"{rel}/{n}" for n in sorted(names, key=natural)[:-1]]
    return found


def top_level_findings(root: Path, flavour: str) -> list[tuple[str, str]]:
    """Rule 1: top-level entries this flavour does not sanction."""
    dirs = DIRS | (PIPELINE_DIRS if flavour == "pipeline" else set())
    files = FILES + (PIPELINE_FILES if flavour == "pipeline" else ())
    out = []
    for entry in sorted(root.iterdir()):
        name = entry.name
        if name.startswith("."):
            continue
        if entry.is_dir() and name not in dirs:
            out.append((name, "top-level directory AGENTS.md does not sanction"))
        elif entry.is_file() and not any(fnmatch(name, g) for g in files):
            out.append((name, "loose top-level file; move it under scripts/, outputs/ or logs/"))
    return out


def output_findings(root: Path) -> list[tuple[str, str]]:
    """Rules 2 and 3: outputs/ entries that are not steps, and older versions in place."""
    # The versions rule looks one level into each sanctioned folder and no deeper: what
    # a step keeps inside its own directory is that step's business, and recursing
    # re-reported every file of a directory already flagged above.
    out, scan = [], []
    outputs, deliverables = root / "outputs", root / "deliverables"
    if outputs.is_dir():
        scan.append(outputs)
        for entry in sorted(outputs.iterdir()):
            if not entry.is_dir() or entry.name == SUPERSEDED:
                continue
            if STEP.match(entry.name):
                scan.append(entry)
            else:
                out.append((f"outputs/{entry.name}", "not a step directory outputs/NN_<step>/"))
    if deliverables.is_dir():
        scan += [
            deliverables,
            *(d for d in sorted(deliverables.iterdir()) if d.is_dir() and d.name != SUPERSEDED),
        ]
    for folder in scan:
        rel = folder.relative_to(root).as_posix()
        why = f"older version; move it to {rel}/{SUPERSEDED}/"
        out += [(p, why) for p in versioned_twins(folder, rel)]
    return out


def main() -> int:
    """Print violations; exit 1, or 2 as a Stop hook."""
    stop_hook = "--stop-hook" in sys.argv
    if stop_hook and json.loads(sys.stdin.read() or "{}").get("stop_hook_active"):
        return 0
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    root = Path(args[0]).resolve() if args else ROOT
    flavour, exempt = read_config(root)
    bad = [
        (p, why)
        for p, why in top_level_findings(root, flavour) + output_findings(root)
        if not any(p == e or p.startswith(e + "/") for e in exempt)
    ]
    if not bad:
        print(f"layout: clean ({flavour} flavour, {len(exempt)} exemption(s))")
        return 0
    stream = sys.stderr if stop_hook else sys.stdout
    print(f"layout: {len(bad)} problem(s), see AGENTS.md 'Where things live'\n", file=stream)
    for path, why in bad:
        print(f"  {path}: {why}", file=stream)
    print(
        "\nMove, never delete; scratch belongs in the system temp dir. Exempting is the user's"
        "\ncall, never a way to pass. With their say-so, paste into .labtemplate.yml:",
        file=stream,
    )
    print("\nlayout_exempt:", file=stream)
    # An exemption covers everything under it, so only the top-most paths are needed.
    paths = [p for p, _ in bad]
    for path in dict.fromkeys(paths):
        if not any(path.startswith(q + "/") for q in paths):
            print(f"  - {path}", file=stream)
    return 2 if stop_hook else 1


if __name__ == "__main__":
    sys.exit(main())
