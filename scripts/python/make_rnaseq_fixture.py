#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side builder of tests/fixtures/rnaseq; needs Docker and the network
# ///
r"""Build tests/fixtures/rnaseq: a real nf-core/rnaseq Salmon-only run, copied, never hand-made.

Runs the pinned release on nf-core's test profile with the cluster's keys (Salmon only, no
alignment, index built from FASTA and GTF), then copies COPIED_FILES from
ctx.gather.rnaseq with MANIFEST.sha256 and a README. Standard library only: it runs in the
`fetch` pixi environment, which does not hold ctx.

Usage:
    pixi run -e fetch python scripts/python/make_rnaseq_fixture.py \\
        --out tests/fixtures/rnaseq --scratch <dir outside the repo> [--profile test,docker,arm64]
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src" / "ctx" / "gather" / "rnaseq.py"
ARM = platform.machine() in {"arm64", "aarch64"}
PROFILE = "test,docker,arm64" if ARM else "test,docker"
# The cluster's keys, plus nulls that stop the test profile shipping a prebuilt Salmon index,
# a transcriptome and a GFP spike-in, and bbsplit off, so the index is built as on the cluster.
PARAMS = {
    "outdir": "out",
    "pseudo_aligner": "salmon",
    "skip_alignment": True,
    "save_reference": True,
    "salmon_index": None,
    "transcript_fasta": None,
    "additional_fasta": None,
    "skip_bbsplit": True,
}
NULLED = ("salmon_index", "transcript_fasta", "additional_fasta")
# The test profile allows 15 GB; fetchngs failed at 36 GB on a 24 GB Mac (D-006).
CAP = "process {\n    resourceLimits = [cpus: 4, memory: 8.GB, time: 1.h]\n}\n"
CANNOT = (
    "GENCODE and `--gencode`; single-end 75-101 nt human reads; the cluster's Nextflow "
    "26.04.4 with Apptainer; the test GTF has transcripts with empty IDs, which the cluster "
    "GTF does not"
)


def definitions() -> tuple[str, tuple[str, ...]]:
    """RNASEQ_VERSION and COPIED_FILES from their one definition, without importing ctx."""
    names = {"RNASEQ_VERSION", "VERSIONS_FILE", "READ_FILES", "COPIED_FILES"}
    tree = ast.parse(SOURCE.read_text())
    tree.body = [
        n
        for n in tree.body
        if isinstance(n, ast.Assign) and n.targets[0].id in names  # type: ignore[attr-defined]
    ]
    values: dict[str, object] = {}
    exec(compile(tree, SOURCE, "exec"), values)  # noqa: S102 - our own constants, no imports
    return values["RNASEQ_VERSION"], values["COPIED_FILES"]  # type: ignore[return-value]


def sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(scratch: Path, version: str, profile: str) -> list[str]:
    """Run the pipeline into scratch/run/out and return the command as run."""
    work = scratch / "run"
    if work.exists():
        raise FileExistsError(f"{work} exists; the pipeline must run into a fresh outdir")
    work.mkdir(parents=True)
    # JSON is valid YAML, and it writes null and true exactly.
    (work / "params.yml").write_text(json.dumps(PARAMS, indent=2) + "\n")
    (work / "cap.config").write_text(CAP)
    command = [
        "nextflow", "run", "nf-core/rnaseq", "-r", version, "-profile", profile,
        "-params-file", "params.yml", "-c", "cap.config", "-work-dir", str(work / "work"),
    ]  # fmt: skip
    env = {**os.environ, "NXF_SYNTAX_PARSER": "v2", "NXF_ANSI_LOG": "false",
           "NXF_HOME": str(scratch / "nxf_home")}  # fmt: skip
    subprocess.run(command, cwd=work, env=env, check=True)
    return command


def copy(source: Path, out: Path, patterns: tuple[str, ...]) -> list[tuple[str, int]]:
    """Copy the files matching patterns; every pattern must match. Return what was left out."""
    chosen = {p for pattern in patterns for p in source.glob(pattern) if p.is_file()}
    empty = [pattern for pattern in patterns if not any(source.glob(pattern))]
    if empty:
        raise FileNotFoundError(f"no file in {source} matches {empty}")
    omitted = []
    for path in sorted(p for p in source.rglob("*") if p.is_file()):
        relative = path.relative_to(source)
        if path not in chosen:
            omitted.append((relative.as_posix(), path.stat().st_size))
            continue
        target = out / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    return omitted


def write_readme(out: Path, command: list[str], profile: str, omitted: list) -> None:
    """Provenance and scope of the fixture, then the manifest over every file."""
    text = (out / "pipeline_info" / "nf_core_rnaseq_software_mqc_versions.yml").read_text()
    found = dict(re.findall(r"^\s+(nf-core/rnaseq|Nextflow|salmon): (\S+)$", text, re.MULTILINE))
    (params_file,) = (out / "pipeline_info").glob("params_*.json")
    recorded = json.loads(params_file.read_text())
    nulls = ", ".join(f"`{k}` = {json.dumps(recorded.get(k))}" for k in NULLED)
    lines = [
        "# nf-core/rnaseq fixture: Salmon only",
        "",
        "A real nf-core/rnaseq run on its test profile, copied as written, never edited. Built by "
        "`scripts/python/make_rnaseq_fixture.py`.",
        "",
        "| Field | Value |",
        "| --- | --- |",
        f"| Command | `NXF_SYNTAX_PARSER=v2 {' '.join(command[:-1])} <scratch>/work` |",
        f"| Params file | `{json.dumps(PARAMS)}` |",
        f"| Resource cap | `{CAP.strip().replace(chr(10), ' ')}` |",
        f"| nf-core/rnaseq | {found.get('nf-core/rnaseq', 'not found')} |",
        f"| Nextflow | {found.get('Nextflow', 'not found')} |",
        f"| Salmon | {found.get('salmon', 'not found')} |",
        f"| Profile | {profile} |",
        f"| Host | {platform.system()} {platform.release()} {platform.machine()} |",
        f"| Date (UTC) | {datetime.now(UTC).date().isoformat()} |",
        f"| Params as recorded by the run | {nulls} |",
        "",
        "## Left out",
        "",
        "| File | Bytes |",
        "| --- | --- |",
        *(f"| `{name}` | {size} |" for name, size in omitted),
        "",
        "## Scope",
        "",
        "Can test: the layout and names ctx reads, fractional length-scaled counts, a sample "
        "with two runs merged, single- and paired-end samples, auto and reverse strandedness.",
        "",
        f"Cannot test: {CANNOT}.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines))
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    manifest = "".join(f"{sha256(p)}  {p.relative_to(out).as_posix()}\n" for p in files)
    (out / "MANIFEST.sha256").write_text(manifest)


def main(argv: list[str] | None = None) -> int:
    """Run the pipeline once and build the fixture from it."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True, help="outside the repository")
    parser.add_argument("--profile", default=PROFILE, help="emulate_amd64 if an arm64 image fails")
    args = parser.parse_args(argv)
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError(f"{args.out} is not empty; a fixture is rebuilt into a fresh dir")
    if args.scratch.resolve().is_relative_to(ROOT):
        raise ValueError(f"{args.scratch} is inside the repository")
    version, patterns = definitions()
    command = run(args.scratch, version, args.profile)
    omitted = copy(args.scratch / "run" / "out", args.out, patterns)
    write_readme(args.out, command, args.profile, omitted)
    print(f"{sum(1 for p in args.out.rglob('*') if p.is_file())} files in {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
