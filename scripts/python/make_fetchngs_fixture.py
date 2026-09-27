#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side builder of tests/fixtures/fetchngs; needs Docker and the network
# ///
r"""Build tests/fixtures/fetchngs/<mode>: a real nf-core/fetchngs run, copied, never hand-made.

Runs the pinned fetchngs release the way the course does (`--input <outdir>/ids.csv --outdir
<outdir> --nf_core_pipeline rnaseq`) in a scratch directory, copies the outdir without its
HTML reports, deletes the scratch run so no path written into the outdir stays reachable, and
writes MANIFEST.sha256 and a README. Standard library only: it runs in the `fetch` pixi
environment, which does not hold ctx.

Usage:
    pixi run -e fetch python scripts/python/make_fetchngs_fixture.py --mode download \\
        --out tests/fixtures/fetchngs/download --scratch <dir outside the repo>
"""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROFILE = "docker,arm" if platform.machine() in {"arm64", "aarch64"} else "docker"
# nf-core/fetchngs's own test list at the commit its 1.12.0 test profile names.
NF_CORE_LIST_URL = (
    "https://raw.githubusercontent.com/nf-core/test-datasets/"
    "2732b911c57e607fa7aea5ba0c3d91b25bafb662/testdata/v1.12.0/sra_ids_test.csv"
)
NF_CORE_LIST_SHA256 = "bf6f1b2fc582c63f55914f27856f301e6d08c9f039cad5d326b9b810ce2c5fbf"
# Two FTP runs (md5 files) and two with no ENA FASTQ (automatic sra-tools fallback, no md5),
# one of each layout; all tiny, none with a GEO accession (ENA filereport, 2026-09-26).
DOWNLOAD_IDS = ("SRR13191702", "ERR1160846", "SRR14709033", "SRR14593545")
MODES = {
    "download": {
        "skip": False,
        "can_test": "verified and no_md5 files, paired and single-end, "
        "FTP and the sra-tools fallback, sheet paths that point nowhere after the move",
    },
    "metadata_only": {
        "skip": True,
        "can_test": "a metadata-only outdir (not_downloaded), a "
        "multi-run experiment (SRX9504942, 4 runs), GEO and DDBJ input IDs",
    },
}
CANNOT = (
    "runs with more than two ENA files, Aspera, a -resume into the same outdir, any "
    "cohort's sample titles, and any fetchngs release but the pinned one"
)


def fetchngs_version() -> str:
    """The pinned release, read from its one definition in ctx.gather.fetchngs."""
    source = (ROOT / "src" / "ctx" / "gather" / "fetchngs.py").read_text()
    return re.search(r'^FETCHNGS_VERSION = "([^"]+)"$', source, re.MULTILINE).group(1)


def sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def input_list(mode: str) -> bytes:
    """The ID list for a mode; the nf-core list is checked against its pinned hash."""
    if mode == "download":
        return "".join(f"{i}\n" for i in DOWNLOAD_IDS).encode()
    with urllib.request.urlopen(NF_CORE_LIST_URL, timeout=60) as response:
        data = response.read()
    if hashlib.sha256(data).hexdigest() != NF_CORE_LIST_SHA256:
        raise ValueError(f"{NF_CORE_LIST_URL} changed: sha256 differs from the pinned one")
    return data


def run(mode: str, scratch: Path) -> list[str]:
    """Run fetchngs into scratch/run_<mode>/out and return the command as run."""
    work = scratch / f"run_{mode}"
    if work.exists():
        raise FileExistsError(f"{work} exists; fetchngs must run into a fresh outdir")
    (work / "out").mkdir(parents=True)
    (work / "out" / "ids.csv").write_bytes(input_list(mode))
    command = [
        "nextflow", "run", "nf-core/fetchngs", "-r", fetchngs_version(), "-profile", PROFILE,
        "--input", "out/ids.csv", "--outdir", "out", "--nf_core_pipeline", "rnaseq",
        "-work-dir", str(work / "work"),
        # fasterq-dump asks for 36 GB by default; a laptop has less, and these runs are tiny.
        "--max_cpus", "4", "--max_memory", "8.GB",
        # nf-validation 1.1.3 reads the undefined params.monochromeLogs when colours are on,
        # and Nextflow warns; plain logs avoid it and read cleanly in a file (D-006).
        "--monochrome_logs", "true",
    ]  # fmt: skip
    if MODES[mode]["skip"]:
        command.append("--skip_fastq_download")
    env = {**os.environ, "NXF_SYNTAX_PARSER": "v1", "NXF_ANSI_LOG": "false",
           "NXF_HOME": str(scratch / "nxf_home")}  # fmt: skip
    subprocess.run(command, cwd=work, env=env, check=True)
    return command


def copy_outdir(source: Path, out: Path) -> list[tuple[str, int, str]]:
    """Copy the outdir without its HTML reports; return what was left out."""
    omitted = []
    for path in sorted(p for p in source.rglob("*") if p.is_file()):
        relative = path.relative_to(source)
        if relative.parts[0] == "pipeline_info" and path.suffix == ".html":
            omitted.append((relative.as_posix(), path.stat().st_size, sha256(path)))
            continue
        target = out / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    return omitted


def versions(out: Path) -> dict[str, str]:
    """Workflow versions from the software-versions file fetchngs writes."""
    text = (out / "pipeline_info" / "nf_core_fetchngs_software_mqc_versions.yml").read_text()
    return dict(re.findall(r"^\s+(nf-core/fetchngs|Nextflow): (\S+)$", text, re.MULTILINE))


def write_readme(out: Path, mode: str, command: list[str], omitted: list) -> None:
    """Provenance and scope of the fixture, then the manifest over every file."""
    found = versions(out)
    shown = " ".join(c if not c.startswith("/") else "<scratch>/" + Path(c).name for c in command)
    source = "the four IDs in DOWNLOAD_IDS" if mode == "download" else NF_CORE_LIST_URL
    lines = [
        f"# nf-core/fetchngs fixture: {mode}",
        "",
        "A real fetchngs run copied as written, never edited. Built by "
        "`scripts/python/make_fetchngs_fixture.py`.",
        "",
        "| Field | Value |",
        "| --- | --- |",
        f"| Command | `NXF_SYNTAX_PARSER=v1 {shown}` (run in a scratch directory) |",
        f"| nf-core/fetchngs | {found.get('nf-core/fetchngs', 'not found')} |",
        f"| Nextflow | {found.get('Nextflow', 'not found')} |",
        f"| Profile | {PROFILE} |",
        f"| Host | {platform.system()} {platform.release()} {platform.machine()} |",
        f"| Date (UTC) | {datetime.now(UTC).date().isoformat()} |",
        f"| Input list | `ids.csv`, from {source}, sha256 `{sha256(out / 'ids.csv')}` |",
        "",
        "## Left out",
        "",
        "| File | Bytes | SHA-256 |",
        "| --- | --- | --- |",
        *(f"| `{name}` | {size} | `{digest}` |" for name, size, digest in omitted),
        "",
        "## Scope",
        "",
        f"Can test: {MODES[mode]['can_test']}.",
        "",
        f"Cannot test: {CANNOT}.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines))
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    manifest = "".join(f"{sha256(p)}  {p.relative_to(out).as_posix()}\n" for p in files)
    (out / "MANIFEST.sha256").write_text(manifest)


def main(argv: list[str] | None = None) -> int:
    """Build one fixture mode."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mode", choices=sorted(MODES), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True, help="outside the repository")
    args = parser.parse_args(argv)
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError(f"{args.out} is not empty; a fixture is rebuilt into a fresh dir")
    command = run(args.mode, args.scratch)
    omitted = copy_outdir(args.scratch / f"run_{args.mode}" / "out", args.out)
    # The scratch run goes: every path fetchngs wrote into the sheet now points nowhere.
    shutil.rmtree(args.scratch / f"run_{args.mode}")
    write_readme(args.out, args.mode, command, omitted)
    print(f"{args.mode}: {sum(1 for p in args.out.rglob('*') if p.is_file())} files in {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
