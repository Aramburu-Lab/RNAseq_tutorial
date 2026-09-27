#!/usr/bin/env python3
# /// codebase-meta
# exempt: a guided walk through the modules built so far (M0, M1); runs ctx and fetchngs
# ///
r"""Walk through the tutorial as a student would, with the modules built so far: M0 and M1.

Each step prints what it teaches, runs the real command, and times it. Everything is written
under --root, a course folder outside the repository. The study accession is an argument, so
no cohort is fixed in code (the tier-1 dataset is decision D5).

Usage:
    pixi run python scripts/python/tutorial_test_run.py --study <ERP/PRJEB accession> \\
        --root ~/ctx_tutorial_test_run [--download-runs 1] [--open]
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import platform
import subprocess
import sys
import time
import urllib.request
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CTX = Path(sys.executable).parent / "ctx"
PROFILE = "docker,arm" if platform.machine() in {"arm64", "aarch64"} else "docker"
ENA = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run"
    "&fields=run_accession,library_layout,library_selection,read_count,base_count,fastq_bytes"
    "&format=tsv"
)


def fetchngs_version() -> str:
    """The pinned release, from its one definition in ctx.gather.fetchngs."""
    from ctx.gather.fetchngs import FETCHNGS_VERSION  # noqa: PLC0415 - only needs the pin

    return FETCHNGS_VERSION


def step(number: str, title: str, teaches: str) -> float:
    """Print a step banner and return its start time."""
    print(f"\n{'=' * 78}\nStep {number}: {title}\n{'-' * 78}\n{teaches}\n", flush=True)
    return time.perf_counter()


def done(started: float, *paths: Path) -> None:
    """Report a step's wall time and what to open."""
    print(f"\n  done in {time.perf_counter() - started:.1f} s", flush=True)
    for path in paths:
        print(f"  open: {path}")


def run(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    """Show a command as the student would type it, then run it."""
    print("  $ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, env={**os.environ, **(env or {})}, check=True)


def find_accession(study: str) -> list[dict[str, str]]:
    """M1a: read the study's run table from ENA and budget the transfer."""
    with urllib.request.urlopen(ENA.format(study), timeout=120) as response:
        runs = list(csv.DictReader(io.StringIO(response.read().decode()), delimiter="\t"))
    total = sum(int(b) for r in runs for b in r["fastq_bytes"].split(";") if b)
    lengths = Counter(
        round(int(r["base_count"]) / int(r["read_count"])) for r in runs if r["read_count"]
    )
    print(f"  runs: {len(runs)}")
    print(f"  library layout: {dict(Counter(r['library_layout'] for r in runs))}")
    print(f"  library selection: {dict(Counter(r['library_selection'] for r in runs))}")
    print(f"  read length (bases per read): {dict(sorted(lengths.items()))}")
    print(f"  total FASTQ: {total / 1e9:.1f} GB. Budget the transfer before you download.")
    return sorted(runs, key=lambda r: sum(int(b) for b in r["fastq_bytes"].split(";") if b))


def fetch(root: Path, slug: str, ids: list[str], skip_download: bool) -> Path:
    """M1b: run nf-core/fetchngs into a fresh outdir under data/raw/."""
    outdir = root / "data" / "raw" / slug
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"{outdir} is not empty; fetchngs must run into a fresh outdir")
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "ids.csv").write_text("".join(f"{i}\n" for i in ids))
    relative = outdir.relative_to(root)
    # --monochrome_logs true: fetchngs 1.12.0's validation plugin reads an undefined
    # params.monochromeLogs when colours are on, and Nextflow warns (D-006).
    command = [
        "pixi", "run", "--manifest-path", str(REPO / "pixi.toml"), "-e", "fetch",
        "nextflow", "run", "nf-core/fetchngs", "-r", fetchngs_version(), "-profile", PROFILE,
        "--input", f"{relative}/ids.csv", "--outdir", str(relative),
        "--nf_core_pipeline", "rnaseq", "--max_cpus", "4", "--max_memory", "8.GB",
        "--monochrome_logs", "true", "-work-dir", str(root / "work" / slug),
    ]  # fmt: skip
    if skip_download:
        command.append("--skip_fastq_download")
    run(command, cwd=root, env={"NXF_SYNTAX_PARSER": "v1", "NXF_ANSI_LOG": "false"})
    return outdir


def main(argv: list[str] | None = None) -> int:
    """Run M0 and M1 end to end."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--study", required=True, help="ENA study accession (ERP or PRJEB)")
    parser.add_argument("--root", type=Path, default=Path.home() / "ctx_tutorial_test_run")
    parser.add_argument("--download-runs", type=int, default=1, help="smallest runs to download")
    parser.add_argument("--open", action="store_true", help="open the pages when done (macOS)")
    args = parser.parse_args(argv)
    root = args.root.expanduser().resolve()
    if root.is_relative_to(REPO):
        raise ValueError(f"{root} is inside the repository; pick a course folder outside it")
    root.mkdir(parents=True, exist_ok=True)
    slug = "tier1"
    pages = []

    t = step(
        "0",
        "M0 Setup: is this computer ready?",
        "ctx records the software versions, the threads and memory it may use, and how\n"
        "long a run takes. Every later report carries the same reproduction block.",
    )
    run([str(CTX), "--root", str(root), "setup"], cwd=root)
    pages += sorted((root / "outputs" / "00_setup").glob("*.html"))
    done(t, *pages[-1:])

    t = step(
        "1a",
        "M1a Find the accession",
        "A paper points to ArrayExpress or GEO; the reads live at ENA/SRA. Read the run\n"
        "table first: layout, library selection (poly-A or not), read length, size.",
    )
    runs = find_accession(args.study)
    done(t)

    t = step("1b", "M1b Fetch, part 1: metadata only",
             "nf-core/fetchngs resolves the study to its runs and writes a samplesheet, without\n"
             "downloading any read. Minutes, and a first look at what you would get.")  # fmt: skip
    meta = fetch(root, f"{slug}_meta", [args.study], skip_download=True)
    done(t, meta / "samplesheet" / "samplesheet.csv")

    downloaded = None
    if args.download_runs > 0:
        chosen = [r["run_accession"] for r in runs[: args.download_runs]]
        size = sum(
            int(b) for r in runs[: args.download_runs] for b in r["fastq_bytes"].split(";") if b
        )
        title = (
            f"M1b Fetch, part 2: download the {len(chosen)} smallest run(s), {size / 1e6:.0f} MB"
        )
        t = step("1b", title,
                 "Each FASTQ is checked against ENA's md5 as it lands. A second copy sits in\n"
                 "work/ until you delete it.")  # fmt: skip
        downloaded = fetch(root, f"{slug}_dl", chosen, skip_download=False)
        done(t, downloaded / "fastq")

    t = step(
        "1c",
        "M1c Make it usable: ctx gather",
        "gather checks the download against itself and the md5 sums, and writes the one\n"
        "sample sheet every later module reads. It never writes into the fetchngs folder.",
    )
    for outdir in filter(None, (meta, downloaded)):
        run([str(CTX), "--root", str(root), "gather", "--dataset", outdir.name,
             "--fetchngs-dir", str(outdir)], cwd=root)  # fmt: skip
    pages += sorted((root / "outputs" / "01_gather").glob("*.html"))
    done(t, *sorted((root / "outputs" / "01_gather").glob("*.html")))

    print(f"\n{'=' * 78}\nDone. Pages to open, in order:")
    for page in pages:
        print(f"  {page}")
    print("\nNot built yet: M2-M8 (QC, alignment, quantification, DE, function, RNA vs protein,")
    print("non-coding RNA), the hub page and the notebooks. Delete work/ to free disk.")
    if args.open and sys.platform == "darwin":
        subprocess.run(["open", *map(str, pages)], check=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
