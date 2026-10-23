#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side builder of tests/fixtures/rnaseq and rnaseq_star; Docker, network
# ///
r"""Build tests/fixtures/rnaseq or rnaseq_star: a real nf-core/rnaseq run, copied, never edited.

Runs the pinned release on nf-core's test profile with the cluster's keys (Salmon only, no
alignment, index built from FASTA and GTF), then copies COPIED_FILES from
ctx.gather.rnaseq with MANIFEST.sha256 and a README. Standard library only: it runs in the
`fetch` pixi environment, which does not hold ctx.

Usage:
    pixi run -e fetch python scripts/python/make_rnaseq_fixture.py \\
        --out tests/fixtures/rnaseq --scratch <dir outside the repo> [--profile test,docker,arm64]

With ``--from-run`` it copies from the run already in ``<scratch>/run`` instead of running
again, as when COPIED_FILES grows and the run that made the fixture is still at hand.

With ``--aligned`` it builds tests/fixtures/rnaseq_star: the test profile with
``aligner: star_salmon`` beside ``pseudo_aligner: salmon`` and Preseq on, as the cluster aligns
(D-113, NF.7), and copies COPIED_FILES and ALIGNED_FILES, never a BAM. ``--results <outdir>``
copies from an outdir already run that way instead of running again:

    pixi run -e fetch python scripts/python/make_rnaseq_fixture.py --aligned \\
        --out tests/fixtures/rnaseq_star --scratch <dir> [--results <outdir of such a run>]
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
    # Stated, not left to the pipeline: under NXF_SYNTAX_PARSER=v2 the -c override below cannot
    # see a param only the pipeline sets, and its publishDir falls back to symlinks (D-025).
    "publish_dir_mode": "copy",
}
NULLED = ("salmon_index", "transcript_fasta", "additional_fasta")
# The aligned fixture's keys: the test profile's own references and indices, as the run of
# 2026-10-23 that made it was given them; only the route and the publish mode are set.
ALIGNED_PARAMS = {"outdir": "out", "aligner": "star_salmon", "pseudo_aligner": "salmon",
                  "publish_dir_mode": "copy", "skip_preseq": False}  # fmt: skip
# The tools whose versions the aligned README records, as the versions file names them.
ALIGNED_TOOLS = ("nf-core/rnaseq", "Nextflow", "star", "salmon", "samtools", "rseqc",
                 "qualimap", "preseq", "bioconductor-dupradar", "subread", "trimgalore",
                 "fastqc")  # fmt: skip
# The test profile allows 15 GB; fetchngs failed at 36 GB on a 24 GB Mac (D-006).
CAP = "process {\n    resourceLimits = [cpus: 4, memory: 8.GB, time: 1.h]\n}\n"
# nf-core/rnaseq 3.27.0 publishes only TRIMGALORE's *.txt reports
# (conf/modules/qc_trim_filter.config), and MultiQC 1.35 reads Trim Galore 2.x only from
# *_trimming_report.json (D-025). The same three entries, the last widened to the JSON report;
# the cluster block in TODO.md carries this text, held equal by
# tests/test_repo.py::test_publish_override_single_source. One Groovy line is split in two
# literals only to keep this file within the line length.
PUBLISH = (
    """process {
    withName: '.*:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE' {
        publishDir = [
            [
                path: { "${params.outdir}/${task.ext.publish_prefix}fastqc/trim" },
                mode: params.publish_dir_mode,
                pattern: "*.{html,zip}"
            ],
            [
                path: { params.save_trimmed ? "${params.outdir}/${task.ext.publish_prefix}"""
    """${params.trimmer}" : params.outdir },
                mode: params.publish_dir_mode,
                pattern: "*.fq.gz",
                saveAs: { filename -> params.save_trimmed ? filename : null }
            ],
            [
                path: { "${params.outdir}/${task.ext.publish_prefix}${params.trimmer}" },
                mode: params.publish_dir_mode,
                pattern: "*_trimming_report.{txt,json}"
            ]
        ]
    }
}
"""
)
CANNOT = (
    "GENCODE and `--gencode`; single-end 75-101 nt human reads; the cluster's Nextflow "
    "26.04.4 with Apptainer; the test GTF has transcripts with empty IDs, which the cluster "
    "GTF does not"
)


def definitions(aligned: bool = False) -> tuple[str, tuple[str, ...]]:
    """RNASEQ_VERSION and the files to copy, from their one definition, without importing ctx.

    The files are COPIED_FILES, and ALIGNED_FILES after them with ``aligned``.
    """
    names = {"RNASEQ_VERSION", "VERSIONS_FILE", "MATRIX_FILE", "READ_FILES", "COPIED_FILES",
             "ALIGNED_MATRIX", "MULTIQC_DATA", "MULTIQC_TABLES", "ALIGNED_FILES"}  # fmt: skip
    tree = ast.parse(SOURCE.read_text())
    tree.body = [
        n
        for n in tree.body
        if isinstance(n, ast.Assign) and n.targets[0].id in names  # type: ignore[attr-defined]
    ]
    values: dict[str, object] = {}
    exec(compile(tree, SOURCE, "exec"), values)  # noqa: S102 - our own constants, no imports
    extra = values["ALIGNED_FILES"] if aligned else ()
    files = (*values["COPIED_FILES"], *extra)  # type: ignore[misc]
    return values["RNASEQ_VERSION"], files  # type: ignore[return-value]


def sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command_for(work: Path, version: str, profile: str) -> list[str]:
    """The pipeline command, as run() runs it."""
    return [
        "nextflow", "run", "nf-core/rnaseq", "-r", version, "-profile", profile,
        "-params-file", "params.yml", "-c", "cap.config", "-c", "publish.config",
        "-work-dir", str(work / "work"),
    ]  # fmt: skip


def aligned_command(work: Path, version: str, profile: str) -> list[str]:
    """The aligned fixture's command: the test profile's own resources, no cap."""
    return [
        "nextflow", "run", "nf-core/rnaseq", "-r", version, "-profile", profile,
        "-params-file", "params.yml", "-c", "publish.config", "-work-dir", str(work / "work"),
    ]  # fmt: skip


def run(scratch: Path, version: str, profile: str, aligned: bool = False) -> list[str]:
    """Run the pipeline into scratch/run/out and return the command as run."""
    work = scratch / "run"
    if work.exists():
        raise FileExistsError(f"{work} exists; the pipeline must run into a fresh outdir")
    work.mkdir(parents=True)
    # JSON is valid YAML, and it writes null and true exactly.
    params = ALIGNED_PARAMS if aligned else PARAMS
    (work / "params.yml").write_text(json.dumps(params, indent=2) + "\n")
    (work / "cap.config").write_text(CAP)
    (work / "publish.config").write_text(PUBLISH)
    command = (aligned_command if aligned else command_for)(work, version, profile)
    env = {**os.environ, "NXF_SYNTAX_PARSER": "v2", "NXF_ANSI_LOG": "false",
           "NXF_HOME": str(scratch / "nxf_home")}  # fmt: skip
    subprocess.run(command, cwd=work, env=env, check=True)
    return command


def copy(source: Path, out: Path, patterns: tuple[str, ...]) -> list[tuple[str, int]]:
    """Copy the files matching patterns; every pattern must match. Return what was left out."""
    chosen = {p for pattern in patterns for p in source.glob(pattern) if p.is_file()}
    # copyfile follows links, which hid the symlinked reports of D-025's rehearsal.
    linked = sorted(p.relative_to(source).as_posix() for p in chosen if p.is_symlink())
    if linked:
        raise ValueError(f"published as symlinks into work/, so any other copy dangles: {linked}")
    if bams := sorted(p.name for p in chosen if p.suffix == ".bam"):
        raise ValueError(f"a pattern copies an alignment, which never leaves the cluster: {bams}")
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
    # The run's own date, so a fixture copied again from the same run keeps it.
    run_date = re.search(r"params_(\d{4}-\d{2}-\d{2})", params_file.name).group(1)
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
        "| Publish override | `publish.config`: TRIMGALORE also publishes "
        "`*_trimming_report.json` (D-025) |",
        f"| nf-core/rnaseq | {found.get('nf-core/rnaseq', 'not found')} |",
        f"| Nextflow | {found.get('Nextflow', 'not found')} |",
        f"| Salmon | {found.get('salmon', 'not found')} |",
        f"| Profile | {profile} |",
        f"| Host | {platform.system()} {platform.release()} {platform.machine()} |",
        f"| Run date (from the params file) | {run_date} |",
        f"| Copied (UTC) | {datetime.now(UTC).date().isoformat()} |",
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


def launches(out: Path) -> list[str]:
    """One line per launch of the run: its params file and the keys that shape the fixture."""
    keys = ("aligner", "pseudo_aligner", "skip_alignment", "skip_preseq", "publish_dir_mode")
    return [f"| `{p.name}` | " + ", ".join(f"`{k}` = {json.dumps(json.loads(p.read_text())
                                                              .get(k))}" for k in keys) + " |"
            for p in sorted((out / "pipeline_info").glob("params_*.json"))]  # fmt: skip


def write_aligned_readme(out: Path, command: str, profile: str, omitted: list) -> None:
    """Provenance and scope of the aligned fixture, then the manifest over every file.

    What was left out is summed per top-level folder: the run publishes over a thousand files.
    """
    text = (out / "pipeline_info" / "nf_core_rnaseq_software_mqc_versions.yml").read_text()
    # A tool can run at two versions in one release (samtools 1.21 inside STAR_ALIGN, 1.24 in
    # SAMTOOLS_*), so every version is kept.
    seen: dict[str, set[str]] = {}
    for tool, version in re.findall(r"^\s+(\S+): (\S+)$", text, re.MULTILINE):
        seen.setdefault(tool, set()).add(version)
    found = {tool: ", ".join(sorted(versions)) for tool, versions in seen.items()}
    log = next(out.glob("multiqc/*/multiqc_report_data/multiqc.log")).read_text()
    multiqc = re.search(r"This is MultiQC v(\S+)", log)
    folders: dict[str, list[int]] = {}
    for name, size in omitted:
        folders.setdefault(name.split("/")[0], []).append(size)
    lines = [
        "# nf-core/rnaseq fixture: STAR and Salmon",
        "",
        "A real nf-core/rnaseq run on its test profile with `aligner: star_salmon`, copied as "
        "written, never edited. Built by `scripts/python/make_rnaseq_fixture.py --aligned`.",
        "",
        "| Field | Value |",
        "| --- | --- |",
        f"| Command | `{command}` |",
        "| Publish override | `publish.config`: TRIMGALORE also publishes "
        "`*_trimming_report.json` (D-025) |",
        *(f"| {tool} | {found.get(tool, 'not found')} |" for tool in ALIGNED_TOOLS),
        f"| MultiQC (multiqc.log) | {multiqc.group(1) if multiqc else 'not found'} |",
        f"| Profile | {profile} |",
        f"| Host | {platform.system()} {platform.release()} {platform.machine()} |",
        f"| Copied (UTC) | {datetime.now(UTC).date().isoformat()} |",
        "",
        "## Launches, as their params files record them",
        "",
        "| Params file | Keys |",
        "| --- | --- |",
        *launches(out),
        "",
        "## Left out, per folder",
        "",
        "| Folder | Files | Bytes |",
        "| --- | --- | --- |",
        *(f"| `{name}/` | {len(sizes)} | {sum(sizes)} |" for name, sizes in folders.items()),
        "",
        "## Scope",
        "",
        "Can test: the files M2's read stage and M3 read from an aligned outdir (nf-core's own "
        "MultiQC tables, both Salmon matrices, params and traces), two launches of one outdir, "
        "a sample with two runs merged, single- and paired-end samples, auto and reverse "
        "strandedness.",
        "",
        "Cannot test: GENCODE and `--gencode`; human reads, genes and biotypes; an unstranded "
        "or a failed strand check; the cluster's Nextflow with Apptainer.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines))
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    manifest = "".join(f"{sha256(p)}  {p.relative_to(out).as_posix()}\n" for p in files)
    (out / "MANIFEST.sha256").write_text(manifest)


def aligned_main(args: argparse.Namespace) -> int:
    """Build tests/fixtures/rnaseq_star, from a run here or from ``--results``."""
    version, patterns = definitions(aligned=True)
    if args.results is not None:
        source, command = args.results, args.command or "not recorded"
    else:
        work = args.scratch / "run"
        command = " ".join(run(args.scratch, version, args.profile, aligned=True))
        source = work / "out"
    omitted = copy(source, args.out, patterns)
    write_aligned_readme(args.out, command, args.profile, omitted)
    print(f"{sum(1 for p in args.out.rglob('*') if p.is_file())} files in {args.out}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Run the pipeline once and build the fixture from it."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True, help="outside the repository")
    parser.add_argument("--profile", default=PROFILE, help="emulate_amd64 if an arm64 image fails")
    parser.add_argument(
        "--from-run", action="store_true",
        help="copy from the run already in --scratch/run, as when COPIED_FILES grows",
    )  # fmt: skip
    parser.add_argument("--aligned", action="store_true", help="tests/fixtures/rnaseq_star")
    parser.add_argument("--results", type=Path, help="with --aligned: an outdir already run")
    parser.add_argument("--command", help="with --results: the command that run was launched by")
    args = parser.parse_args(argv)
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError(f"{args.out} is not empty; a fixture is rebuilt into a fresh dir")
    if args.scratch.resolve().is_relative_to(ROOT):
        raise ValueError(f"{args.scratch} is inside the repository")
    if args.aligned:
        return aligned_main(args)
    version, patterns = definitions()
    if args.from_run:
        if not (args.scratch / "run" / "out").is_dir():
            raise FileNotFoundError(f"no run in {args.scratch / 'run' / 'out'}")
        command = command_for(args.scratch / "run", version, args.profile)
    else:
        command = run(args.scratch, version, args.profile)
    omitted = copy(args.scratch / "run" / "out", args.out, patterns)
    write_readme(args.out, command, args.profile, omitted)
    print(f"{sum(1 for p in args.out.rglob('*') if p.is_file())} files in {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
