"""wo13-run: plant each mutation of S.3 and K.12, run the test that names it, restore the file.

A mutation is caught when its test fails. Run through wo13_run_env.sh, with nothing else
running in this checkout: each mutation edits a tracked file for the length of its test.
"""

import shutil
import subprocess
import sys
from pathlib import Path

WT = Path(__file__).resolve().parents[3]
PIPE, GDC, CLI = "src/ctx/pipeline.py", "src/ctx/gather/gdc.py", "src/ctx/cli.py"
KIT, DS = "scripts/python/build_kit.py", "scripts/python/kit_datasets.py"
TWO = "tests/test_pipeline.py::test_two_projects_reach_gather_each_folder"
TOGETHER = "tests/test_pipeline.py::test_check_names_every_projects_folder_problem_together"
REFUSE = "tests/test_pipeline.py::test_check_refuses_gdc_folders_gather_would_refuse"
NEWEST = "tests/test_kit_datasets.py::test_each_version_names_the_newest_commit_holding_it"
REBUILD = "tests/test_kit_datasets.py::test_a_rebuild_reads_the_bytes_the_run_read"
MUTATIONS = [
    (
        "run plans only the first --gdc-dir",
        PIPE,
        'for folder in options.gdc_dirs for a in ("--gdc-dir"',
        'for folder in options.gdc_dirs[:1] for a in ("--gdc-dir"',
        TWO,
    ),
    (
        "ctx run keeps the last --gdc-dir alone",
        CLI,
        "gdc_dirs = tuple(gdc_dir or ())",
        "gdc_dirs = tuple((gdc_dir or [])[-1:])",
        TWO,
    ),
    ("a doubled project passes", GDC, "elif project in by:", "elif False:", TOGETHER),
    (
        "a project with no folder unnamed",
        GDC,
        "for p in projects if p not in by]",
        "for p in projects if False]",
        TOGETHER,
    ),
    (
        "the check stops at the first folder problem",
        PIPE,
        "return [Problem(text) for text in problems]",
        "return [Problem(text) for text in problems[:1]]",
        TOGETHER,
    ),
    (
        "--gdc-dir on a file with no [gdc] passes",
        PIPE,
        "if options.gdc_dirs else []",
        "if False else []",
        REFUSE,
    ),
    (
        "miRNA folders for two projects pass",
        PIPE,
        "if options.mirna_dir is not None and len(spec.gdc.projects) > 1:",
        "if False:",
        REFUSE,
    ),
    (
        "the oldest commit holding a version kept",
        DS,
        "blobs.setdefault(meta.split()[3], (paths[-1], commit))",
        "blobs[meta.split()[3]] = (paths[-1], commit)",
        NEWEST,
    ),
    ("merges skipped (no -m)", DS, '"--full-history", "-m",', '"--full-history",', NEWEST),
    (
        "the kit commit's copy not compared",
        DS,
        "None if head.get(path) == blob else commit",
        "commit",
        NEWEST,
    ),
    (
        "the command still names conf/datasets/<file>",
        DS,
        "mapped = Path(version.path).name if version.commit else version.path",
        "mapped = version.path",
        REBUILD,
    ),
    (
        "git show after the gather command",
        KIT,
        '[*fetches, *([restore(earlier)] if earlier else []), m["command"]]',
        '[*fetches, m["command"], *([restore(earlier)] if earlier else [])]',
        REBUILD,
    ),
    (
        "the commit left out of LEFT_OUT.json",
        KIT,
        '"dataset_file_commit": earlier.commit if earlier else None',
        '"dataset_file_commit": None',
        REBUILD,
    ),
]


def main() -> int:
    """Each mutation planted, tested and restored; exit 1 if any survived."""
    survived = []
    for name, path, old, new, test in MUTATIONS:
        target = WT / path
        original = target.read_text()
        assert original.count(old) == 1, (name, old)
        target.write_text(original.replace(old, new))
        try:
            done = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p",
                                   "no:cacheprovider", test], cwd=WT, capture_output=True,
                                  text=True, check=False)  # fmt: skip
        finally:
            target.write_text(original)
            for cache in WT.glob("**/__pycache__"):
                if ".pixi" not in cache.parts:
                    shutil.rmtree(cache, ignore_errors=True)
        caught = done.returncode != 0
        print(f"{'caught' if caught else 'SURVIVED'}: {name}", flush=True)
        if not caught:
            survived.append(name)
    print(f"{len(MUTATIONS) - len(survived)} of {len(MUTATIONS)} caught")
    return 1 if survived else 0


if __name__ == "__main__":
    sys.exit(main())
