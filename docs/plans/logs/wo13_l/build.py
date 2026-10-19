"""Package wo13-l: build every module page of one dataset in a scratch root of this worktree.

Usage: <main env python> docs/plans/logs/wo13_l/build.py <label> <dataset[:step,...]> ...

The root, ``outputs/scratch_roots/<root>_<label>``, links the real root's ``data`` and each
entry of its ``resources`` but ``manifest.json``, which is copied, so nothing of the real root
is written. The code is ``outputs/scratch_src/<label>/src`` where it exists (``git archive`` of
main's src for the before pages), else this worktree's src. Setup runs without a dataset
file, so it fetches nothing; genecard runs
``--offline``, reading the response cache and writing none.
"""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
WT = Path(__file__).absolute().parents[4]
HOME = Path.home()
R = {"tier1": HOME / "ctx_tier1", "breast": HOME / "ctx_breast", "kidney": HOME / "ctx_kidney",
     "lung": HOME / "ctx_lung"}  # fmt: skip
STEPS = {
    "tier1": ("tier1", {
        "setup": [],
        "gather": ["--fetchngs-dir", "{data}/raw/tier1_meta", "--rnaseq-dir",
                   "{data}/raw/tier1_rnaseq"],
        "qc": ["--rnaseq-dir", "{data}/raw/tier1_rnaseq"],
        "quant": ["--rnaseq-dir", "{data}/raw/tier1_rnaseq", "--salmon-dir",
                  "{data}/raw/tier1_salmon"],
        "de": [], "function": [], "proteo": [], "ncrna": [],
        "genecard": ["--offline", "--genes", "ISG15,IFI6,STAT1,OAS2,CD274,AIM2,GBP5"],
        "report": []}),
    "cptac_brca": ("breast", {
        "setup": [], "gather": ["--gdc-dir", "{data}/raw/cptac_brca"], "de": [], "proteo": [],
        "report": []}),
    "tcga_brca": ("breast", {
        "gather": ["--gdc-dir", "{data}/raw/tcga_brca", "--mirna-dir", "{data}/raw/tcga_brca_mirna",
                   "--isoform-dir", "{data}/raw/tcga_brca_isoform"],
        "de": [], "function": [], "ncrna": [], "genecard": ["--offline"], "clinical": [],
        "report": []}),
    "cptac_ccrcc": ("kidney", {
        "setup": [], "gather": ["--gdc-dir", "{data}/raw/cptac_ccrcc_star"], "de": [],
        "proteo": [], "report": []}),
    "tcga_lung": ("lung", {
        "setup": [],
        "gather": ["--gdc-dir", "{data}/raw/tcga_luad", "--gdc-dir", "{data}/raw/tcga_lusc"],
        "qc": ["--strata", "project", "--peers", "project"], "de": [], "report": []}),
}  # fmt: skip


def scratch(which: str, label: str) -> Path:
    """The scratch root, its data and resources linked from the real root."""
    real = R[which]
    root = WT / "outputs" / "scratch_roots" / f"{which}_{label}"
    (root / "outputs").mkdir(parents=True, exist_ok=True)
    (root / "resources").mkdir(exist_ok=True)
    if not (root / "data").exists():
        (root / "data").symlink_to(real / "data")
    for entry in (real / "resources").iterdir():
        mine = root / "resources" / entry.name
        if mine.exists() or mine.is_symlink():
            continue
        if entry.name == "manifest.json":
            shutil.copy2(entry, mine)
        else:
            mine.symlink_to(entry)
    return root


def main() -> None:
    """Run each dataset's steps in turn."""
    label, *datasets = sys.argv[1:]
    for one in datasets:
        dataset, _, only = one.partition(":")
        build(label, dataset, only.split(",") if only else [])


def build(label: str, dataset: str, only: list[str]) -> None:
    """Run the dataset's steps, each timed, into one log per dataset and label."""
    which, steps = STEPS[dataset]
    root = scratch(which, label)
    src = WT / "outputs" / "scratch_src" / label / "src"
    env = dict(os.environ, PATH=BIN + os.pathsep + os.environ.get("PATH", ""),
               PYTHONPATH=str(src if src.is_dir() else WT / "src"),
               PYTHONDONTWRITEBYTECODE="1")  # fmt: skip
    log = WT / "docs" / "plans" / "logs" / "wo13_l" / f"build_{dataset}_{label}.log"
    file = WT / "conf" / "datasets" / f"{dataset}.toml"
    for step, extra in steps.items():
        if only and step not in only:
            continue
        args = [BIN + "/python", "-m", "ctx", "--root", str(root), step]
        if step != "setup":
            args += ["--dataset", dataset, "--dataset-file", str(file)]
        args += [a.format(data=root / "data") for a in extra]
        began = time.time()
        with log.open("a") as out:
            out.write(f"== {step} {time.strftime('%Y-%m-%dT%H:%M:%S')} src={env['PYTHONPATH']}\n")
            out.flush()
            done = subprocess.run(args, cwd=root, env=env, stdout=out, stderr=out, check=False)
            out.write(f"exit {done.returncode} {time.time() - began:.1f} s\n")
        print(dataset, step, done.returncode, f"{time.time() - began:.1f} s", flush=True)


main()
