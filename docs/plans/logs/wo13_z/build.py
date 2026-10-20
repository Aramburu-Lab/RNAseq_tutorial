"""wo13-z: rebuild every module page of the real roots in scratch roots, for their sizes.

Run: python build.py <side> <group> [<group> ...], side ``main`` (an export of main's src, its
write refusal lifted so a page over the cap is written and measured) or ``z`` (this
package's src, copied outside git, run as it is). Both sides run from a copy of src outside git,
so their panels name no commit (as Y's cap_compare did); side ``checkout`` runs this worktree's
src in git, as a reader's checkout does. Every write lands in
``outputs/scratch_z/<side>/<group>``: the real root's ``data/`` and its resources other than
``cache/`` are linked, its ``cache/``, ``00_setup`` and gather tables copied, never written.
"""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

WT = Path(__file__).resolve().parents[4]
BIN = "/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.pixi/envs/default/bin"
HOME = Path("/Users/rytis")
SCRATCH = WT / "outputs" / "scratch_z"


def conf(slug: str) -> str:
    """The dataset file of a slug: the worktree's, or the real root's own for a GEO root."""
    own = {"s2_gse327166": HOME / "ctx_s2_gse327166" / "s2_gse327166.toml",
           "s2_gse327906": HOME / "ctx_s2_gse327906" / "s2_gse327906.toml",
           "r2_gse130605": HOME / "ctx_u_gse130605" / "r2_gse130605.toml"}  # fmt: skip
    return str(own.get(slug, WT / "conf" / "datasets" / f"{slug}.toml"))


T1 = HOME / "ctx_tier1" / "data" / "raw"
GENES = {"tcga_kirc": "CA9,NDUFA4L2,VEGFA,MKI67", "tcga_luad": "NKX2-1,SFTPC,EGFR,MKI67",
         "tcga_lusc": "SOX2,TP63,KRT5,MKI67"}  # fmt: skip
# Each group: its real root, the slugs whose gather tables are copied, and the steps in the
# order the real root ran them, with the arguments its pages' Reproduce sections recorded.
GROUPS: dict[str, tuple[str, tuple[str, ...], list[tuple[str, str, list[str]]]]] = {
    "tier1": ("ctx_tier1", (), [
        ("gather", "tier1", ["--fetchngs-dir", f"{T1}/tier1_meta", "--rnaseq-dir",
                             f"{T1}/tier1_rnaseq"]),
        ("qc", "tier1", ["--rnaseq-dir", f"{T1}/tier1_rnaseq"]),
        ("quant", "tier1", ["--rnaseq-dir", f"{T1}/tier1_rnaseq", "--salmon-dir",
                            f"{T1}/tier1_salmon"]),
        ("de", "tier1", []), ("function", "tier1", []), ("proteo", "tier1", []),
        ("ncrna", "tier1", []),
        ("genecard", "tier1", ["--online", "--genes", "ISG15,IFI6,STAT1,OAS2,CD274,AIM2,GBP5"]),
    ]),
    "breast": ("ctx_breast", ("tcga_brca", "cptac_brca"), [
        ("de", "tcga_brca", []), ("function", "tcga_brca", []), ("ncrna", "tcga_brca", []),
        ("clinical", "tcga_brca", []), ("genecard", "tcga_brca", ["--online"]),
        ("de", "cptac_brca", []), ("proteo", "cptac_brca", []),
    ]),
    "kidney": ("ctx_kidney", ("tcga_kirc", "cptac_ccrcc"), [
        ("de", "tcga_kirc", []), ("function", "tcga_kirc", []), ("ncrna", "tcga_kirc", []),
        ("clinical", "tcga_kirc", []), ("genecard", "tcga_kirc", ["--genes", GENES["tcga_kirc"]]),
        ("de", "cptac_ccrcc", []), ("proteo", "cptac_ccrcc", []),
    ]),
    "lung": ("ctx_lung", ("tcga_luad", "tcga_lusc", "tcga_lung"), [
        ("qc", "tcga_lung", ["--strata", "project", "--peers", "project"]),
        ("de", "tcga_lung", []),
        *[(m, s, a) for s in ("tcga_luad", "tcga_lusc") for m, a in (
            ("de", []), ("ncrna", []), ("clinical", []), ("genecard", ["--genes", GENES[s]]))],
    ]),
    # The two commands Y.1 named, for the checkout side: this worktree's src, as committed.
    "step1_qc": ("ctx_tier1", (), [(m, "tier1", a) for m, a in (
        ("gather", ["--fetchngs-dir", f"{T1}/tier1_meta", "--rnaseq-dir", f"{T1}/tier1_rnaseq"]),
        ("qc", ["--rnaseq-dir", f"{T1}/tier1_rnaseq"]))]),
    "brca_clinical": ("ctx_breast", ("tcga_brca",), [("clinical", "tcga_brca", [])]),
    "hnsc": ("ctx_tcga_hnsc", ("tcga_hnsc",), [
        (m, "tcga_hnsc", []) for m in ("de", "function", "ncrna", "clinical")]),
    **{slug: (root, (slug,), [(m, slug, []) for m in ("qc", "de", "function")])
       for slug, root in (("s2_gse327166", "ctx_s2_gse327166"),
                          ("s2_gse327906", "ctx_s2_gse327906"),
                          ("r2_gse130605", "ctx_u_gse130605"))},
}  # fmt: skip

# main's write refusal lifted: write() calls _sized through the module, so this alone changes.
LIFT = (
    "import sys\nfrom ctx.report import builder\n"
    "builder._sized = lambda html, path: html.encode()\n"
    "from ctx.cli import app\napp(sys.argv[1:], standalone_mode=False)\n"
)


def prepare(root: Path, real: Path, slugs: tuple[str, ...]) -> None:
    """The scratch root: data linked, resources linked but cache copied, setup and gather copied."""
    if root.exists():
        sys.exit(f"{root} exists: move it aside first")
    (root / "outputs" / "01_gather").mkdir(parents=True)
    (root / "data").symlink_to(real / "data")
    (root / "resources").mkdir()
    for entry in (real / "resources").iterdir():
        if entry.name == "cache":
            shutil.copytree(entry, root / "resources" / "cache", symlinks=True)
        else:
            (root / "resources" / entry.name).symlink_to(entry)
    shutil.copytree(real / "outputs" / "00_setup", root / "outputs" / "00_setup")
    gather = real / "outputs" / "01_gather"
    for slug in slugs:
        days = sorted(p.name[:10] for p in gather.glob(f"*_gather_{slug}_overview_v1.html"))
        for path in gather.glob(f"{days[-1]}_gather_{slug}_*"):
            shutil.copy2(path, root / "outputs" / "01_gather" / path.name)


def run(side: str, group: str) -> None:
    """Every step of one group on one side, each step's exit, wall time and pages printed."""
    real_name, slugs, steps = GROUPS[group]
    root = SCRATCH / side / group
    prepare(root, HOME / real_name, slugs)
    env = {**os.environ, "PATH": f"{BIN}:{os.environ['PATH']}", "PYTHONDONTWRITEBYTECODE": "1",
           "PYTHONPATH": str(WT / "src" if side == "checkout" else SCRATCH / f"{side}_src" / "src"),
           "CTX_QC_REGISTRY": str(WT / "conf" / "qc_metrics.toml")}  # fmt: skip
    head = ["python", "-c", LIFT] if side == "main" else ["python", "-m", "ctx"]
    for module, slug, args in steps:
        log = root / f"{module}_{slug}.log"
        argv = [*head, "--root", str(root), module, "--dataset", slug,
                "--dataset-file", conf(slug), *args]  # fmt: skip
        started = time.monotonic()
        with log.open("w") as out:
            code = subprocess.call(argv, stdout=out, stderr=subprocess.STDOUT, env=env, cwd=WT)
        print(f"{side} {group} {module} {slug}: exit {code}, {time.monotonic() - started:.0f} s",
              flush=True)  # fmt: skip
    for page in sorted(root.glob("outputs/*/*_overview_v1.html")):
        print(f"  {page.stat().st_size:>11,} B {page.relative_to(root)}", flush=True)


if __name__ == "__main__":
    for name in sys.argv[2:]:
        run(sys.argv[1], name)
