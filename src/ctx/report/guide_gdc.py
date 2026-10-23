"""The guide's TCGA part: any GDC project from fetch to `ctx run`, on a laptop or Kebnekaise.

Every launcher line is a line of ``scripts/bash/fetch_gdc.sbatch``'s Usage block (``USAGE``),
filled with the dataset's project and slug where the dataset file has ``[gdc]``, else left as
placeholders for a student to fill; the ``ctx`` lines are built from option names the CLI
parses; the dataset file starts as ``conf/dataset_template_gdc.toml``. tests/test_guide.py holds
each line to its source. The page text is in ``content/guide.toml`` under ``[gdc]``.
"""

import functools
import shlex
from pathlib import Path

from ctx.dataset import Dataset, load
from ctx.report import text
from ctx.report.builder import Code, Group, ResultTable

PAGE = "guide"
REPO = Path(__file__).resolve().parents[3]
TEMPLATE = REPO / "conf" / "dataset_template_gdc.toml"
SCRIPT = "fetch_gdc.sbatch"
# (script, its Usage lines) per block, each line as the block holds it.
USAGE = {
    "laptop": (SCRIPT, (
        'export CTX_REPO=<your clone of the course repository> RAW="$HOME/ctx_<slug>/data/raw"',
        'git clone <course repository> "$CTX_REPO"',
        'pixi install --manifest-path "$CTX_REPO/pixi.toml" -e fetch --frozen',
        'CTX_SITE=laptop DATA=0 bash "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" <project> "$RAW/<slug>"',  # noqa: E501
    )),
    "laptop_download": (SCRIPT, (
        "# its last line is the size at GDC: read it, then download",
        'CTX_SITE=laptop bash "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" <project> "$RAW/<slug>"',
    )),
    "cluster": (SCRIPT, (
        'RAW="$CTX_STORE/ctx_<slug>/data/raw"',
        'DATA=0 bash "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" <project> "$RAW/<slug>"',
    )),
    "cluster_download": (SCRIPT, (
        "# its last line is the size at GDC: read it, then download",
        'sbatch -A <compute project> "$CTX_REPO/scripts/bash/fetch_gdc.sbatch" <project> "$RAW/<slug>"',  # noqa: E501
    )),
    "cluster_done": (SCRIPT, (
        "# once that job has ended: every check passed, or the problems it names",
        "tail -n 3 slurm-RNAseq-Tutorial_NA_gdc-download_NA-<job id>.out",
    )),
}  # fmt: skip
CLUSTER = "<user>@kebnekaise.hpc2n.umu.se"
STORE = "<your folder on the storage project>"
KINDS = ("mirna", "isoform")


def names(spec: Dataset) -> tuple[str, str]:
    """The project and slug a block names: the dataset's own, else the Usage's placeholders."""
    if spec.gdc is None:
        return "<project>", "<slug>"
    return spec.gdc.project, spec.slug


def filled(line: str, spec: Dataset) -> str:
    """A Usage line with the dataset's project and slug in place of their placeholders."""
    project, slug = names(spec)
    return line.replace("<project>", project).replace("<slug>", slug)


def dataset_file(spec: Dataset) -> str:
    """The dataset file as the blocks name it: in the clone's conf/datasets, by its slug."""
    if spec.gdc is not None and spec.path.resolve().is_relative_to(REPO):
        return shlex.quote(spec.path.resolve().relative_to(REPO).as_posix())
    return f"conf/datasets/{names(spec)[1]}.toml"


def course_lines(spec: Dataset) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    """The dataset file from the template; setup and the check; then the run, from the clone."""
    slug = names(spec)[1]
    root, given = f"~/ctx_{slug}", dataset_file(spec)
    raw = f"{root}/data/raw/{slug}"
    folders = ("--gdc-dir", raw, *(a for k in KINDS for a in (f"--{k}-dir", f"{raw}_{k}")))
    run = ("pixi run -e laptop ctx --root", root, "run", "--dataset", slug, "--dataset-file",
           given, *folders)  # fmt: skip
    start = ('cd "$CTX_REPO"', "pixi install -e laptop --frozen",
             f"cp -n {TEMPLATE.relative_to(REPO).as_posix()} {given}")  # fmt: skip
    check = (f'cd "$CTX_REPO" && pixi run -e laptop ctx --root {root} setup --dataset-file {given}',
             " ".join((*run, "--check")))  # fmt: skip
    return start, check, ('cd "$CTX_REPO"', " ".join(run))


def copy_line(spec: Dataset) -> str:
    """The copy of every folder the job wrote, from the store to the laptop's ctx root."""
    slug = names(spec)[1]
    return f"/usr/bin/rsync -av {CLUSTER}:{STORE}/ctx_{slug}/data/raw/ ~/ctx_{slug}/data/raw/"


def downloads(spec: Dataset) -> ResultTable:
    """Every file ``ctx setup`` downloads for a TCGA project, with its size and host."""
    word = functools.partial(text.word, PAGE)
    resources = spec.resources if spec.gdc is not None else load(TEMPLATE).resources
    rows = tuple((r.name, f"{r.bytes:,}", r.urls[0].split("/")[2]) for r in resources)
    return ResultTable(word("gdc.downloads"), text.words(PAGE, "gdc.columns"), rows,
                       numeric=(1,))  # fmt: skip


def group(spec: Dataset, n: int) -> Group:
    """The TCGA part, numbered ``n``: find, fetch on either computer, then into the course."""
    word = functools.partial(text.word, PAGE)

    def code(key: str, lines: tuple[str, ...]) -> Code:
        """One block, titled from the page text, its lines ending as a file's do."""
        return Code(word(f"gdc.{key}"), "\n".join(lines) + "\n")

    def usage(name: str) -> tuple[str, ...]:
        """One block's Usage lines, filled for this dataset."""
        return tuple(filled(line, spec) for line in USAGE[name][1])

    start, check, run = course_lines(spec)
    items = (
        word("gdc.find"), word("gdc.pairs"),
        code("laptop", usage("laptop")), code("laptop_download", usage("laptop_download")),
        word("gdc.cluster_login"), code("cluster", usage("cluster")),
        code("cluster_download", usage("cluster_download")),
        code("cluster_done", usage("cluster_done")), code("copy", (copy_line(spec),)),
        word("gdc.template"), code("start", start), word("gdc.edit"), downloads(spec),
        code("check", check), word("gdc.check_reads"), code("run", run), word("gdc.genecard"),
    )  # fmt: skip
    return Group(text.fill(PAGE, "gdc.group", n=n), items, lead=word("gdc.lead"))
