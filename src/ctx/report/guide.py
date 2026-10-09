"""`ctx guide`: the custom-data page, from accessions to a study's own course pages.

Every command on the page comes from one place, and a test holds each to it: the launchers'
lines are their Usage blocks' (``USAGE``, filled with the dataset's accessions and slug), the
``ctx`` lines are built here from option names the CLI parses, and the copy's includes are
``rnaseq.COPIED_FILES``. The dataset file fills Parts 1 and 2; the page text is in
``content/guide.toml``.
"""

import functools
import re
import shlex
from datetime import date
from pathlib import Path

from ctx import runtime
from ctx.dataset import Dataset
from ctx.gather.fetchngs import FETCHNGS_VERSION
from ctx.gather.rnaseq import COPIED_FILES, RNASEQ_VERSION
from ctx.io.naming import MODULES, stem
from ctx.report import text
from ctx.report.builder import (
    Code,
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    ResultTable,
    write,
)

PAGE = "guide"
MODULE = "report"
# The repository this page is built from, for the dataset file's path as a reader types it.
REPO = Path(__file__).resolve().parents[3]
DEFAULT = REPO / "conf" / "datasets" / "four_runs.toml"
# (script, its Usage lines) per block, each line as the block holds it; tests/test_guide.py holds
# them to the scripts. A line with <accession> or <run name, no spaces> is filled below.
USAGE = {
    "fetch_laptop": ("fetchngs.sbatch", (
        'export CTX_REPO=<your clone of the course repository> RUN="$HOME/ctx_runs/<run name, no spaces>"',  # noqa: E501
        'pixi install --manifest-path "$CTX_REPO/pixi.toml" -e fetch --frozen',
        'mkdir -p "$RUN" && cd "$RUN"',
        "printf '%s\\n' <accession> <accession> > ids.txt",
        'CTX_SITE=laptop SKIP_DOWNLOAD=1 bash "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt',  # noqa: E501
    )),
    # Apart from the metadata run, so a pasted block never downloads before the size is read.
    "fetch_laptop_download": ("fetchngs.sbatch", (
        "# its last line is the size at ENA: read it, then download",
        'CTX_SITE=laptop bash "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt',
    )),
    "setup": ("kebnekaise_env.sh", (
        "export CTX_STORE=<your folder on the storage project>",
        'export CTX_REPO="$CTX_STORE/RNAseq_tutorial"',
        'mkdir -p "$CTX_STORE"',
        'git clone <course repository> "$CTX_REPO"',
        "curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v0.76.2 PIXI_NO_PATH_UPDATE=1 bash",
        'PIXI_CACHE_DIR="$CTX_STORE/pixi_cache" "${PIXI_HOME:-$HOME/.pixi}/bin/pixi" install --manifest-path "$CTX_REPO/pixi.toml" -e fetch --frozen',  # noqa: E501
        'source "$CTX_REPO/scripts/bash/kebnekaise_env.sh"',
        "nextflow -version",
    )),
    "fetch_cluster": ("fetchngs.sbatch", (
        'RUN="$CTX_STORE/<run name, no spaces>"',
        'mkdir -p "$RUN" && cd "$RUN"',
        "printf '%s\\n' <accession> <accession> > ids.txt",
        'DRY_RUN=1 bash "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt',
        'SKIP_DOWNLOAD=1 sbatch -A <compute project> "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt',  # noqa: E501
    )),
    "fetch_cluster_size": ("fetchngs.sbatch", (
        "# once that job has ended: the size of the download, read before it is submitted",
        "grep 'bytes of FASTQ at ENA' slurm-ctx_fetchngs-<job id>.out",
    )),
    "fetch_cluster_download": ("fetchngs.sbatch", (
        'sbatch -A <compute project> "$CTX_REPO/scripts/bash/fetchngs.sbatch" "$RUN" ids.txt',
    )),
    "count_laptop": ("rnaseq.sbatch", (
        'REFS="$HOME/ctx_<run name, no spaces>/resources/gencode"',
        'mkdir -p "$REFS" && cd "$REFS"',
        "B=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50",
        'curl -fL --remote-name-all "$B/MD5SUMS" "$B/GRCh38.primary_assembly.genome.fa.gz" "$B/gencode.v50.primary_assembly.annotation.gtf.gz"',  # noqa: E501
        'cd "$RUN" && CTX_SITE=laptop bash "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" "$REFS/GRCh38.primary_assembly.genome.fa.gz" "$REFS/gencode.v50.primary_assembly.annotation.gtf.gz" true',  # noqa: E501
    )),
    "count_cluster": ("rnaseq.sbatch", (
        'mkdir -p "$RUN/refs" && cd "$RUN/refs"',
        "B=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50",
        'curl -fL --remote-name-all "$B/MD5SUMS" "$B/GRCh38.primary_assembly.genome.fa.gz" "$B/gencode.v50.primary_assembly.annotation.gtf.gz"',  # noqa: E501
        'cd "$RUN"',
        'DRY_RUN=1 bash "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" refs/GRCh38.primary_assembly.genome.fa.gz refs/gencode.v50.primary_assembly.annotation.gtf.gz true',  # noqa: E501
        'sbatch -A <compute project> "$CTX_REPO/scripts/bash/rnaseq.sbatch" "$RUN" refs/GRCh38.primary_assembly.genome.fa.gz refs/gencode.v50.primary_assembly.annotation.gtf.gz true',  # noqa: E501
        "# keep nfcore_rnaseq/work until ctx gather on the laptop has read the copied results",
    )),
}  # fmt: skip
CLUSTER = "<user>@kebnekaise.hpc2n.umu.se"
CLUSTER_RUN = "<run folder on Kebnekaise>"
# The outdirs each launcher writes inside its run folder (fetchngs.sbatch, rnaseq.sbatch).
META, RNASEQ = "nfcore_fetchngs_meta/results", "nfcore_rnaseq/results"
# A study of one's own, as the template's [dataset] slug names it.
MINE = "mine"
# The resource Part 2's laptop block saves where ctx setup finds it.
GTF_RESOURCE = "annotation_gtf"
# Measured, the one place the page's sizes come from. The Salmon index's peak resident memory,
# the largest of Step-1's and the breast studies' four cluster builds (GENCODE 50, genome as
# decoys; docs/plans/2026-10-07_custom_data_guide_plan_proposal.md). The laptop run folder of
# U.10 (2026-10-08) before the index finished, with the references saved into results/.
NEEDS = {"index_gib": 16.8, "folder_gb": 61, "work_gb": 37, "results_gb": 22}


def filled(line: str, spec: Dataset) -> str:
    """A Usage line with the dataset's accessions and slug in place of their placeholders."""
    if spec.accessions:
        line = line.replace("<accession> <accession>", " ".join(spec.accessions))
    return line.replace("<run name, no spaces>", spec.slug)


def ctx_line(root: str, command: str, *options: str) -> str:
    """A ``ctx`` command run through the repository's laptop environment, from its folder."""
    return " ".join(["pixi run -e laptop ctx --root", root, command, *options])


def copy_lines(slug: str) -> tuple[str, ...]:
    """The two copies from Kebnekaise: the metadata outdir, then ``COPIED_FILES`` alone.

    macOS's own rsync is named by its path: a newer one earlier on PATH reads the same flags.
    """
    includes = " ".join(f"--include='{p}'" for p in COPIED_FILES)
    raw = f"~/ctx_{slug}/data/raw"
    return (
        f"/usr/bin/rsync -av {CLUSTER}:{CLUSTER_RUN}/{META}/ {raw}/{slug}_meta/",
        f"/usr/bin/rsync -avmL --include='*/' {includes} --exclude='*' "
        f"{CLUSTER}:{CLUSTER_RUN}/{RNASEQ}/ {raw}/{slug}_rnaseq/",
    )  # fmt: skip


def shown(path: Path) -> str:
    """The dataset file as a reader types it: relative to the repository, else from home."""
    path = path.resolve()
    if path.is_relative_to(REPO):
        return shlex.quote(path.relative_to(REPO).as_posix())
    return runtime.retyped([str(path)])


def course_lines(spec: Dataset, dataset_file: str, meta: str, counts: str) -> tuple[str, ...]:
    """Setup, the check, the run and nothing else, from the repository's folder."""
    root = f"~/ctx_{spec.slug}"
    run = ("--dataset", spec.slug, "--dataset-file", dataset_file, "--fetchngs-dir", meta,
           "--rnaseq-dir", counts)  # fmt: skip
    return (
        'cd "$CTX_REPO"',
        "pixi install -e laptop --frozen && pixi install -e reads --frozen",
        ctx_line(root, "setup", "--dataset-file", dataset_file),
        ctx_line(root, "run", *run, "--check"),
        ctx_line(root, "run", *run),
    )


def mine_lines() -> tuple[str, ...]:
    """A study of one's own: the uncurated sheet its table starts from, and the template."""
    meta = f'"$RUN/{META}"'
    return (
        'cd "$CTX_REPO"',
        ctx_line(f"~/ctx_{MINE}", "gather", "--dataset", MINE, "--fetchngs-dir", meta),
        f"cp conf/dataset_template.toml ~/ctx_{MINE}/{MINE}.toml",
    )


def kit_lines(slug: str) -> tuple[str, ...]:
    """The kit builder on the study's own root (U.7)."""
    return (
        'cd "$CTX_REPO"',
        f"pixi run -e laptop python scripts/python/build_kit.py --out ~/ctx_{slug}_kit "
        f"--real {slug}=~/ctx_{slug} --zip",
    )


def blocks(spec: Dataset, dataset_file: Path) -> tuple[Group, ...]:
    """The accession, then the three parts, each a group of code blocks under its lead."""
    word = functools.partial(text.word, PAGE)

    def code(key: str, lines: tuple[str, ...]) -> Code:
        """One block, titled from the page text, its lines ending as a file's do."""
        return Code(word(key), "\n".join(lines) + "\n")

    def usage(name: str) -> tuple[str, ...]:
        """One script's Usage lines, filled for this dataset."""
        return tuple(filled(line, spec) for line in USAGE[name][1])

    given = shown(dataset_file)
    raw = f"~/ctx_{spec.slug}/data/raw/{spec.slug}"
    part3: list[Code | ResultTable | str] = [
        code("part3.laptop", course_lines(spec, given, f'"$RUN/{META}"', f'"$RUN/{RNASEQ}"')),
        code("part3.copy", copy_lines(spec.slug)),
        code("part3.cluster", course_lines(spec, given, f"{raw}_meta", f"{raw}_rnaseq")),
        code("part3.mine", mine_lines()),
        code("part3.kit", kit_lines(spec.slug)),
        word("part3.site"),
    ]
    if caution := spec.info.get("caution"):
        part3.insert(0, str(caution))
    part3.insert(0, downloads(spec))
    parts: tuple[tuple[Code | ResultTable | str, ...], ...] = (
        (code("part1.laptop", usage("fetch_laptop")),
         code("part1.laptop_download", usage("fetch_laptop_download")),
         code("part1.setup", usage("setup")), code("part1.cluster", usage("fetch_cluster")),
         code("part1.cluster_size", usage("fetch_cluster_size")),
         code("part1.cluster_download", usage("fetch_cluster_download"))),
        (text.fill(PAGE, "part2.needs", **NEEDS), code("part2.laptop", usage("count_laptop")),
         code("part2.cluster", usage("count_cluster"))),
        tuple(part3),
    )  # fmt: skip
    numbered = tuple(
        Group(text.fill(PAGE, f"part{n}.group", n=n), items, lead=word(f"part{n}.lead"))
        for n, items in enumerate(parts, start=1)
    )
    return (find(), *numbered)


def find() -> Group:
    """Before Part 1: the accession levels, the routes to a study's runs, the class exercise."""
    word = functools.partial(text.word, PAGE)
    levels = ResultTable(word("find.levels"), text.words(PAGE, "find.columns"),
                         tuple(tuple(row) for row in text.raw(PAGE, "find.rows")))  # fmt: skip
    notes = tuple(word(f"find.{key}") for key in ("archives", "routes", "arrayexpress", "budget",
                                                  "exercise"))  # fmt: skip
    return Group(word("find.group"), (levels, *notes), lead=word("find.lead"))


def downloads(spec: Dataset) -> ResultTable:
    """Every file ``ctx setup`` fetches for the dataset, with its size and where it comes from.

    The GTF is found in place, never fetched again, on the laptop route: Part 2 saved it where
    the dataset file's ``annotation_gtf`` names it.
    """
    word = functools.partial(text.word, PAGE)
    rows = tuple(
        (r.name, f"{r.bytes:,}", r.urls[0].split("/")[2],
         word("part3.found") if r.name == GTF_RESOURCE else word("part3.fetched"))
        for r in spec.resources
    )  # fmt: skip
    return ResultTable(word("part3.downloads"), text.words(PAGE, "part3.columns"), rows,
                       numeric=(1,))  # fmt: skip


def _keys(spec: Dataset) -> tuple[KeyNumber, ...]:
    """Accessions, the two pinned releases and the reference release the Usage downloads."""
    word = functools.partial(text.word, PAGE)
    (release,) = set(re.findall(r"release_(\d+)", "\n".join(USAGE["count_laptop"][1])))
    return (
        KeyNumber(word("keys.runs"), f"{len(spec.accessions):,}"),
        KeyNumber(word("keys.fetchngs"), FETCHNGS_VERSION),
        KeyNumber(word("keys.rnaseq"), RNASEQ_VERSION),
        KeyNumber(word("keys.genome"), text.fill(PAGE, "keys.genome_value", release=release)),
    )


def run(
    settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str, spec: Dataset
) -> Path:
    """Write the guide page for one dataset file under ``outputs/09_report/``.

    Returns:
        The page.
    """
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    out.mkdir(parents=True, exist_ok=True)
    page_stem = stem(MODULE, PAGE, 1, day=date.today())
    usage = runtime.measure(settings.started)
    budget = runtime.budget(MODULE)
    page = ModulePage(
        title=text.fill(PAGE, "title", dataset=spec.name),
        **text.skeleton(PAGE),
        parameters=cards,
        results=blocks(spec, spec.path),
        key_numbers=_keys(spec),
        kicker=text.word(PAGE, "kicker"),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children,
            profile=f"{settings.profile}, {settings.threads} threads",
            budget=text.fill(PAGE, "budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(),
        ),
    )  # fmt: skip
    return write(page, out / f"{page_stem}.html")
