"""`ctx run`: every module a dataset declares, in course order, each run as its own command.

Each module runs in its own process exactly as its command line would run it, so its page keeps
its own Reproduce block and its own peak memory, and the run stops at the first module that
fails, naming it. ``--check`` runs nothing and writes nothing: it reads the dataset file, the
outdirs and the resources, curates the runs in memory and lists every problem at once, each
with the lines to add where lines fix it. Those lines come from ``conf/dataset_template.toml``,
filled from the curated runs.
"""

import re
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ctx import editions, runtime
from ctx.dataset import Contrast, Dataset, DatasetError, load
from ctx.dataset.loader import NEEDS
from ctx.de import results as m5
from ctx.gather.rnaseq import MATRIX_FILE
from ctx.report.text import fill

# Course order, which a new module joins in its place; report reads every page, so it is last.
# Not naming.MODULES, whose numbers name output folders, not the order a dataset is run in.
# genecard is not here: it calls HPA and Ensembl, and one command runs nothing that fetches.
ORDER = ("gather", "qc", "quant", "de", "function", "proteo", "ncrna", "clinical", "report")
# Every module, those ``ctx run`` leaves to the reader too, and the ones each route never runs:
# GDC and count-matrix sheets hold no reads, so quant and process refuse both
# (``ctx.quant.run``, ``ctx.process.run``) and qc runs on a matrix's counts, not on GDC's.
MODULES = (*ORDER, "process", "genecard")
NOT_ON_ROUTE = {"fetchngs": frozenset(), "counts": frozenset({"quant", "process"}),
                "gdc": frozenset({"qc", "quant", "process"})}  # fmt: skip
TEMPLATE = Path(__file__).resolve().parents[2] / "conf" / "dataset_template.toml"
# The curated columns a suggested contrast compares, in the order tried.
VARIABLES = ("condition", "cell_type")
GENE_SETS = "msigdb_hallmark"


@dataclass(frozen=True)
class Options:
    """What `ctx run` was given beside the root."""

    dataset: str
    dataset_file: Path
    fetchngs_dir: Path | None = None
    rnaseq_dir: Path | None = None
    gdc_dirs: tuple[Path, ...] = ()
    salmon_dir: Path | None = None
    mirna_dir: Path | None = None
    isoform_dir: Path | None = None
    counts: Path | None = None

    def dataset_root(self, root: Path) -> Path:
        """Where gather reads the fetchngs outdir when ``--fetchngs-dir`` is not given."""
        return root / "data" / "raw" / self.dataset


@dataclass(frozen=True)
class Problem:
    """One thing a run would meet, and the lines that fix it, if lines do.

    A table the file lacks only leaves its module out of the run (``blocks`` False); anything
    else stops ``ctx run`` before its first module.
    """

    text: str
    lines: str = ""
    blocks: bool = True


def steps(root: Path, spec: Dataset, options: Options) -> list[tuple[str, tuple[str, ...]]]:
    """Each module the dataset declares, with its own options, in ``ORDER``: ``plan``'s first."""
    return plan(root, spec, options)[0]


def plan(
    root: Path, spec: Dataset, options: Options
) -> tuple[list[tuple[str, tuple[str, ...]]], list[str]]:
    """The modules to run, with their own options, in ``ORDER``, and each left out, with why.

    gather and report always run; qc on a fetchngs route; quant with an rnaseq outdir and an
    ``annotation_gtf``, and the salmon folder when given; gather on the GDC route takes every
    ``--gdc-dir``, one per project, and the miRNA and isoform folders when given; de, function,
    proteo, ncrna and clinical with the tables their own commands need (``loader.NEEDS``), and
    de, ncrna and a proteo that is not cross-study only when gather is given counts; a module
    that reads M5's results (``m5.required``) only when de runs too or they are under the root.
    A count matrix (``--counts``) goes to gather alone; qc then judges its counts, and quant,
    which reads an rnaseq outdir, does not run. A module outside the selected edition
    (``ctx.editions``) never runs, and is named among those left out.
    """
    base = ("--dataset", options.dataset, "--dataset-file", str(options.dataset_file))
    rnaseq = ("--rnaseq-dir", str(options.rnaseq_dir)) if options.rnaseq_dir else ()
    fetchngs = ("--fetchngs-dir", str(options.fetchngs_dir)) if options.fetchngs_dir else ()
    gdc = bool(options.gdc_dirs)
    salmon = ("--salmon-dir", str(options.salmon_dir)) if options.salmon_dir else ()
    mirna = tuple(arg for flag, folder in (("--mirna-dir", options.mirna_dir),
                                           ("--isoform-dir", options.isoform_dir))
                  if folder is not None for arg in (flag, str(folder)))  # fmt: skip
    gdc_args = (*(a for folder in options.gdc_dirs for a in ("--gdc-dir", str(folder))), *mirna)
    matrix = ("--counts", str(options.counts)) if options.counts else ()
    counts = gdc or options.rnaseq_dir is not None or bool(matrix)
    has_gtf = any(r.name == "annotation_gtf" for r in spec.resources)
    # These read gather's counts, which this run's gather writes only from a counts input; M5's
    # results, which function, ncrna and cross-study M7 read, come from de or the reader's own.
    reads_counts = {"de", "ncrna", "proteo"}
    if spec.proteo and spec.proteo.design == "cross_dataset":
        reads_counts.discard("proteo")
    skipped = NOT_ON_ROUTE["gdc" if gdc else "counts" if matrix else "fetchngs"]
    planned: dict[str, tuple[str, ...] | None] = {
        "gather": base + (gdc_args if gdc else matrix or fetchngs + rnaseq),
        "qc": None if "qc" in skipped else base + rnaseq,
        "quant": None if "quant" in skipped or not (rnaseq and has_gtf) else base + rnaseq + salmon,
        **{
            module: None if spec.lacks(module) or (module in reads_counts and not counts) else base
            for module in NEEDS
        },
        "report": base,
    }
    edition = editions.active()
    outside = [m for m, args in planned.items() if args is not None and not edition.has_module(m)]
    planned |= dict.fromkeys(outside)
    left_out = [] if planned["de"] else _without_results(root, spec, planned)
    left_out += [fill("run", "outside_edition", module=m, edition=edition.name) for m in outside]
    return [(module, args) for module in ORDER if (args := planned[module]) is not None], left_out


def _without_results(root: Path, spec: Dataset, planned: dict[str, tuple[str, ...] | None]
                     ) -> list[str]:  # fmt: skip
    """Each planned module whose M5 results are not under the root, taken out, and why."""
    reasons = []
    for module in ("function", "proteo", "ncrna"):
        if planned[module] is not None and (absent := m5.absent(root, spec, module)):
            planned[module] = None
            reasons.append(fill("run", "left_out", module=module, slug=spec.slug,
                                contrasts=", ".join(absent)))  # fmt: skip
    return reasons


def planned(root: Path, options: Options) -> list[tuple[str, tuple[str, ...]]]:
    """``steps`` of the dataset file the options name, read again."""
    return steps(root, load(options.dataset_file), options)


def route(spec: Dataset) -> str:
    """The route a dataset file's sheet comes by: ``gdc``, ``counts`` or ``fetchngs``."""
    return "gdc" if spec.gdc else "counts" if spec.counts else "fetchngs"


def route_modules(spec: Dataset) -> set[str]:
    """Every module the dataset's route can run in the selected edition (``ctx.editions``).

    All but ``NOT_ON_ROUTE``'s for that route. Not ``steps``: it also leaves out a module the
    run's inputs cannot feed, while a figure redraws from a reader's own results.
    """
    edition = editions.active()
    return {m for m in MODULES if edition.has_module(m)} - NOT_ON_ROUTE[route(spec)]


def execute(root: Path, profile: str | None, planned: Sequence[tuple[str, tuple[str, ...]]]) -> int:
    """Run each module as its own ``ctx`` command, stopping at the first that fails.

    Returns:
        0, or the failing module's exit code (1 if it gave none).
    """
    head = ["--root", str(root), *(["--profile", profile] if profile else [])]
    for i, (module, args) in enumerate(planned):
        argv = [*head, module, *args]
        # As the module's own page will print it (runtime.retypable).
        print(f"ctx run: {i + 1} of {len(planned)}: {runtime.retypable(argv)}", flush=True)
        done = subprocess.run([sys.executable, "-m", "ctx", *argv], check=False)
        if done.returncode:
            rest = ", ".join(m for m, _ in planned[i + 1 :]) or "none"
            print(f"ctx run: {module} failed with exit {done.returncode}; not run: {rest}",
                  file=sys.stderr, flush=True)  # fmt: skip
            return done.returncode
    return 0


def block(heading: str, /, **values: object) -> str:
    """A commented block of the template, ``# --- <name>`` to the next bare ``#``, filled.

    Raises:
        DatasetError: If the template or the block is not there.
    """
    if not TEMPLATE.is_file():
        raise DatasetError(f"no {TEMPLATE}: the lines to add come from the course repository")
    lines = TEMPLATE.read_text().splitlines()
    try:
        start = lines.index(f"# --- {heading}") + 1
    except ValueError as error:
        raise DatasetError(f"{TEMPLATE} has no block {heading}") from error
    body = []
    for line in lines[start:]:
        if not line.startswith("# "):
            break
        body.append(line[2:])
    return "\n".join(body).format(**values) + "\n"


def realised(sheet: pl.DataFrame, contrast: Contrast, spec: Dataset) -> tuple[int, int]:
    """(reference, tested) samples a contrast compares once its filters and replicates apply."""
    from ctx.de.deseq import select_samples  # noqa: PLC0415 - heavy; only after the curation

    used, _ = select_samples(sheet, contrast, spec.order)
    tested = int((used["group"] == "tested").sum())
    return used.height - tested, tested


def suggested(sheet: pl.DataFrame, spec: Dataset) -> str:
    """A first contrast from the two most frequent levels of a curated column, or empty."""
    for variable in VARIABLES:
        counts = sheet.filter(pl.col(variable) != "").group_by(variable).len()
        top = counts.sort(["len", variable], descending=[True, False])[variable].head(2)
        if top.len() == 2:  # noqa: PLR2004 - a pair of levels
            reference, tested = sorted(top.to_list())
            name = re.sub(r"[^a-z0-9]+", "_", f"{tested}_{reference}".lower()).strip("_")
            draft = Contrast(name, variable, reference, tested, expected_n=(0, 0))
            n = realised(sheet, draft, spec)
            return block("[[contrasts]]", name=name, variable=variable, reference=reference,
                         tested=tested, reference_n=n[0], tested_n=n[1])  # fmt: skip
    return ""


def _contrasts(sheet: pl.DataFrame, spec: Dataset) -> list[Problem]:
    """Every contrast's realised n against its declared n, and the tables M5 and M6 need."""
    if not spec.contrasts:
        lines = suggested(sheet, spec)
        return [
            Problem(
                "no [[contrasts]]: de and function have nothing to compare", lines, blocks=False
            )
        ]
    problems = []
    for c in spec.contrasts:
        n = realised(sheet, c, spec)
        if n != c.expected_n:
            declared = f"{c.name}: n = {list(c.expected_n)} declared, {list(n)} found"
            problems.append(Problem(f"{declared}; set n = [{n[0]}, {n[1]}] in its [[contrasts]]"))
    first = spec.contrasts[0].name
    if spec.de is None:
        lines = block("[de]", contrast=first)
        problems.append(Problem("no [de]: de names no featured contrast", lines, blocks=False))
    if spec.function is None:
        lines = block("[function]", contrast=first)
        if not any(r.name == GENE_SETS for r in spec.resources):
            lines = _resource(GENE_SETS) + "\n" + lines
        problems.append(Problem("no [function]: function tests no gene sets", lines, blocks=False))
    return problems


def _resource(name: str) -> str:
    """The template's ``[[resources]]`` entry of that name, as its lines."""
    text = TEMPLATE.read_text()
    found = re.search(rf'^\[\[resources\]\]\nname = "{name}"\n.*?(?=\n\n)', text, re.M | re.S)
    if found is None:
        raise DatasetError(f"{TEMPLATE} has no [[resources]] entry {name!r}")
    return found.group(0) + "\n"


def _inputs(root: Path, spec: Dataset, options: Options) -> list[Problem]:
    """Outdirs that are not folders, and resources not fetched under the root."""
    problems = []
    default = None if options.gdc_dirs or spec.gdc or options.counts else options.dataset_root(root)
    given = (("--fetchngs-dir", options.fetchngs_dir or default),
             ("--rnaseq-dir", options.rnaseq_dir),
             *(("--gdc-dir", folder) for folder in options.gdc_dirs),
             ("--salmon-dir", options.salmon_dir), ("--mirna-dir", options.mirna_dir),
             ("--isoform-dir", options.isoform_dir))  # fmt: skip
    for flag, folder in given:
        if folder is not None and not folder.is_dir():
            problems.append(Problem(f"{flag} {folder} is not a folder"))
    if options.counts is not None and not options.counts.exists():
        problems.append(Problem(f"--counts {options.counts} is not a file or a folder"))
    problems += _pairing(options) + _projects(spec, options)
    counts = options.rnaseq_dir
    if counts is not None and counts.is_dir() and not (matrix := counts / MATRIX_FILE).is_file():
        problems.append(Problem(f"--rnaseq-dir holds no count matrix, {matrix}: nf-core/rnaseq "
                                "has not finished there, or the copy left it out"))  # fmt: skip
    fetch = runtime.retypable(["--root", str(root), "setup", "--dataset-file",
                               str(options.dataset_file)])  # fmt: skip
    problems += [Problem(f"resource {r.name} is not fetched to {root / r.dest}; run {fetch}")
                 for r in spec.resources if not (root / r.dest).is_file()]  # fmt: skip
    if spec.sdrf is None and spec.gdc is None:
        route = f"route 1 of {TEMPLATE.name} is your own table"
        problems.append(Problem(f"no [sdrf]: no run has a donor or condition; {route}"))
    if (table := spec.table()) is not None and not table.is_file():
        problems.append(Problem(f"no {table}: start it from the sheet ctx gather writes without "
                                "--dataset-file, one row per run"))  # fmt: skip
    return problems


def _pairing(options: Options) -> list[Problem]:
    """The folders a module would refuse or never read, as gather refuses them, before any run."""
    problems = []
    if (options.mirna_dir is None) != (options.isoform_dir is None):
        problems.append(Problem("--mirna-dir and --isoform-dir go together: gather reads both"))
    if (options.mirna_dir or options.isoform_dir) and not options.gdc_dirs:
        problems.append(Problem("--mirna-dir and --isoform-dir take --gdc-dir"))
    if options.salmon_dir is not None and options.rnaseq_dir is None:
        problems.append(Problem("--salmon-dir takes --rnaseq-dir: quant reads both"))
    others = (options.fetchngs_dir, options.rnaseq_dir, *options.gdc_dirs, options.salmon_dir)
    if options.counts is not None and any(f is not None for f in others):
        problems.append(Problem("--counts is a route of its own: no other input folder"))
    return problems


def _projects(spec: Dataset, options: Options) -> list[Problem]:
    """The GDC folders gather would refuse, every one at once, before any run.

    A folder of an undeclared or doubled project, each project with no folder, and miRNA
    folders for several projects. Folders are matched to projects as gather matches them
    (``gdc.by_project``); a folder that is not there is ``_inputs``'s to name.
    """
    if spec.gdc is None:
        given = "--gdc-dir takes a dataset file with a [gdc] table"
        return [Problem(f"{given}; {options.dataset_file} has none")] if options.gdc_dirs else []
    from ctx.gather import gdc  # noqa: PLC0415 - heavy; after the runtime profile

    folders = [folder for folder in options.gdc_dirs if folder.is_dir()]
    _, problems = gdc.by_project(folders, spec.gdc.projects)
    if options.mirna_dir is not None and len(spec.gdc.projects) > 1:
        problems.append("--mirna-dir reads one project; this dataset file names several")
    return [Problem(text) for text in problems]


def _curatable(root: Path, spec: Dataset, options: Options) -> bool:
    """Whether every input curation reads is there; ``_inputs`` names any that is not."""
    if spec.sdrf is None or spec.gdc is not None:
        return False
    table = spec.table() or root / spec.resource(spec.sdrf.resource).dest
    if options.counts is not None:
        return table.is_file() and options.counts.exists() and spec.counts is not None
    source = options.fetchngs_dir or options.dataset_root(root)
    return table.is_file() and source.is_dir()


def _curated(
    root: Path, spec: Dataset, options: Options
) -> tuple[pl.DataFrame | None, list[Problem]]:
    """The sheet as gather would write it, curated in memory, or every reason it is not."""
    from ctx.gather import counts_run, fetchngs, harmonise  # noqa: PLC0415 - after the profile
    from ctx.gather.counts import CountsError  # noqa: PLC0415 - same reason
    from ctx.gather.curation import CurationError  # noqa: PLC0415 - same reason

    if options.counts is not None:
        try:
            return counts_run.check(root, spec, options.counts.absolute()).sheet, []
        except CountsError as error:
            return None, [Problem(f"counts: {c}: {k}: {d}") for c, k, d in error.problems]
        except DatasetError as error:
            return None, [Problem(str(error))]
    source = (options.fetchngs_dir or options.dataset_root(root)).resolve()
    try:
        checked = fetchngs.read(source)
        return harmonise.sheet(checked, harmonise.curated_runs(root, spec, checked, {})), []
    except CurationError as error:
        return None, [Problem(f"curation: {c}: {k}: {d}") for c, k, d in error.problems]
    except (fetchngs.FetchngsError, DatasetError) as error:
        return None, [Problem(str(error))]


def check(root: Path, options: Options) -> tuple[list[Problem], list[str]]:
    """Every problem a run would meet, and the modules it would run, without running any.

    Returns:
        The problems, each with its lines to add where lines fix it, and the planned modules.
    """
    try:
        spec = load(options.dataset_file)
    except (DatasetError, OSError) as error:
        return [Problem(str(error))], []
    if spec.slug != options.dataset:
        return [Problem(f"{options.dataset_file} is for {spec.slug!r}, not --dataset "
                        f"{options.dataset!r}")], []  # fmt: skip
    problems = _inputs(root, spec, options)
    scheduled, left_out = plan(root, spec, options)
    modules = [m for m, _ in scheduled]
    problems += [Problem(reason, blocks=False) for reason in left_out]
    if options.rnaseq_dir and options.salmon_dir and "quant" not in modules:
        problems.append(Problem("--salmon-dir is read by quant alone, which needs an "
                                "annotation_gtf resource in the dataset file"))  # fmt: skip
    if _curatable(root, spec, options):
        sheet, found = _curated(root, spec, options)
        problems += found
        if sheet is not None:
            problems += _contrasts(sheet, spec)
    # The folders are read only once nothing stops the run: a refused folder is named first.
    if spec.gdc is not None and options.gdc_dirs and not any(p.blocks for p in problems):
        from ctx.gather import gdc_check  # noqa: PLC0415 - heavy; reads the GDC folders

        problems += gdc_check.check(spec, options, root)
    return problems, modules
