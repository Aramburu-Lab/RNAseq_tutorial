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

from ctx import runtime
from ctx.dataset import Contrast, Dataset, DatasetError, load
from ctx.gather.rnaseq import MATRIX_FILE

# Course order, which a new module joins in its place; report reads every page, so it is last.
# Not naming.MODULES, whose numbers name output folders, not the order a dataset is run in.
# genecard is not here: it calls HPA and Ensembl, and one command runs nothing that fetches.
ORDER = ("gather", "qc", "quant", "de", "function", "proteo", "ncrna", "clinical", "report")
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
    gdc_dir: Path | None = None

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


def steps(spec: Dataset, options: Options) -> list[tuple[str, tuple[str, ...]]]:
    """Each module the dataset declares, with its own options, in ``ORDER``.

    gather and report always run; qc on a fetchngs route; quant with an rnaseq outdir and an
    ``annotation_gtf``; de with contrasts, ``[de]`` and counts; function with ``[function]``
    too; proteo with ``[proteo]`` and counts, and de's tables when it compares fold changes;
    ncrna with ``[ncrna]`` and de's tables; clinical with ``[clinical]``.
    """
    base = ("--dataset", options.dataset, "--dataset-file", str(options.dataset_file))
    rnaseq = ("--rnaseq-dir", str(options.rnaseq_dir)) if options.rnaseq_dir else ()
    fetchngs = ("--fetchngs-dir", str(options.fetchngs_dir)) if options.fetchngs_dir else ()
    gdc = options.gdc_dir is not None
    counts = gdc or options.rnaseq_dir is not None
    has_gtf = any(r.name == "annotation_gtf" for r in spec.resources)
    # Across studies M7 reads M5's fold changes; across samples, gather's counts alone.
    fold_changes = bool(spec.contrasts and spec.de)
    proteo = spec.proteo and counts and (spec.proteo.design == "across_samples" or fold_changes)
    planned: dict[str, tuple[str, ...] | None] = {
        "gather": base + (("--gdc-dir", str(options.gdc_dir)) if gdc else fetchngs + rnaseq),
        "qc": None if gdc else base + rnaseq,
        "quant": base + rnaseq if rnaseq and has_gtf else None,
        "de": base if spec.contrasts and spec.de and counts else None,
        "function": base if spec.function and spec.contrasts and spec.de and counts else None,
        "proteo": base if proteo else None,
        "ncrna": base if spec.ncrna and spec.contrasts and spec.de and counts else None,
        "clinical": base if spec.clinical else None,
        "report": base,
    }
    return [(module, args) for module in ORDER if (args := planned[module]) is not None]


def planned(options: Options) -> list[tuple[str, tuple[str, ...]]]:
    """``steps`` of the dataset file the options name, read again."""
    return steps(load(options.dataset_file), options)


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
    default = None if options.gdc_dir else options.dataset_root(root)
    given = (("--fetchngs-dir", options.fetchngs_dir or default),
             ("--rnaseq-dir", options.rnaseq_dir), ("--gdc-dir", options.gdc_dir))  # fmt: skip
    for flag, folder in given:
        if folder is not None and not folder.is_dir():
            problems.append(Problem(f"{flag} {folder} is not a folder"))
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


def _curatable(root: Path, spec: Dataset, options: Options) -> bool:
    """Whether every input curation reads is there; ``_inputs`` names any that is not."""
    if spec.sdrf is None or spec.gdc is not None:
        return False
    table = spec.table() or root / spec.resource(spec.sdrf.resource).dest
    source = options.fetchngs_dir or options.dataset_root(root)
    return table.is_file() and source.is_dir()


def _curated(
    root: Path, spec: Dataset, options: Options
) -> tuple[pl.DataFrame | None, list[Problem]]:
    """The sheet as gather would write it, curated in memory, or every reason it is not."""
    from ctx.gather import fetchngs, harmonise  # noqa: PLC0415 - heavy; after the runtime profile
    from ctx.gather.curation import CurationError  # noqa: PLC0415 - same reason

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
    if _curatable(root, spec, options):
        sheet, found = _curated(root, spec, options)
        problems += found
        if sheet is not None:
            problems += _contrasts(sheet, spec)
    return problems, [m for m, _ in steps(spec, options)]
