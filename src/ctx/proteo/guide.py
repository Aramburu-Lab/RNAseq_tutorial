"""The guide's M7 part: from a PRIDE project or an own export to ``ctx proteo`` beside the RNA.

Every line comes from one place: the ``prepare`` lines from ``prepare.parser``'s subcommands,
the ``[proteo]`` block from the dataset template (``pipeline.block``), the ``ctx`` lines from
the CLI's option names; tests/test_proteo_study.py holds each to its source. The words are
``content/guide.toml``'s ``[mass_spec]``.
"""

import functools

from ctx import pipeline
from ctx.dataset import Dataset
from ctx.report import text
from ctx.report.builder import Code, Group

PAGE = "guide"
PREPARE = "pixi run -e laptop python -m ctx.proteo.prepare"
ACCESSION = "<PXD_accession>"
FILE = "<file_name>"
EXPORT = "<your_export_or_its_archive>"


def lines(spec: Dataset, given: str) -> dict[str, tuple[str, ...]]:
    """Each block's lines, filled with the dataset's slug and file as the reader types them."""
    root = f"~/ctx_{spec.slug}"
    beside = f"{root}/{spec.slug}_proteo.csv"
    ctx = f"pixi run -e laptop ctx --root {root}"
    return {
        "list": (
            'cd "$CTX_REPO"',
            f"{PREPARE} files {ACCESSION} --match txt",
            f"{PREPARE} files {ACCESSION} --match pg_matrix",
        ),
        "fetch": (
            'cd "$CTX_REPO"',
            f"{PREPARE} fetch {ACCESSION} {FILE} --root {root} --slug {spec.slug}",
        ),
        "samples": (
            'cd "$CTX_REPO"',
            f"{PREPARE} samples --format maxquant --table {EXPORT} "
            f"--member proteinGroups.txt --out {beside}",
            f"{PREPARE} samples --format maxquant --table {EXPORT} "
            f'--quantity "Reporter intensity corrected" --out {beside}',
            f"{PREPARE} samples --format diann --table {EXPORT} --out {beside}",
        ),
        "block": tuple(pipeline.block("[proteo]", slug=spec.slug).splitlines()),
        "run": (
            'cd "$CTX_REPO"',
            f"{ctx} setup --dataset-file {given}",
            f"{ctx} proteo --dataset {spec.slug} --dataset-file {given}",
        ),
    }


def part(spec: Dataset, given: str, n: int) -> Group:
    """The M7 part as a group of notes and code blocks under its lead."""
    word = functools.partial(text.word, PAGE)
    found = lines(spec, given)
    items: list[Code | str] = [word("mass_spec.formats"), word("mass_spec.missing")]
    for key in ("list", "fetch", "samples", "block", "run"):
        items.append(Code(word(f"mass_spec.{key}"), "\n".join(found[key]) + "\n"))
    items += [word("mass_spec.rows"), word("mass_spec.check")]
    return Group(text.fill(PAGE, "mass_spec.group", n=n), tuple(items), lead=word("mass_spec.lead"))
