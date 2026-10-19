"""`ctx run --check` on a GDC dataset: what the folders say a TCGA project's dataset file must hold.

The folders are read as gather reads them, every file but no MD5 (gather checks those when it
runs), and the sheet is built in memory by ``gdc_run.build_sheet``, so each contrast's n is the
one ``ctx de`` will count. Each problem carries the lines that fix it, from
``conf/dataset_template_gdc.toml`` or measured from the files; nothing is written. The Table 3
counts are the CDR's own, printed so the student can hold them against Liu et al. 2018's paper.
"""

import difflib
import json
import re
import tomllib
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import polars as pl

from ctx.clinical import harmonise
from ctx.dataset import Contrast, Dataset, DatasetError, GdcSpec
from ctx.gather import gdc, gdc_mirna, gdc_run
from ctx.ncrna import targets
from ctx.pipeline import TEMPLATE, Options, Problem, realised
from ctx.qc import registry

GDC_TEMPLATE = TEMPLATE.with_name("dataset_template_gdc.toml")
# The registration time a [sanity] table records, as conf/datasets/tcga_kirc.toml writes it.
REGISTERED = re.compile(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}")
ZEROS = (0, 0)
# Table 3's column order, so the counts read as the paper prints them.
TABLE3 = ("OS", "PFI", "DFI", "DSS")
# The table M10's endpoint sentence reads Table 3's Use marks from.
TABLE3_USE = "table3_use"
# The gene or miRNA a [sanity] metric names, as M5 and M8 record them.
GENE_METRIC = re.compile(r"^gene\.[^.]+\.(.+)\.(?:lfc|padj)$")
MIRNA_METRIC = re.compile(r"^(?:mirna\.(.+)\.(?:lfc|padj)|pair\.([^.]+)\..+)$")
FAMILIES = "targetscan_families"


def table(heading: str, **values: int) -> str:
    """A table of the GDC template, from its heading to the next blank line, n set to ``values``.

    Raises:
        DatasetError: If the template or the table is not there.
    """
    if not GDC_TEMPLATE.is_file():
        raise DatasetError(f"no {GDC_TEMPLATE}: the lines to add come from the course repository")
    found = re.search(rf"^{re.escape(heading)}\n.*?(?=\n\n|\Z)", GDC_TEMPLATE.read_text(),
                      re.M | re.S)  # fmt: skip
    if found is None:
        raise DatasetError(f"{GDC_TEMPLATE} has no table {heading}")
    lines = found.group(0)
    if values:
        lines = lines.replace("n = [0, 0]", f"n = [{values['reference']}, {values['tested']}]")
    return lines + "\n"


def template_contrast() -> Contrast:
    """The GDC template's tumour-against-normal contrast, its n unset."""
    (entry,) = tomllib.loads(table("[[contrasts]]"))["contrasts"]
    return Contrast(entry["name"], entry["variable"], entry["reference"], entry["tested"],
                    expected_n=ZEROS, design=entry["design"], paired=entry["paired"])  # fmt: skip


def _sheet(path: Path, spec: GdcSpec) -> tuple[pl.DataFrame, set[str]]:
    """The STAR sheet gather would write for one folder, and its genes' symbols; no MD5s."""
    metric = registry.load(registry.resolve(None)).metrics[gdc.STRAND_METRIC]
    folder = gdc.read_folder(path, spec)
    counts = gdc.read_counts(folder, spec)
    notes = gdc.annotations(path / gdc.ANNOTATIONS, folder.files)
    sheet = gdc_run.build_sheet(folder, counts, notes, spec, metric,
                                gdc.library_selection(path, spec))  # fmt: skip
    return sheet, set(counts.genes["gene_name"].drop_nulls())


def _metrics(rules: Iterable[Mapping[str, Any]]) -> list[str]:
    """Every metric the conditions of a [sanity]-shaped table's rules name."""
    return [str(c.get("metric", "")) for r in rules for c in r.get("conditions", [])]


def _unknown(where: str, names: Iterable[str], known: set[str], what: str) -> list[Problem]:
    """A problem per name that is not among ``known``, with the closest names that are."""
    problems = []
    for name in dict.fromkeys(names):
        if name not in known:
            close = difflib.get_close_matches(name, sorted(known), n=3)
            hint = f"; the closest are {', '.join(close)}" if close else ""
            problems.append(Problem(f"{where}: {name!r} is not {what}{hint}"))
    return problems


def _names(root: Path, spec: Dataset, genes: set[str]) -> list[Problem]:
    """Gene symbols and miRNA names the dataset file names, against the ones the files spell.

    A misspelt name is never an error later: the page leaves the gene out, and a rule on it is
    NOT MEASURED. Genes are the counts files' own symbols; miRNAs are TargetScan's human
    miRBase names, the ones M8 reads, once ``ctx setup`` has fetched the family file.
    """
    symbols = [m.group(1) for m in map(GENE_METRIC.match, _metrics(spec.sanity.get("rules", [])))
               if m]  # fmt: skip
    problems = _unknown("[clinical] genes", spec.clinical.genes if spec.clinical else (), genes,
                        "a gene symbol of the counts files")  # fmt: skip
    problems += _unknown("[[sanity.rules]]", symbols, genes, "a gene symbol of the counts files")
    families = [r for r in spec.resources if r.name == FAMILIES]
    if spec.ncrna is None or not families or not (root / families[0].dest).is_file():
        return problems
    known = set(targets.read_names(root / families[0].dest).values())
    metrics = _metrics(spec.ncrna.sanity.get("rules", []))
    found = [m.group(1) or m.group(2) for m in map(MIRNA_METRIC.match, metrics) if m]
    what = "a human miRBase name in TargetScan's family file"
    return problems + _unknown("[ncrna] mirnas and rules", [*spec.ncrna.mirnas, *found], known,
                               what)  # fmt: skip


def _mirna_sheet(options: Options, spec: GdcSpec, star: Path, rna: pl.DataFrame) -> pl.DataFrame:
    """The miRNA sheet gather would write, beside the STAR sheet of the same project."""
    assert options.mirna_dir is not None and options.isoform_dir is not None  # noqa: S101 - checked by the caller
    given = ((gdc_mirna.MIRNA, options.mirna_dir), (gdc_mirna.ISOFORM, options.isoform_dir))
    folders = [gdc_mirna.read_kind(p.resolve(), spec, kind) for kind, p in given]
    found = gdc_mirna.read_matrices(*folders)
    notes = gdc.annotations(star / gdc.ANNOTATIONS, found.files)
    return gdc_mirna.sheet(found, notes, rna, spec)


def _release(path: Path, spec: GdcSpec) -> list[Problem]:
    """The release the folder's status names against the dataset file's."""
    found = json.loads((path / gdc.STATUS).read_text()).get("data_release")
    if found == spec.release:
        return []
    return [Problem(f"[gdc] release is {spec.release!r}; {path.name}/{gdc.STATUS} names "
                    f"{found!r}", f'release = "{found}"')]  # fmt: skip


def _cdr(path: Path, spec: Dataset, gdc_spec: GdcSpec) -> list[Problem]:
    """The CDR's type, and its Table 3 counts against the ones the dataset file copied."""
    if not gdc_spec.cdr_file or not gdc_spec.cdr_type:
        return []
    try:
        rows = gdc.cdr(path / gdc_spec.cdr_file, gdc_spec.cdr_type)
    except gdc.GdcError as error:
        guess = gdc_spec.project.split("-", 1)[-1]
        try:
            gdc.cdr(path / gdc_spec.cdr_file, guess)
        except gdc.GdcError:
            return [Problem(f"[gdc] cdr_type: {error}")]
        return [Problem(f"[gdc] cdr_type: {error}; the CDR has rows of type {guess!r}",
                        f'cdr_type = "{guess}"')]  # fmt: skip
    flags = harmonise.flags(rows)
    counted = {e: flags[e] for e in TABLE3}
    shown = "; ".join(f"{e} {n[0]} events, {n[1]} censored" for e, n in counted.items())
    copied = spec.clinical.published.get("table3", {}) if spec.clinical else {}
    given = {e: tuple(copied.get(e, ZEROS)) for e in counted}
    if all(n == ZEROS for n in given.values()):
        return [Problem(f"[clinical.published.table3]: copy Liu et al. 2018 Table 3's row for "
                        f"{gdc_spec.cdr_type} from the paper; the CDR's own counts, which that "
                        f"row must equal, are {shown}")]  # fmt: skip
    differ = [e for e in counted if given[e] != counted[e]]
    if not differ:
        return []
    return [Problem(f"[clinical.published.table3]: {', '.join(differ)} differ from the CDR's own "
                    f"counts, {shown}; read Table 3's {gdc_spec.cdr_type} row again",
                    blocks=False)]  # fmt: skip


def _uses(spec: Dataset) -> list[Problem]:
    """Table 3's Use marks as M10 reads them, with the template's lines where any is blank.

    Blank, M10 still runs and says the dataset file holds no copy; a mark M10 cannot read
    stops its page. Either stops the run here, as Table 3's zeros do.
    """
    assert spec.gdc is not None  # noqa: S101 - check() calls this for [gdc] alone
    if spec.clinical is None or not spec.gdc.cdr_file:
        return []
    kind = spec.gdc.cdr_type or spec.gdc.project
    blank, fault = harmonise.unread_uses(spec.clinical.published.get(TABLE3_USE, {}))
    if blank:
        lines = table(f"[clinical.published.{TABLE3_USE}]").replace("XXXX", kind)
        return [Problem(f"[clinical.published.{TABLE3_USE}]: {', '.join(blank)} blank; copy Liu "
                        f"et al. 2018 Table 3's Use marks for {kind} from the paper, or M10 says "
                        "the dataset file holds no copy", lines)]  # fmt: skip
    return [Problem(f"[clinical.published.{TABLE3_USE}]: {fault}")] if fault else []


def _n(sheet: pl.DataFrame, spec: Dataset, contrasts: list[Contrast], where: str) -> list[Problem]:
    """Each contrast's declared n against the pairs the chosen files give."""
    problems = []
    for c in contrasts:
        n = realised(gdc.chosen(sheet), c, spec)
        if n != c.expected_n:
            declared = f"{where} {c.name}: n = {list(c.expected_n)} declared"
            fix = f"set n = [{n[0]}, {n[1]}]"
            problems.append(Problem(f"{declared}, {list(n)} in the files; {fix}"))
    return problems


def _missing(spec: Dataset) -> list[Problem]:
    """Tables the GDC template holds that the dataset file lacks, and rules not registered."""
    problems = []
    absent = (("[clinical]", spec.clinical, "clinical has no survival to compute"),
              ("[de]", spec.de, "de names no featured contrast"),
              ("[function]", spec.function, "function tests no gene sets"),
              ("[ncrna]", spec.ncrna, "ncrna reads nothing"))  # fmt: skip
    for heading, parsed, what in absent:
        if parsed is None:
            problems.append(Problem(f"no {heading}: {what}", table(heading), blocks=False))
    if spec.ncrna is not None and not spec.ncrna.mirnas:
        problems.append(Problem("[ncrna] mirnas is empty: M8 tests no miRNA's targets",
                                blocks=False))  # fmt: skip
    rules = spec.sanity.get("rules", [])
    registered = str(spec.sanity.get("registered_on", ""))
    if not rules:
        text = "no [[sanity.rules]]: the hub decides no rule; register each, with its source"
        problems.append(Problem(f"{text}, before the first ctx run", blocks=False))
    elif not REGISTERED.fullmatch(registered):
        problems.append(Problem(f"[sanity] registered_on is {registered!r}: the time the rules "
                                "were added, as YYYY-MM-DD HH:MM", blocks=False))  # fmt: skip
    return problems


def check(spec: Dataset, options: Options, root: Path) -> list[Problem]:
    """Every problem a GDC dataset's folders show, each with its lines where lines fix it.

    Folders gather would refuse for their project are ``pipeline._projects``'s to name; this
    reads only folders matched to a declared project.
    """
    assert spec.gdc is not None  # noqa: S101 - pipeline.check calls this for [gdc] alone
    by, unmatched = gdc.by_project([f for f in options.gdc_dirs if f.is_dir()], spec.gdc.projects)
    if unmatched or not by:
        return []
    problems = _missing(spec) + _uses(spec)
    sheets, genes = [], set()
    for project in spec.gdc.projects:
        one, path = spec.gdc.one(project), by[project]
        if release := _release(path, one):
            problems += release
            continue
        try:
            sheet, symbols = _sheet(path, one)
            sheets.append(sheet)
            genes |= symbols
            problems += _cdr(path, spec, one)
        except (gdc.GdcError, DatasetError) as error:
            problems.append(Problem(str(error)))
    if len(sheets) != len(spec.gdc.projects):
        return problems
    problems += _names(root, spec, genes)
    combined = pl.concat(sheets, how="diagonal_relaxed")
    sheet = combined.with_columns(run_accession=gdc.run_accession())
    if not spec.contrasts:
        n = realised(gdc.chosen(sheet), template_contrast(), spec)
        problems.append(Problem("no [[contrasts]]: de compares nothing",
                                table("[[contrasts]]", reference=n[0], tested=n[1])))  # fmt: skip
    problems += _n(sheet, spec, list(spec.contrasts), "[[contrasts]]")
    mirna = spec.ncrna.mirna if spec.ncrna else None
    if mirna is not None and options.mirna_dir and options.isoform_dir and len(sheets) == 1:
        try:
            small = _mirna_sheet(options, spec.gdc, by[spec.gdc.project], sheets[0])
            problems += _n(small, spec, [mirna], "[ncrna.mirna]")
        except gdc.GdcError as error:
            problems.append(Problem(str(error)))
    return problems
