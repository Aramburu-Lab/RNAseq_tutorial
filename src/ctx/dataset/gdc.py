"""``[gdc]`` and ``[clinical]``: a GDC project folder and what M10 builds from it."""

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any

from ctx.dataset.checks import _need

# The count column a GDC STAR - Counts file holds per library protocol (ctx.gather.gdc).
PROTOCOLS = ("unstranded", "stranded_first", "stranded_second")

# TCGA-CDR's four endpoints (Liu et al. 2018), each an event flag and a time in days.
ENDPOINTS = ("OS", "DSS", "DFI", "PFI")

# The patient-table variables M10 tests expression against and maps for confounding, with the
# tests each takes (ctx.clinical.analyses); a dataset file lists those its patients carry.
CLINICAL_VARIABLES = ("er", "stage", "grade", "histology", "menopause", "age", "sex")
# The two-level variables M10's box can split patients by (ctx.clinical.figures.box).
CLINICAL_BOXES = ("er", "stage_group")
# [gdc]'s optional text keys.
GDC_OPTIONAL = ("primary_site", "cdr_file", "cdr_type", "mirna_workflow", "attribution",
                "biospecimens")  # fmt: skip
# Projects whose own attribution policy the course quotes: the key under ui.toml [gdc.attribution].
ATTRIBUTIONS = ("tcga",)


@dataclass(frozen=True)
class GdcSpec:
    """A GDC project folder: what ``ctx.gather.gdc`` checks it against and how it chooses.

    Attributes:
        project: The GDC project ID the metadata's files must belong to.
        release: The GDC data release, as ``gdc_status.json`` names it.
        workflow: The expression workflow the files come from, e.g. ``STAR - Counts``.
        gene_model: The annotation each counts file's first line names, e.g. ``GENCODE v36``.
        protocol: The library protocol kept, one of ``PROTOCOLS``; files of another are dropped,
            and counts come from this column.
        tumour: The GDC sample type of the tumours; one is chosen per case.
        normal: The GDC sample type of the normal tissue; one is chosen per case.
        primary_site: The case primary site to keep, where the folder's metadata carries it;
            None keeps every file of the project.
        cdr_file: The TCGA-CDR workbook in the folder, or None for a project without one.
        cdr_type: The CDR's ``type`` of this project's patients.
        mirna_workflow: The workflow of the project's miRNA and isoform files, e.g. ``BCGSC
            miRNA Profiling``, which ``ctx gather --mirna-dir`` checks; None for a project
            gathered without them.
        attribution: The project's own attribution policy, a key of ``ATTRIBUTIONS``, which
            pages quote beside the GDC's citation; None where none was read from its source.
        projects: Every project of the dataset, one folder each, in the file's order;
            ``project`` names them together. One unless ``[gdc] projects`` lists several.
        exclude_cases: Cases whose files no analysis takes, by the dataset file's own rule.
        exclude_reason: Why, with its source; required with ``exclude_cases``.
        biospecimens: A PDC biospecimen record in the folder, ``biospecimenPerStudy`` as the
            fetch saved it, giving each protein aliquot's case and sample type; None for none.
    """

    project: str
    release: str
    workflow: str
    gene_model: str
    protocol: str
    tumour: str
    normal: str
    primary_site: str | None = None
    cdr_file: str | None = None
    cdr_type: str | None = None
    mirna_workflow: str | None = None
    attribution: str | None = None
    projects: tuple[str, ...] = ()
    exclude_cases: tuple[str, ...] = ()
    exclude_reason: str = ""
    biospecimens: str | None = None

    def one(self, project: str) -> "GdcSpec":
        """This spec for one of its projects' folders: the project and that project alone."""
        return replace(self, project=project, projects=(project,))


@dataclass(frozen=True)
class ClinicalSpec:
    """What M10 builds for a GDC dataset: the build list, the endpoint and published numbers.

    Attributes:
        genes: The build list, gene symbols the static figures and the default selection hold.
        endpoint: The default endpoint, one of ``ENDPOINTS``, as TCGA-CDR Table 3 recommends.
        published: Numbers a page prints beside its own, each table with its ``source``.
        variables: Patient-table variables the page tests and maps, from
            ``CLINICAL_VARIABLES``, those the cohort's records carry; "er" adds the ER check.
        box: The two-level variable of the expression box, from ``CLINICAL_BOXES``.
        therapy: The GDC treatment type the treatment table counts, e.g. ``Hormone Therapy``.
    """

    genes: tuple[str, ...]
    endpoint: str
    published: Mapping[str, Any] = field(default_factory=dict)
    variables: tuple[str, ...] = ("stage", "age")
    box: str = "stage_group"
    therapy: str = "Chemotherapy"


def _projects(raw: Mapping[str, Any], problems: list[str]) -> tuple[str, ...]:
    """[gdc] project, one string, or projects, a list of two or more, never both."""
    if ("project" in raw) == ("projects" in raw):
        problems.append("gdc: give project (one) or projects (several), not both")
        return ()
    if "project" in raw:
        one = _need(raw, "project", str, "gdc", problems)
        return (one,) if one is not None else ()
    listed = raw["projects"]
    names = isinstance(listed, list) and all(isinstance(x, str) for x in listed)
    if not names or not len(set(listed)) == len(listed) > 1:
        problems.append(f"gdc.projects: two or more distinct project IDs, found {listed!r}")
        return ()
    return tuple(listed)


def _excluded(raw: Mapping[str, Any], problems: list[str]) -> tuple[tuple[str, ...], str]:
    """[gdc] exclude_cases with the exclude_reason it needs."""
    cases, reason = raw.get("exclude_cases", []), raw.get("exclude_reason", "")
    if not isinstance(cases, list) or not all(isinstance(c, str) for c in cases):
        problems.append(f"gdc.exclude_cases: a list of case IDs, found {cases!r}")
        return (), ""
    if len(set(cases)) < len(cases):
        problems.append("gdc.exclude_cases: a case is listed twice")
    if cases and (not isinstance(reason, str) or not reason):
        problems.append("gdc.exclude_reason: exclude_cases needs the reason, with its source")
    return tuple(cases), str(reason)


def _gdc(raw: Mapping[str, Any], problems: list[str]) -> GdcSpec | None:
    """Parse [gdc]."""
    projects = _projects(raw, problems)
    keys = ("release", "workflow", "gene_model", "protocol", "tumour", "normal")
    fields = [" and ".join(projects) or None, *(_need(raw, k, str, "gdc", problems) for k in keys)]
    optional = {k: raw.get(k) for k in GDC_OPTIONAL}
    excluded, reason = _excluded(raw, problems)
    if len(projects) > 1 and (optional["cdr_file"] or optional["mirna_workflow"]):
        problems.append("gdc.projects: a CDR or miRNA folder is read for one project alone")
    problems += [f"gdc.{k}: expected str, found {v!r}" for k, v in optional.items()
                 if v is not None and not isinstance(v, str)]  # fmt: skip
    if fields[4] is not None and fields[4] not in PROTOCOLS:
        problems.append(f"gdc.protocol: {fields[4]!r}, not one of {PROTOCOLS}")
    if (optional["cdr_file"] is None) != (optional["cdr_type"] is None):
        problems.append("gdc.cdr_file and gdc.cdr_type are given together or not at all")
    if optional["attribution"] not in (None, *ATTRIBUTIONS):
        problems.append(f"gdc.attribution: {optional['attribution']!r}, not one of {ATTRIBUTIONS}")
    if None in fields:
        return None
    project, release, workflow, model, protocol, tumour, normal = (str(f) for f in fields)
    return GdcSpec(project, release, workflow, model, protocol, tumour, normal,
                   primary_site=optional["primary_site"], cdr_file=optional["cdr_file"],
                   cdr_type=optional["cdr_type"],
                   mirna_workflow=optional["mirna_workflow"],
                   attribution=optional["attribution"], projects=projects,
                   exclude_cases=excluded, exclude_reason=reason,
                   biospecimens=optional["biospecimens"])  # fmt: skip


def _clinical(raw: Mapping[str, Any], problems: list[str]) -> ClinicalSpec | None:
    """Parse [clinical]."""
    genes = _need(raw, "genes", list, "clinical", problems)
    endpoint = _need(raw, "endpoint", str, "clinical", problems)
    if endpoint is not None and endpoint not in ENDPOINTS:
        problems.append(f"clinical.endpoint: {endpoint!r}, not one of {ENDPOINTS}")
    published = raw.get("published", {})
    problems += [f"clinical.published.{k}: no source" for k, v in published.items()
                 if not isinstance(v, dict) or not isinstance(v.get("source"), str)]  # fmt: skip
    variables = _need(raw, "variables", list, "clinical", problems)
    box = _need(raw, "box", str, "clinical", problems)
    therapy = _need(raw, "therapy", str, "clinical", problems)
    if variables is not None and (bad := [v for v in variables if v not in CLINICAL_VARIABLES]):
        problems.append(f"clinical.variables: {bad}, not in {CLINICAL_VARIABLES}")
    if box is not None and box not in CLINICAL_BOXES:
        problems.append(f"clinical.box: {box!r}, not one of {CLINICAL_BOXES}")
    if box == "er" and variables is not None and "er" not in variables:
        problems.append("clinical.box: 'er' needs 'er' in clinical.variables")
    if None in (genes, endpoint, variables, box, therapy):
        return None
    if not genes or not all(isinstance(g, str) for g in genes):
        problems.append("clinical.genes: a non-empty list of gene symbols")
    return ClinicalSpec(tuple(genes), endpoint, published, tuple(variables), box, therapy)
