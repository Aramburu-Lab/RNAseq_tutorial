"""The one reader of a dataset file (``conf/datasets/<slug>.toml``).

A dataset file says where a cohort's metadata comes from and how to read it: the SDRF join,
which columns become curation columns and how their values are spelled, the batch rule, and
the counts the curated sheet must reproduce. Every schema error is reported at once.
"""

import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

# Curation columns every dataset with an SDRF must fill (harmonise.CURATION).
REQUIRED_COLUMNS = ("donor", "cell_type", "condition")
TABLES = {
    "dataset", "resources", "guard", "sdrf", "columns", "derived", "order", "expect",
    "contrasts", "de", "function", "sanity", "gdc", "clinical", "proteo", "ncrna", "counts",
}  # fmt: skip
# The count column a GDC STAR - Counts file holds per library protocol (ctx.gather.gdc).
PROTOCOLS = ("unstranded", "stranded_first", "stranded_second")
# An own table's endings: comma-separated, or tab-separated as an SDRF is (ctx.gather.curation).
TABLE_SUFFIXES = (".csv", ".tsv", ".txt")
# TCGA-CDR's four endpoints (Liu et al. 2018), each an event flag and a time in days.
ENDPOINTS = ("OS", "DSS", "DFI", "PFI")
# M7's designs: one cohort measured at both layers, RNA and protein from two studies, or one
# cohort's tumours against their paired normals at both layers.
PROTEO_DESIGNS = ("across_samples", "cross_dataset", "tumour_normal")
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


class DatasetError(ValueError):
    """A dataset file that cannot be used; the message lists every problem."""


@dataclass(frozen=True)
class Resource:
    """A file ``ctx setup`` fetches: its mirrors, size, checksum and place under the root.

    ``tls_ciphers`` is an OpenSSL cipher string added to Python's default TLS 1.2 list for
    this resource alone, for a server that offers no cipher the default list holds;
    certificate and host checks are unchanged.
    """

    name: str
    urls: tuple[str, ...]
    bytes: int
    sha256: str
    dest: str
    tls_ciphers: str = ""


@dataclass(frozen=True)
class Column:
    """A curated column: the SDRF header it comes from, its unit, meaning and value map.

    ``values`` maps each SDRF value, whitespace-collapsed, to its curated spelling; None
    passes values through unchanged.
    """

    name: str
    source: str
    unit: str
    meaning: str
    values: Mapping[str, str] | None


@dataclass(frozen=True)
class Sdrf:
    """How SDRF rows meet the fetchngs runs: the join pair and the cross-checked pairs.

    The table is a fetched ``[[resources]]`` entry, or ``file``: the student's own table, CSV
    or TSV, named relative to the dataset file, for a study with no SDRF.
    """

    resource: str
    join: tuple[str, str]
    cross_check: tuple[tuple[str, str], ...]
    file: str = ""


@dataclass(frozen=True)
class Split:
    """A seeded split of individuals into two arbitrary groups, balanced within strata."""

    seed: int
    strata: str


@dataclass(frozen=True)
class Contrast:
    """One comparison: which samples, which two levels, which design, and the n it must have.

    Attributes:
        name: Identifier, used in output names.
        variable: Column holding the compared levels; ignored for a split.
        reference: Baseline level; the LFC is tested over reference.
        tested: The other level.
        expected_n: (reference, tested) samples after every step; a mismatch stops the fit.
        design: pydeseq2 formula over ``group`` and sheet columns, e.g. ``~donor + group``.
        filters: Column to allowed values; a sample must match every one.
        paired: Keep only donors with a sample at both levels.
        split: Replace the compared levels by a seeded split of individuals.
        label: The contrast in plain words, as titles, leads, captions and axes name it; empty
            reads as its two levels (Rytis 2026-10-01: never the identifier on a page).
    """

    name: str
    variable: str
    reference: str
    tested: str
    expected_n: tuple[int, int]
    design: str = "~group"
    filters: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    paired: bool = False
    split: Split | None = None
    label: str = ""


@dataclass(frozen=True)
class PcaSpec:
    """A PCA figure and its PERMANOVA metrics.

    Attributes:
        name: Identifier, used in output and metric names.
        filters: Column to allowed values.
        factors: Sheet columns to test; the first colours the points.
        symbol: Column setting the point shape, or None.
        strata: Permute labels only within this column (the individual), or None.
        one_per: Keep one sample per value of this column (the first run), or None.
    """

    name: str
    filters: Mapping[str, tuple[str, ...]]
    factors: tuple[str, ...]
    symbol: str | None = None
    strata: str | None = None
    one_per: str | None = None


@dataclass(frozen=True)
class DeSpec:
    """What the M5 page draws beyond each contrast's tables."""

    featured: str
    pca: tuple[PcaSpec, ...] = ()
    lessons: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    scatters: tuple[tuple[str, str], ...] = ()
    box_contrast: str | None = None
    box_genes: tuple[str, ...] = ()


@dataclass(frozen=True)
class FunctionSpec:
    """Which gene sets M6 tests, their licence line, and what the page features."""

    gene_sets: str
    licence: str
    featured_set: str
    running: tuple[str, str]


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
class CountsSpec:
    """Where a count matrix from a repository comes from, for ``ctx gather --counts``.

    Attributes:
        source: The repository, accession and file, as the outputs' READMEs name them.
        url: Where the file was downloaded from; empty when it is the student's own.
        annotation: The annotation release the source says its gene IDs come from, quoted
            from the source; empty where it names none.
    """

    source: str
    url: str = ""
    annotation: str = ""


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


@dataclass(frozen=True)
class GeneClass:
    """A gene class M7 compares: the union of named sets from one GMT resource."""

    resource: str
    sets: tuple[str, ...]


@dataclass(frozen=True)
class ProteoSpec:
    """What M7 reads for a dataset: the protein table, how its samples meet RNA, the classes.

    Attributes:
        design: One of ``PROTEO_DESIGNS``.
        protein: The ``[[resources]]`` entry holding the protein table.
        label: How a page names the protein source; the resource name where empty.
        member: The table's file inside that resource, when the resource is an archive.
        sample_prefix: Text before a case ID in the protein table's sample names.
        plex: The annotation column naming each sample's TMT plex.
        tracks: Annotation columns drawn beside the paired heatmap.
        amplified: The annotation column marking the worked example's amplification.
        genes: The worked example's genes.
        contrast: The M5 contrast whose fold changes meet the protein's (cross dataset).
        tested: A regular expression naming the protein table's tested-group columns.
        reference: The same for the reference group.
        classes: Gene classes by name.
        published: Numbers a page prints beside its own, each table with its ``source``.
        sanity: The ``[proteo.sanity]`` table: ``registered_at`` and ``[[rules]]``; not
            ``registered``, the key whose commits mark the Step-1 [sanity] registration.
        descriptors: A TSV protein table's leading row-descriptor columns, the gene symbol
            first; the sample columns follow them (tumour against normal).
        plex_map: The ``[[resources]]`` workbook naming each aliquot's TMT plex.
        plex_rows: A regular expression the plex rows' folder names match, the whole proteome's
            and not another enrichment's.
        unit: The protein table's values as a page names them; empty for a log2 ratio.
        positive_fdr: The BH level a per-gene correlation needs to count as positive.
    """

    design: str
    protein: str
    label: str = ""
    member: str | None = None
    sample_prefix: str = ""
    plex: str | None = None
    tracks: tuple[str, ...] = ()
    amplified: str | None = None
    genes: tuple[str, ...] = ()
    contrast: str | None = None
    tested: str | None = None
    reference: str | None = None
    classes: Mapping[str, GeneClass] = field(default_factory=dict)
    published: Mapping[str, Any] = field(default_factory=dict)
    sanity: Mapping[str, Any] = field(default_factory=dict)
    descriptors: tuple[str, ...] = ()
    plex_map: str | None = None
    plex_rows: str = ""
    unit: str = ""
    positive_fdr: float = 0.05


@dataclass(frozen=True)
class NcrnaSpec:
    """What M8 builds: the contrast it stratifies, the pairs its select offers, its own rules.

    Attributes:
        featured: The contrast whose M5 results and samples M8 reads.
        pairs: How many lncRNA and neighbour pairs the page's select offers.
        sanity: M8's pre-registered rules, a ``[sanity]``-shaped table with its
            ``registered_on`` time; apart from ``[sanity]``, whose record is fixed
            (tests/test_sanity_record.py), which finds its registration by any line starting
            ``registered = ``, so this one is named otherwise.
        mirna: The miRNA contrast of a GDC dataset (``[ncrna.mirna]``), fitted on the mature
            miRNA counts gather wrote; None for a dataset without small-RNA libraries.
        mirnas: The mature miRNAs, by miRBase name, whose targets the page tests and pairs.
    """

    featured: str
    pairs: int
    sanity: Mapping[str, Any] = field(default_factory=dict)
    mirna: Contrast | None = None
    mirnas: tuple[str, ...] = ()


@dataclass(frozen=True)
class Dataset:
    """A parsed dataset file."""

    path: Path
    slug: str
    info: Mapping[str, str]
    resources: tuple[Resource, ...]
    sdrf: Sdrf | None
    columns: tuple[Column, ...]
    batch_bins: Mapping[str, tuple[int, ...]]
    order: Mapping[str, tuple[str, ...]]
    expect_runs: int | None
    expect_counts: Mapping[str, Mapping[str, int]]
    contrasts: tuple[Contrast, ...] = ()
    de: DeSpec | None = None
    function: FunctionSpec | None = None
    sanity: Mapping[str, Any] = field(default_factory=dict)
    gdc: GdcSpec | None = None
    clinical: ClinicalSpec | None = None
    proteo: ProteoSpec | None = None
    ncrna: NcrnaSpec | None = None
    counts: CountsSpec | None = None

    @property
    def name(self) -> str:
        """The name pages show: ``[dataset] label`` where given, else the slug.

        The slug names files, so it is lowercase with underscores; the label need not be.
        """
        return str(self.info.get("label", self.slug))

    @property
    def accessions(self) -> tuple[str, ...]:
        """``[dataset] accessions``: what fetchngs downloads, in the file's order."""
        return tuple(self.info.get("accessions", ()))

    def table(self) -> Path | None:
        """The own table of ``[sdrf] file``, beside the dataset file, never the working folder."""
        return self.path.parent / self.sdrf.file if self.sdrf and self.sdrf.file else None

    def resource(self, name: str) -> Resource:
        """The resource with this name."""
        return next(r for r in self.resources if r.name == name)

    def parent(self) -> "Dataset | None":
        """The dataset this one is a subset of (``[dataset] subset_of``), or None for a whole one.

        Its file sits beside this one, named by its slug, as ``conf/datasets`` holds both; a
        page names it by its own ``name``, its label: M3's subset note named the slug, "tier1",
        where the parent's label is "Step-1" (kit v6, 2026-10-07).

        Raises:
            DatasetError: If that file is missing, unusable or names another slug.
        """
        slug = self.info.get("subset_of")
        if not slug:
            return None
        path = self.path.with_name(f"{slug}.toml")
        if not path.is_file():
            raise DatasetError(f"{self.path.name}: subset_of {slug!r}, but no {path} beside it")
        found = load(path)
        if found.slug != slug:
            raise DatasetError(f"{path}: slug {found.slug!r}, not the subset_of {slug!r}")
        return found


def _need(table: Mapping[str, Any], key: str, kind: type, where: str, problems: list[str]) -> Any:
    """Read one required key of a given type, recording a problem instead of raising."""
    value = table.get(key)
    # bool is an int to isinstance; a count written as true is still a mistake.
    if not isinstance(value, kind) or (isinstance(value, bool) and kind is not bool):
        problems.append(f"{where}.{key}: expected {kind.__name__}, found {value!r}")
        return None
    return value


def _resources(raw: list[Any], problems: list[str]) -> tuple[Resource, ...]:
    """Parse [[resources]]."""
    out = []
    for i, entry in enumerate(raw):
        where = f"resources[{i}]"
        fields = [
            _need(entry, "name", str, where, problems),
            _need(entry, "urls", list, where, problems),
            _need(entry, "bytes", int, where, problems),
            _need(entry, "sha256", str, where, problems),
            _need(entry, "dest", str, where, problems),
        ]
        if None not in fields:
            name, urls, size, digest, dest = fields
            if Path(dest).is_absolute() or ".." in Path(dest).parts:
                problems.append(f"{where}.dest: {dest!r} must stay under the root")
            ciphers = entry.get("tls_ciphers", "")
            if not isinstance(ciphers, str):
                problems.append(f"{where}.tls_ciphers: expected str, found {ciphers!r}")
                ciphers = ""
            out.append(Resource(name, tuple(urls), size, digest, dest, ciphers))
    return tuple(out)


def _columns(raw: Mapping[str, Any], problems: list[str]) -> tuple[Column, ...]:
    """Parse [columns.<name>]."""
    out = []
    for name, entry in raw.items():
        where = f"columns.{name}"
        source = _need(entry, "source", str, where, problems)
        meaning = _need(entry, "meaning", str, where, problems)
        values = entry.get("values")
        if values is not None and not all(isinstance(v, str) for v in values.values()):
            problems.append(f"{where}.values: every mapped value must be a string")
        if source is not None and meaning is not None:
            out.append(Column(name, source, entry.get("unit", "none"), meaning, values))
    return tuple(out)


def _sdrf(raw: Mapping[str, Any], resources: tuple[Resource, ...], problems: list[str]) -> Sdrf:
    """Parse [sdrf]: a fetched resource, or the student's own table beside the dataset file."""
    if ("resource" in raw) == ("file" in raw):
        problems.append("sdrf: give resource (a fetched SDRF) or file (your own table), not both")
    key = "file" if "file" in raw else "resource"
    named = _need(raw, key, str, "sdrf", problems)
    join = _need(raw, "join", list, "sdrf", problems) or ["", ""]
    pairs = raw.get("cross_check", [])
    if key == "resource" and named is not None and named not in {r.name for r in resources}:
        problems.append(f"sdrf.resource: no [[resources]] entry named {named!r}")
    if key == "file" and named is not None and not named.endswith(TABLE_SUFFIXES):
        problems.append(f"sdrf.file: {named!r} is read by its ending, one of {TABLE_SUFFIXES}")
    if len(join) != 2 or any(len(p) != 2 for p in pairs):  # noqa: PLR2004 - pairs
        problems.append("sdrf.join and each sdrf.cross_check entry must be [sdrf header, column]")
    resource, file = (named or "", "") if key == "resource" else ("", named or "")
    return Sdrf(resource, (join[0], join[-1]), tuple((a, b) for a, b in pairs), file)


def _values(raw: Mapping[str, Any]) -> dict[str, tuple[str, ...]]:
    """A filter table: column to allowed values."""
    return {k: tuple(str(x) for x in v) for k, v in raw.items()}


def _contrasts(raw: list[Any], problems: list[str]) -> tuple[Contrast, ...]:
    """Parse [[contrasts]]."""
    out = []
    for i, entry in enumerate(raw):
        where = f"contrasts[{i}]"
        fields = [_need(entry, key, str, where, problems)
                  for key in ("name", "variable", "reference", "tested")]  # fmt: skip
        n = _need(entry, "n", list, where, problems)
        split = entry.get("split")
        if n is not None and (len(n) != 2 or not all(isinstance(x, int) for x in n)):  # noqa: PLR2004 - a pair
            problems.append(f"{where}.n: expected [reference n, tested n], found {n!r}")
        if None in fields or n is None:
            continue
        name, variable, reference, tested = (str(f) for f in fields)
        out.append(Contrast(
            name, variable, reference, tested, expected_n=(n[0], n[-1]),
            design=entry.get("design", "~group"),
            filters=_values(entry.get("filters", {})), paired=bool(entry.get("paired", False)),
            split=Split(int(split["seed"]), str(split["strata"])) if split else None,
            label=str(entry.get("label", "")),
        ))  # fmt: skip
    names = [c.name for c in out]
    problems += [
        f"contrasts: {n!r} appears twice" for n in sorted(set(names)) if names.count(n) > 1
    ]
    return tuple(out)


def _de(raw: Mapping[str, Any], names: set[str], problems: list[str]) -> DeSpec:
    """Parse [de] and check it names only declared contrasts."""
    pca = tuple(
        PcaSpec(
            p["name"],
            _values(p.get("filters", {})),
            tuple(p["factors"]),
            p.get("symbol"),
            p.get("strata"),
            p.get("one_per"),
        )  # fmt: skip
        for p in raw.get("pca", [])
    )
    spec = DeSpec(
        featured=raw.get("featured", ""),
        pca=pca,
        lessons={k: tuple(v) for k, v in raw.get("lessons", {}).items()},
        scatters=tuple((a, b) for a, b in raw.get("scatters", [])),
        box_contrast=raw.get("box_contrast"),
        box_genes=tuple(raw.get("box_genes", [])),
    )
    used = [spec.featured, *(c for v in spec.lessons.values() for c in v)]
    used += [c for pair in spec.scatters for c in pair] + (
        [spec.box_contrast] if spec.box_contrast else []
    )
    problems += [
        f"de: {c!r} is not a declared contrast" for c in dict.fromkeys(used) if c not in names
    ]
    return spec


def _function(raw: Mapping[str, Any], names: set[str], resources: tuple[Resource, ...],
              problems: list[str]) -> FunctionSpec | None:  # fmt: skip
    """Parse [function]."""
    fields = [_need(raw, k, str, "function", problems)
              for k in ("gene_sets", "licence", "featured_set")]  # fmt: skip
    running = raw.get("running", [])
    if len(running) != 2:  # noqa: PLR2004 - a pair
        problems.append("function.running: expected [contrast, gene set]")
        return None
    if fields[0] is not None and fields[0] not in {r.name for r in resources}:
        problems.append(f"function.gene_sets: no [[resources]] entry named {fields[0]!r}")
    if running[0] not in names:
        problems.append(f"function.running: {running[0]!r} is not a declared contrast")
    if None in fields:
        return None
    gene_sets, licence, featured = (str(f) for f in fields)
    return FunctionSpec(gene_sets, licence, featured, running=(running[0], running[1]))


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


def _counts(raw: Mapping[str, Any], problems: list[str]) -> CountsSpec | None:
    """Parse [counts]: the matrix's source is required, its URL and annotation are not."""
    source = _need(raw, "source", str, "counts", problems)
    optional = {k: raw.get(k, "") for k in ("url", "annotation")}
    problems += [f"counts.{k}: expected str, found {v!r}" for k, v in optional.items()
                 if not isinstance(v, str)]  # fmt: skip
    problems += [f"counts.{k}: unknown key" for k in sorted(set(raw) - {"source", *optional})]
    if source is None or not all(isinstance(v, str) for v in optional.values()):
        return None
    return CountsSpec(source, str(optional["url"]), str(optional["annotation"]))


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


def _proteo(raw: Mapping[str, Any], names: set[str], resources: tuple[Resource, ...],
            problems: list[str]) -> ProteoSpec | None:  # fmt: skip
    """Parse [proteo]: each design's own keys, its resources and its contrast checked."""
    design = _need(raw, "design", str, "proteo", problems)
    protein = _need(raw, "protein", str, "proteo", problems)
    known = {r.name for r in resources}
    if design is not None and design not in PROTEO_DESIGNS:
        problems.append(f"proteo.design: {design!r}, not one of {PROTEO_DESIGNS}")
    if protein is not None and protein not in known:
        problems.append(f"proteo.protein: no [[resources]] entry named {protein!r}")
    keys = {"across_samples": ("member", "plex"),
            "cross_dataset": ("contrast", "tested", "reference"),
            "tumour_normal": ("member", "plex", "contrast", "descriptors", "plex_map",
                              "plex_rows", "positive_fdr")}  # fmt: skip
    needed = keys.get(design or "", ())
    problems += [f"proteo.{k}: the {design} design needs it" for k in needed if not raw.get(k)]
    if raw.get("contrast") and raw["contrast"] not in names:
        problems.append(f"proteo.contrast: {raw['contrast']!r} is not a declared contrast")
    if raw.get("plex_map") and raw["plex_map"] not in known:
        problems.append(f"proteo.plex_map: no [[resources]] entry named {raw['plex_map']!r}")
    fdr = raw.get("positive_fdr", 0.05)
    if not isinstance(fdr, float) or not 0 < fdr < 1:
        problems.append(f"proteo.positive_fdr: a level in (0, 1), found {fdr!r}")
        fdr = 0.05
    classes = {}
    for name, entry in raw.get("classes", {}).items():
        if entry.get("resource") not in known or not entry.get("sets"):
            problems.append(f"proteo.classes.{name}: needs a known resource and a list of sets")
            continue
        classes[name] = GeneClass(entry["resource"], tuple(entry["sets"]))
    published = raw.get("published", {})
    problems += [f"proteo.published.{k}: no source" for k, v in published.items()
                 if not isinstance(v, dict) or not isinstance(v.get("source"), str)]  # fmt: skip
    if design is None or protein is None:
        return None
    return ProteoSpec(
        design, protein, label=raw.get("label", ""), member=raw.get("member"),
        sample_prefix=raw.get("sample_prefix", ""), plex=raw.get("plex"),
        tracks=tuple(raw.get("tracks", [])), amplified=raw.get("amplified"),
        genes=tuple(raw.get("genes", [])), contrast=raw.get("contrast"),
        tested=raw.get("tested"), reference=raw.get("reference"), classes=classes,
        published=published, sanity=raw.get("sanity", {}),
        descriptors=tuple(raw.get("descriptors", [])), plex_map=raw.get("plex_map"),
        plex_rows=raw.get("plex_rows", ""), unit=raw.get("unit", ""), positive_fdr=fdr,
    )  # fmt: skip


def _ncrna(raw: Mapping[str, Any], names: set[str], problems: list[str]) -> NcrnaSpec | None:
    """Parse [ncrna]."""
    featured = _need(raw, "featured", str, "ncrna", problems)
    pairs = _need(raw, "pairs", int, "ncrna", problems)
    if featured is not None and featured not in names:
        problems.append(f"ncrna.featured: {featured!r} is not a declared contrast")
    if pairs is not None and pairs < 1:
        problems.append(f"ncrna.pairs: {pairs}, at least 1")
    sanity = raw.get("sanity", {})
    if sanity and not isinstance(sanity.get("registered_on"), str):
        problems.append("ncrna.sanity.registered_on: the time the rules were registered")
    mirna = None
    if "mirna" in raw:
        parsed = _contrasts([raw["mirna"]], problems)
        mirna = parsed[0] if parsed else None
    mirnas = raw.get("mirnas", [])
    if not isinstance(mirnas, list) or not all(isinstance(m, str) for m in mirnas):
        problems.append("ncrna.mirnas: a list of mature miRNA names")
        mirnas = []
    if mirnas and mirna is None:
        problems.append("ncrna.mirnas: needs [ncrna.mirna], the contrast they come from")
    if featured is None or pairs is None:
        return None
    return NcrnaSpec(featured, pairs, sanity, mirna, tuple(mirnas))


def load(path: Path) -> Dataset:
    """Read and check a dataset file.

    Args:
        path: The TOML file.

    Returns:
        The parsed dataset.

    Raises:
        DatasetError: Listing every problem found, not only the first.
    """
    raw = tomllib.loads(path.read_text())
    problems = [f"unknown table [{name}]" for name in sorted(set(raw) - TABLES)]
    info = raw.get("dataset", {})
    slug = _need(info, "slug", str, "dataset", problems)
    if "label" in info:
        _need(info, "label", str, "dataset", problems)
    # The runs `ctx guide` writes into the fetchngs list (Part 1), in the file's order.
    accessions = info.get("accessions", [])
    if not isinstance(accessions, list) or not all(isinstance(a, str) for a in accessions):
        problems.append(f"dataset.accessions: expected a list of accessions, found {accessions!r}")
    resources = _resources(raw.get("resources", []), problems)
    columns = _columns(raw.get("columns", {}), problems)
    sdrf = _sdrf(raw["sdrf"], resources, problems) if "sdrf" in raw else None
    named = {c.name for c in columns}
    if sdrf and (missing := [c for c in REQUIRED_COLUMNS if c not in named]):
        problems.append(f"columns: an SDRF dataset must define {missing}")
    bins = raw.get("derived", {}).get("batch", {}).get("bins", {})
    lengths = [n for levels in bins.values() for n in levels]
    if len(lengths) != len(set(lengths)):
        problems.append("derived.batch.bins: a read length sits in two bins")
    order = {k: tuple(v) for k, v in raw.get("order", {}).items()}
    problems += [f"order.{k}: no [columns.{k}]" for k in order if k not in named]
    expect = raw.get("expect", {})
    counts = expect.get("counts", {})
    known = named | {"batch", "read_length"}
    problems += [f"expect.counts.{k}: not a curated column" for k in counts if k not in known]
    contrasts = _contrasts(raw.get("contrasts", []), problems)
    names = {c.name for c in contrasts}
    de = _de(raw["de"], names, problems) if "de" in raw else None
    function = _function(raw["function"], names, resources, problems) if "function" in raw else None
    gdc = _gdc(raw["gdc"], problems) if "gdc" in raw else None
    clinical = _clinical(raw["clinical"], problems) if "clinical" in raw else None
    proteo = _proteo(raw["proteo"], names, resources, problems) if "proteo" in raw else None
    ncrna = _ncrna(raw["ncrna"], names, problems) if "ncrna" in raw else None
    counts_spec = _counts(raw["counts"], problems) if "counts" in raw else None
    if counts_spec and gdc:
        problems.append("[counts] and [gdc] are two routes into gather; keep one")
    if "clinical" in raw and "gdc" not in raw:
        problems.append("[clinical] needs [gdc]: M10 reads a GDC project's gather outputs")
    if proteo and proteo.design == "tumour_normal" and not (gdc and gdc.biospecimens):
        problems.append("proteo: the tumour_normal design needs [gdc] biospecimens, which "
                        "gather reads for each aliquot's case and sample type")  # fmt: skip
    if problems:
        raise DatasetError(f"{path}:\n  " + "\n  ".join(problems))
    return Dataset(
        path=path,
        slug=slug,
        info=info,
        resources=resources,
        sdrf=sdrf,
        columns=columns,
        batch_bins={level: tuple(v) for level, v in bins.items()},
        order=order,
        expect_runs=expect.get("runs"),
        expect_counts={k: {str(v): int(n) for v, n in c.items()} for k, c in counts.items()},
        contrasts=contrasts,
        de=de,
        function=function,
        sanity=raw.get("sanity", {}),
        gdc=gdc,
        clinical=clinical,
        proteo=proteo,
        ncrna=ncrna,
        counts=counts_spec,
    )
