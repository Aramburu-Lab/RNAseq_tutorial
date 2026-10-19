"""``[proteo]`` and ``[proteo.table]``: what M7 reads for a dataset."""

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ctx.dataset.checks import TABLE_SUFFIXES, _need, _sequence
from ctx.dataset.sources import Resource

# M7's designs: one cohort measured at both layers, RNA and protein from two studies, or one
# cohort's tumours against their paired normals at both layers.
PROTEO_DESIGNS = ("across_samples", "cross_dataset", "tumour_normal")
# [proteo.table]: a protein (or gene) by sample table of any study, read by
# ctx.proteo.quantities through the adapter its format names (ctx.proteo.adapters).
PROTEO_FORMATS = ("table", "maxquant", "diann")
PROTEO_SCALES = ("linear", "log")
PROTEO_LOG_BASES = (2, 10, "e")
PROTEO_ID_KINDS = ("symbol", "uniprot", "ensembl")
PROTEO_ACQUISITIONS = ("tmt", "label_free")
# The missing-value token each tool writes, unless [proteo.table] missing names others: MaxQuant
# writes 0 for a protein not quantified in a sample, DIA-NN leaves the cell empty.
PROTEO_MISSING = {"maxquant": ("0",), "diann": ("",)}


@dataclass(frozen=True)
class GeneClass:
    """A gene class M7 compares: the union of named sets from one GMT resource."""

    resource: str
    sets: tuple[str, ...]


@dataclass(frozen=True)
class ProteoTable:
    """``[proteo.table]``: what a protein (or gene) by sample table holds, and its sample table.

    Attributes:
        format: One of ``PROTEO_FORMATS``: a plain table, or a tool's own export.
        samples: The sample table, a CSV or TSV beside the dataset file: each protein column,
            its RNA sample, its batch.
        scale: ``linear`` for intensities, ``log`` for logarithms of ``log_base``.
        missing: Every cell text that means a missing value; never read as zero.
        id_kind: One of ``PROTEO_ID_KINDS``.
        acquisition: ``tmt`` (the batch is the plex) or ``label_free`` (the batch is the run).
        batch: The sample table's batch column.
        log_base: 2, 10 or ``e``, for ``log`` values.
        id_column: The identifier column; a tool's own by ``id_kind`` where empty.
        file: The protein table beside the dataset file, in place of ``[proteo] protein``.
        quantity: MaxQuant's quantity, the prefix of its sample columns.
        id_map: The ``[[resources]]`` table from UniProt accession to gene symbol.
        id_map_columns: That table's accession and symbol columns.
        sample_column: The sample table's protein column.
        rna_column: The sample table's RNA column; an empty cell is a sample with no RNA.
        reference: A sample table column whose ``yes`` marks each TMT plex's reference channel;
            every value becomes a log2 ratio to its plex's reference, which then leaves.
    """

    format: str
    samples: str
    scale: str
    missing: tuple[str, ...]
    id_kind: str
    acquisition: str
    batch: str
    log_base: int | str = 2
    id_column: str = ""
    file: str = ""
    quantity: str = "LFQ intensity"
    id_map: str = ""
    id_map_columns: tuple[str, str] = ("", "")
    sample_column: str = "protein_sample"
    rna_column: str = "rna_sample"
    reference: str = ""


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
        table: ``[proteo.table]``, a table of any study; None for the CPTAC and Rieckmann
            readers.
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
    table: ProteoTable | None = None


def _plain_file(name: Any, where: str, problems: list[str]) -> None:
    """A file named beside the dataset file: a CSV, TSV or text name with no folder."""
    if not isinstance(name, str) or Path(name).name != name or not name:
        problems.append(f"{where}: {name!r}, a file name beside the dataset file")


def _choice(raw: Mapping[str, Any], key: str, options: tuple[Any, ...], default: Any,
            problems: list[str]) -> Any:  # fmt: skip
    """One key of [proteo.table] that takes one of ``options``."""
    value = raw.get(key, default)
    if value not in options or isinstance(value, bool):
        problems.append(f"proteo.table.{key}: {value!r}, one of {list(options)}")
    return value


def _proteo_table(raw: Mapping[str, Any], known: set[str], problems: list[str]
                  ) -> ProteoTable | None:  # fmt: skip
    """Parse [proteo.table]: each key checked, every problem listed."""
    fmt = _choice(raw, "format", PROTEO_FORMATS, None, problems)
    samples = _need(raw, "samples", str, "proteo.table", problems)
    if samples is not None and not samples.endswith(TABLE_SUFFIXES):
        problems.append(f"proteo.table.samples: {samples!r} ends in one of {TABLE_SUFFIXES}")
    _plain_file(samples, "proteo.table.samples", problems)
    scale = _choice(raw, "scale", PROTEO_SCALES, None, problems)
    base = _choice(raw, "log_base", PROTEO_LOG_BASES, 2, problems)
    kind = _choice(raw, "id_kind", PROTEO_ID_KINDS, None, problems)
    acquisition = _choice(raw, "acquisition", PROTEO_ACQUISITIONS, None, problems)
    batch = _need(raw, "batch", str, "proteo.table", problems)
    missing = raw.get("missing", PROTEO_MISSING.get(fmt or ""))
    if not isinstance(missing, (list, tuple)) or not all(isinstance(m, str) for m in missing):
        problems.append(f"proteo.table.missing: the {fmt} format needs a list of the texts a "
                        f"missing value is written as, found {missing!r}")  # fmt: skip
        missing = ()
    if fmt == "table" and not raw.get("id_column"):
        problems.append("proteo.table.id_column: a plain table names its identifier column")
    if raw.get("file"):
        _plain_file(raw["file"], "proteo.table.file", problems)
    if raw.get("reference") and acquisition != "tmt":
        problems.append("proteo.table.reference: a reference channel belongs to a TMT plex")
    id_map, columns = raw.get("id_map", ""), tuple(raw.get("id_map_columns", ("", "")))
    if kind == "uniprot" and (id_map not in known or len(columns) != 2 or not all(columns)):  # noqa: PLR2004
        problems.append("proteo.table.id_map: UniProt identifiers need a [[resources]] map and "
                        "id_map_columns = [accession column, gene symbol column]")  # fmt: skip
    if None in (fmt, samples, scale, kind, acquisition, batch):
        return None
    return ProteoTable(fmt, samples, scale, tuple(missing), kind, acquisition, batch, base,
                       id_column=raw.get("id_column", ""), file=raw.get("file", ""),
                       quantity=raw.get("quantity", "LFQ intensity"), id_map=id_map,
                       id_map_columns=(columns + ("", ""))[:2],
                       sample_column=raw.get("sample_column", "protein_sample"),
                       rna_column=raw.get("rna_column", "rna_sample"),
                       reference=raw.get("reference", ""))  # fmt: skip


def _proteo_source(raw: Mapping[str, Any], known: set[str], problems: list[str]
                   ) -> tuple[ProteoTable | None, str | None]:  # fmt: skip
    """[proteo.table] and the protein resource: a fetched resource, or a file of one's own."""
    table = _proteo_table(raw["table"], known, problems) if "table" in raw else None
    own = bool(raw.get("table", {}).get("file"))
    if own and "protein" in raw:
        problems.append("proteo: give protein (a fetched resource) or table.file, not both")
    protein = raw.get("protein", "") if own else _need(raw, "protein", str, "proteo", problems)
    if protein and protein not in known:
        problems.append(f"proteo.protein: no [[resources]] entry named {protein!r}")
    return table, protein


def _proteo(raw: Mapping[str, Any], names: set[str], resources: tuple[Resource, ...],
            problems: list[str]) -> ProteoSpec | None:  # fmt: skip
    """Parse [proteo]: each design's own keys, its resources and its contrast checked."""
    design = _need(raw, "design", str, "proteo", problems)
    known = {r.name for r in resources}
    table, protein = _proteo_source(raw, known, problems)
    if design is not None and design not in PROTEO_DESIGNS:
        problems.append(f"proteo.design: {design!r}, not one of {PROTEO_DESIGNS}")
    if "table" in raw and design != "across_samples":
        problems.append("proteo.table: read by the across_samples design alone")
    keys = {"across_samples": () if "table" in raw else ("member", "plex"),
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
        sets = _sequence(entry.get("sets"), f"proteo.classes.{name}.sets", [])
        if entry.get("resource") not in known or not sets:
            problems.append(f"proteo.classes.{name}: needs a known resource and a list of sets")
            continue
        classes[name] = GeneClass(entry["resource"], tuple(sets))
    lists = {k: _sequence(raw.get(k, []), f"proteo.{k}", problems)
             for k in ("tracks", "genes", "descriptors")}  # fmt: skip
    published = raw.get("published", {})
    problems += [f"proteo.published.{k}: no source" for k, v in published.items()
                 if not isinstance(v, dict) or not isinstance(v.get("source"), str)]  # fmt: skip
    if design is None or protein is None:
        return None
    return ProteoSpec(
        design, protein, label=raw.get("label", ""), member=raw.get("member"),
        sample_prefix=raw.get("sample_prefix", ""),
        plex=raw.get("plex") or (table.batch if table else None),
        tracks=tuple(lists["tracks"]), amplified=raw.get("amplified"),
        genes=tuple(lists["genes"]), contrast=raw.get("contrast"),
        tested=raw.get("tested"), reference=raw.get("reference"), classes=classes,
        published=published, sanity=raw.get("sanity", {}),
        descriptors=tuple(lists["descriptors"]), plex_map=raw.get("plex_map"),
        plex_rows=raw.get("plex_rows", ""), unit=raw.get("unit", ""), positive_fdr=fdr,
        table=table,
    )  # fmt: skip
