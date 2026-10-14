"""Human Protein Atlas per-gene responses: the XML for numbers, the JSON for categories.

Read from ``https://v<NN>.proteinatlas.org/<ENSG>.{xml,json}``, the versioned host HPA's
licence page asks a citation to use. The XML's ``entry version`` must match the release the
cache files it under. Standard library only, so notebook 09 can run it in Pyodide.
"""

import json
import re
import urllib.parse
import xml.etree.ElementTree as ET  # noqa: S405 - HPA's own XML, fetched over https
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from ctx.genecard import GenecardError

RELEASE = "v25"
ENDPOINTS = ("xml", "json")
# HPA's immune-cell section adds a pooled total beside its 18 cell types: not a cell type.
POOLED = "Total PBMC"
# HPA's JSON keys, as of v25 (read 2026-10-08).
BLOOD_IM = "Blood concentration - Conc. blood IM [pg/L]"
BLOOD_MS = "Blood concentration - Conc. blood MS [pg/L]"
CATEGORIES = ("RNA tissue specificity", "RNA tissue distribution",
              "RNA single cell type specificity", "RNA blood cell specificity",
              "Secretome location", "Reliability (IH)", "Reliability (IF)")  # fmt: skip
LOCATIONS = ("Subcellular main location", "Subcellular additional location")
RELEASE_PATTERN = re.compile(r"v([0-9]+)")
# The search API's cell-line groups, as HPA's search columns name them: cell_RNA_leukemia.
GROUP_PATTERN = re.compile(r"cell_RNA_[a-z_]+")
# One search row's key per cell line of the group asked for (read 2026-10-09).
LINE_KEY = re.compile(r"Cell line RNA - (.+) \[nTPM\]")
# The dataset file's [dataset] key naming an HPA cell-line group and the lines a card shows.
CELL_LINES = "hpa_cell_lines"


def url(release: str, gene: str, endpoint: str) -> str:
    """The versioned URL of one gene's XML or JSON."""
    return f"https://{release}.proteinatlas.org/{gene}.{endpoint}"


def page(release: str, gene: str) -> str:
    """The gene's versioned HPA page, as the licence page asks a citation to give it."""
    return f"https://{release}.proteinatlas.org/{gene}"


def search_url(release: str, query: str, columns: tuple[str, ...]) -> str:
    """HPA's search API: the rows ``query`` finds, with ``columns``, as uncompressed JSON."""
    terms = f"search={urllib.parse.quote(query)}&format=json&columns={','.join(columns)}"
    return f"https://{release}.proteinatlas.org/api/search_download.php?{terms}&compress=no"


@dataclass(frozen=True)
class Xml:
    """The XML's numbers: consensus tissue nTPM, immune cells, single cell types, IHC.

    Attributes:
        version: The XML's ``entry version``.
        tissues: (tissue, nTPM) of the consensus tissue RNA, in the file's order.
        immune: (cell type, lineage, nTPM) of the immune cell RNA, the pooled total included.
        single_cell: (cell type, nCPM) of the single cell type RNA.
        ihc: (tissue, level) of the tissue IHC summary; empty without one.
        reliability: The IHC reliability, or "".
    """

    version: str
    tissues: tuple[tuple[str, float], ...]
    immune: tuple[tuple[str, str, float], ...]
    single_cell: tuple[tuple[str, float], ...]
    ihc: tuple[tuple[str, str], ...]
    reliability: str


def _ntpm(data: ET.Element) -> float:
    """The nTPM level of one ``data`` element."""
    for level in data.iter("level"):
        if level.get("type") == "normalizedRNAExpression":
            return float(level.get("expRNA", "nan"))
    raise GenecardError("an HPA data element without its nTPM")


def parse_xml(body: bytes, release: str) -> Xml:
    """The numbers of one gene's XML.

    Raises:
        GenecardError: If the XML holds no entry, or its version is not ``release``'s.
    """
    entry = ET.fromstring(body).find("entry")  # noqa: S314 - HPA's own XML over https
    if entry is None:
        raise GenecardError("HPA XML without an entry")
    version, wanted = entry.get("version", ""), RELEASE_PATTERN.fullmatch(release)
    if wanted is None or version != wanted[1]:
        raise GenecardError(f"HPA XML version {version!r} filed under release {release!r}")
    sections = {e.get("assayType"): e for e in entry.findall("rnaExpression")}
    tissues = tuple((d.findtext("tissue", ""), _ntpm(d))
                    for d in _data(sections.get("consensusTissue")))  # fmt: skip
    immune = tuple(_cell(d) for d in _data(sections.get("immuneCell")))
    cells = entry.find("cellTypeExpression")
    single = () if cells is None else tuple(
        (e.get("name", ""), float(e.get("expRNA", "nan")))
        for e in cells.findall("singleCellTypeExpression")
    )  # fmt: skip
    ihc: tuple[tuple[str, str], ...] = ()
    reliability = ""
    for block in entry.findall("tissueExpression"):
        if block.get("technology") == "IHC" and block.get("assayType") == "tissue":
            reliability = block.findtext("verification", "").strip()
            ihc = tuple((d.findtext("tissue", ""), d.findtext("level", "").strip())
                        for d in block.findall("data"))  # fmt: skip
    return Xml(version, tissues, immune, single, ihc, reliability)


def _cell(data: ET.Element) -> tuple[str, str, float]:
    """(cell type, lineage, nTPM) of one immune-cell ``data`` element."""
    cell = data.find("immuneCell")
    if cell is None:
        raise GenecardError("an HPA immune-cell data element without its cell type")
    return (cell.text or "").strip(), cell.get("lineage", ""), _ntpm(data)


def _data(section: ET.Element | None) -> list[ET.Element]:
    """A section's ``data`` elements; none for a section the gene lacks."""
    return [] if section is None else section.findall("data")


@dataclass(frozen=True)
class Facts:
    """The JSON's categories and concentrations, each typed once.

    Attributes:
        categories: Each of ``CATEGORIES`` to HPA's word, or "" where HPA gives none.
        locations: Main and additional subcellular locations.
        blood_im: Immunoassay concentration in blood, pg/L, or None.
        blood_ms: Mass-spectrometry concentration in blood, pg/L, or None.
        interactions: HPA's interaction count, or None.
    """

    categories: dict[str, str]
    locations: tuple[tuple[str, ...], tuple[str, ...]]
    blood_im: float | None
    blood_ms: float | None
    interactions: int | None


def _concentration(value: Any) -> float | None:
    """A concentration as HPA writes it, an integer or a string of digits, in pg/L."""
    if value is None or value == "":
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise GenecardError(f"HPA concentration {value!r}: neither a number nor digits")
    return float(value)


def parse_json(body: bytes) -> Facts:
    """The categories, locations, concentrations and interaction count of one gene's JSON.

    Raises:
        GenecardError: If the JSON is not one gene's object.
    """
    raw = json.loads(body)
    if not isinstance(raw, dict):
        raise GenecardError("HPA JSON is not one gene's object")
    categories = {k: str(raw.get(k) or "") for k in CATEGORIES}
    locations = tuple(tuple(raw.get(k) or ()) for k in LOCATIONS)
    count = raw.get("Interactions")
    return Facts(categories, (locations[0], locations[1]), _concentration(raw.get(BLOOD_IM)),
                 _concentration(raw.get(BLOOD_MS)),
                 None if count is None else int(count))  # fmt: skip


def _rows(body: bytes) -> list[dict[str, Any]]:
    """A search response's rows.

    Raises:
        GenecardError: If the response is not a list of objects.
    """
    rows = json.loads(body)
    if not isinstance(rows, list) or not all(isinstance(r, dict) for r in rows):
        raise GenecardError("an HPA search response that is not a list of rows")
    return rows


@dataclass(frozen=True)
class CellLines:
    """The dataset file's HPA cell-line group and the lines its cards name.

    Attributes:
        group: HPA's search column of the group, as ``cell_RNA_leukemia``.
        lines: The lines a card names, as HPA spells them.
    """

    group: str
    lines: tuple[str, ...]


def cell_lines(info: Mapping[str, Any], where: str) -> CellLines | None:
    """``[dataset] hpa_cell_lines``, checked; None where the dataset file names none.

    Raises:
        GenecardError: If the group is not one of HPA's ``cell_RNA_*`` columns, or the lines
            are not a non-empty list of names.
    """
    raw = info.get(CELL_LINES)
    if raw is None:
        return None
    group = raw.get("group") if isinstance(raw, dict) else None
    lines = raw.get("lines") if isinstance(raw, dict) else None
    if not isinstance(group, str) or not GROUP_PATTERN.fullmatch(group):
        raise GenecardError(f"{where}: {CELL_LINES}.group {group!r} is not a cell_RNA_ column")
    if not (isinstance(lines, list) and lines and all(isinstance(x, str) for x in lines)):
        raise GenecardError(f"{where}: {CELL_LINES}.lines must list the lines by name")
    return CellLines(group, tuple(lines))


@dataclass(frozen=True)
class Lines:
    """One gene's nTPM in every line of one HPA cell-line group, in HPA's order."""

    group: str
    values: tuple[tuple[str, float], ...]


def parse_lines(body: bytes, gene: str, group: str) -> Lines:
    """The group's cell-line nTPM in one gene's search row.

    Raises:
        GenecardError: If the response has no single row of ``gene``, or no cell line:
            HPA leaves out a column it does not know, so a misspelt group gives none.
    """
    rows = [r for r in _rows(body) if r.get("Ensembl") == gene]
    if len(rows) != 1:
        raise GenecardError(f"HPA search for {gene}: {len(rows)} rows of that gene, not one")
    values = []
    for key, value in rows[0].items():
        if match := LINE_KEY.fullmatch(key):
            try:
                values.append((match[1], float(value)))
            except (TypeError, ValueError) as e:
                raise GenecardError(f"HPA search for {gene}: {key} is {value!r}") from e
    if not values:
        raise GenecardError(f"HPA search for {gene}: no cell line in group {group!r}")
    return Lines(group, tuple(values))


def parse_names(body: bytes) -> tuple[tuple[str, str, tuple[str, ...]], ...]:
    """(Ensembl ID, HPA symbol, HPA synonyms) of each row of a search by name."""
    return tuple((str(r.get("Ensembl") or ""), str(r.get("Gene") or ""),
                  tuple(str(s) for s in r.get("Gene synonym") or ()))
                 for r in _rows(body))  # fmt: skip
