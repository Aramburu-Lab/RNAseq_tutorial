"""GTEx V10 median junction reads per tissue, through GTEx API v2, for the junction panel.

V11 has no per-gene junction route: the API's ``DatasetId`` lists gtex_v8, gtex_snrnaseq_pilot
and gtex_v10 (openapi.json, read 2026-10-09), and no V11 junction file is on disk. So the panel
is V10's, labelled V10 wherever drawn, on V10's GENCODE v39. ``medianJunctionExpression`` takes
a versioned ``gencodeId`` and no ``gencodeVersion``, so the v39 version is first read from
``reference/gene`` with ``gencodeVersion=v39``; ``datasetId`` is always named, since the
routes' defaults differ (gtex_v10 here, gencode v26 for ``reference/gene``). Standard library
and polars only.
"""

import json
import re
import urllib.parse
from typing import Any

import polars as pl

from ctx.genecard import GenecardError
from ctx.report import genes

API = "https://gtexportal.org/api/v2"
RELEASE = "V10"
DATASET = "gtex_v10"
GENCODE = "v39"
GENOME = "GRCh38/hg38"
# The API's largest page; a gene whose junctions need two pages is refused, never cut.
PAGE_ROWS = 100000
# As V10 writes them: chr7:98617888-98619642:+, the intron's first and last base (NPTX2's
# first MANE intron, read 2026-10-09).
_JUNCTION = re.compile(r"(chr[0-9A-Za-z_]+):([0-9]+)-([0-9]+):[+-]")


def gene_url(gene: str) -> str:
    """``reference/gene``: the gene's GENCODE v39 record, which holds its versioned ID."""
    query = urllib.parse.urlencode({"geneId": gene, "gencodeVersion": GENCODE,
                                    "genomeBuild": GENOME})  # fmt: skip
    return f"{API}/reference/gene?{query}"


def junction_url(gencode_id: str) -> str:
    """``medianJunctionExpression`` of one versioned v39 ID in every tissue of gtex_v10."""
    query = urllib.parse.urlencode({"gencodeId": gencode_id, "datasetId": DATASET,
                                    "itemsPerPage": PAGE_ROWS})  # fmt: skip
    return f"{API}/expression/medianJunctionExpression?{query}"


def _data(body: bytes, what: str) -> list[dict[str, Any]]:
    """A paginated response's rows, refused when the response spans more than one page.

    Raises:
        GenecardError: If the body is not GTEx's paginated shape, or holds more pages.
    """
    raw = json.loads(body)
    if not isinstance(raw, dict) or not isinstance(raw.get("data"), list):
        raise GenecardError(f"GTEx {what}: not a paginated response")
    pages = (raw.get("paging_info") or {}).get("numberOfPages", 1)
    if pages > 1:
        raise GenecardError(f"GTEx {what}: {pages} pages; only the first was read")
    return list(raw["data"])


def versioned(body: bytes, gene: str) -> str | None:
    """The gene's versioned v39 ID in a ``reference/gene`` response; None where v39 lacks it.

    Raises:
        GenecardError: If the response names another GENCODE release, or two v39 IDs.
    """
    wanted = genes.number(gene)
    found = set()
    for row in _data(body, "reference/gene"):
        if (version := str(row.get("gencodeVersion", ""))) != GENCODE:
            raise GenecardError(f"GTEx reference/gene: GENCODE {version!r}, not {GENCODE}")
        gencode_id = str(row.get("gencodeId", ""))
        if "_PAR_Y" not in gencode_id and genes.number(gencode_id) == wanted:
            found.add(gencode_id)
    if len(found) > 1:
        raise GenecardError(f"GTEx reference/gene: {sorted(found)} for {gene}")
    return found.pop() if found else None


def parse(body: bytes, gencode_id: str) -> pl.DataFrame:
    """One gene's median junction reads: junction, chrom, start, end, tissue, median.

    Raises:
        GenecardError: If a row is of another dataset or gene, or a junction ID is not
            ``chrom:start-end:strand``.
    """
    rows = []
    for row in _data(body, "medianJunctionExpression"):
        seen = (row.get("datasetId"), row.get("gencodeId"))
        if seen != (DATASET, gencode_id):
            raise GenecardError(f"GTEx junctions: a row of {seen}, not {DATASET}, {gencode_id}")
        match = _JUNCTION.fullmatch(str(row.get("junctionId", "")))
        if match is None:
            raise GenecardError(f"GTEx junction ID {row.get('junctionId')!r}")
        rows.append((match[0], match[1], int(match[2]), int(match[3]),
                     str(row["tissueSiteDetailId"]), float(row["median"])))  # fmt: skip
    schema = {"junction": pl.String, "chrom": pl.String, "start": pl.Int64, "end": pl.Int64,
              "tissue": pl.String, "median": pl.Float64}  # fmt: skip
    return pl.DataFrame(rows, schema=schema, orient="row")
