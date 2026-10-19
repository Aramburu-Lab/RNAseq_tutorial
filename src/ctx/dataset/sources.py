"""``[[resources]]``, ``[columns.<name>]`` and ``[sdrf]``: what curation reads, and how."""

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ctx.dataset.checks import TABLE_SUFFIXES, _need

# Curation columns every dataset with an SDRF must fill (harmonise.CURATION).
REQUIRED_COLUMNS = ("donor", "cell_type", "condition")


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
        if values is not None and not isinstance(values, Mapping):
            problems.append(f"{where}.values: expected a table, found {values!r}")
            values = None
        elif values is not None and not all(isinstance(v, str) for v in values.values()):
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
    paired = isinstance(pairs, list) and all(isinstance(p, list) and len(p) == 2 for p in pairs)  # noqa: PLR2004 - pairs
    if key == "resource" and named is not None and named not in {r.name for r in resources}:
        problems.append(f"sdrf.resource: no [[resources]] entry named {named!r}")
    if key == "file" and named is not None and not named.endswith(TABLE_SUFFIXES):
        problems.append(f"sdrf.file: {named!r} is read by its ending, one of {TABLE_SUFFIXES}")
    if len(join) != 2 or not paired:  # noqa: PLR2004 - pairs
        problems.append("sdrf.join and each sdrf.cross_check entry must be [sdrf header, column]")
    resource, file = (named or "", "") if key == "resource" else ("", named or "")
    pairs = pairs if paired else []
    return Sdrf(resource, (join[0], join[-1]), tuple((a, b) for a, b in pairs), file)
