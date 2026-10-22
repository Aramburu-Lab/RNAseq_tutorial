#!/usr/bin/env python3
# /// codebase-meta
# exempt: course tool, not a pipeline step; tabulates each sample of an ENA study as ENA records it
# ///
r"""Every sample of an ENA study with every attribute its submitter recorded, as two TSV tables.

Two ENA APIs, no account and no personal data in any request:
- the portal API's read_run file report lists the study's runs and the sample of each;
- the browser API's XML gives each sample's record, 50 per request, as the sample page shows it.

Tags are kept exactly as ENA writes them (``tissue``, ``source_name``, ``ENA-FIRST-PUBLIC``);
ENA's browser shows them title-cased (``Tissue``, ``Source Name``), the XML does not.

The wide table holds one row per sample: fixed columns (the sample page's names, its runs and
layout from the file report, and its free-text description), then one column per tag in the
order the study first shows it, empty where a sample lacks the tag, repeated tags joined by
``"; "``. The long table holds one row per sample, tag and value, so nothing is lost.
Tab-separated; the csv module quotes a value only when it holds a tab, a newline or a double
quote, which a csv reader would otherwise misread. Each file is written under a temporary name and
renamed when complete.

Usage:
    python scripts/python/ena_samples.py <study> ... --out-dir <folder> [--version N]
        [--label <study>=<label> ...]

With ``--version N`` the files take the deliverable names
``<date>_ena_samples_<study>_vN.tsv`` and ``<date>_ena_samples_<study>_attributes_long_vN.tsv``,
the accession lowercased, as the deliverable rule allows only lowercase; without it,
``<study>_samples.tsv`` and ``<study>_attributes_long.tsv``. ``--label PRJNA000001=lung_tissue``
puts the label before the accession in both, so a name says what the samples are.
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import re
import sys
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
BROWSER_XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
RUN_FIELDS = ("run_accession,sample_accession,secondary_sample_accession,sample_alias,"
              "sample_title,experiment_accession,library_layout")  # fmt: skip
STUDY = re.compile(r"^(PRJ[EDN][A-Z]\d+|[SED]RP\d+)$")
# A label is a description of the deliverable rule, labdocs lint's DELIVERABLE_STEM; copied,
# not imported, because this script runs on the standard library alone.
LABEL = re.compile(r"^[a-z0-9]+(_[a-z0-9]+)*$")
BATCH = 50
TIMEOUT = 60
RETRIES = 2
SERVER_ERROR = 500
FIXED = ("sample_accession", "secondary_accession", "sample_alias", "sample_title",
         "scientific_name", "center_name", "broker_name", "runs", "layout",
         "description")  # fmt: skip
LONG = ("sample_accession", "tag", "value")
# A product name, nothing about the person running it.
AGENT = "ctx-teaching-suite ena_samples.py (urllib)"

Opener = Callable[[str, float], bytes]


class EnaError(RuntimeError):
    """ENA refused or never answered; the message names the accession."""


@dataclass
class Sample:
    """One sample record from ENA's browser XML, with its runs from the file report."""

    accession: str
    secondary: str = ""
    alias: str = ""
    title: str = ""
    scientific_name: str = ""
    center_name: str = ""
    broker_name: str = ""
    description: str = ""
    runs: list[str] = field(default_factory=list)
    layouts: list[str] = field(default_factory=list)
    attributes: list[tuple[str, str]] = field(default_factory=list)


@dataclass
class Requests:
    """HTTP attempts made, so the summary can say how many and how many were retries."""

    made: int = 0
    retried: int = 0


def urlopen(url: str, timeout: float) -> bytes:
    """The body of one GET, with a user agent naming the tool and nothing else."""
    request = urllib.request.Request(url, headers={"User-Agent": AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return bytes(response.read())


def _transient(error: Exception) -> bool:
    """True for a 5xx or a timeout, the failures worth one more try."""
    if isinstance(error, urllib.error.HTTPError):
        return error.code >= SERVER_ERROR
    if isinstance(error, urllib.error.URLError):
        return isinstance(error.reason, TimeoutError)
    return isinstance(error, TimeoutError)


def fetch(url: str, what: str, opener: Opener, counter: Requests, pause: float = 2.0) -> bytes:
    """GET ``url``, retrying a 5xx or timeout up to RETRIES times.

    Args:
        url: The request.
        what: The accession or accessions it asks for, named in any error.
        opener: Performs one GET; ``urlopen`` outside tests.
        counter: Counts every attempt.
        pause: Seconds before the first retry, doubled for the next.

    Returns:
        The response body.

    Raises:
        EnaError: ENA answered 4xx, or failed RETRIES + 1 times.
    """
    for attempt in range(RETRIES + 1):
        counter.made += 1
        try:
            return opener(url, TIMEOUT)
        except (urllib.error.URLError, TimeoutError) as error:
            if not _transient(error):
                raise EnaError(f"{what}: ENA refused the request ({error})") from error
            if attempt == RETRIES:
                raise EnaError(f"{what}: no answer after {RETRIES + 1} tries ({error})") from error
            counter.retried += 1
            time.sleep(pause * 2**attempt)
    raise AssertionError("unreachable")


def parse_runs(study: str, text: str) -> dict[str, Sample]:
    """Samples of a read_run file report, in first-seen order, each with its runs.

    Raises:
        EnaError: The report holds no run, which is how ENA answers an accession it lacks.
    """
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise EnaError(f"{study}: ENA resolves no runs for this accession")
    samples: dict[str, Sample] = {}
    for row in rows:
        sample = samples.setdefault(row["sample_accession"], Sample(row["sample_accession"]))
        sample.runs.append(row["run_accession"])
        if row["library_layout"] not in sample.layouts:
            sample.layouts.append(row["library_layout"])
    return samples


def _text(node: ET.Element, path: str) -> str:
    """The stripped text at ``path`` under ``node``, empty when absent."""
    found = node.find(path)
    return (found.text or "").strip() if found is not None else ""


def parse_sample_xml(body: bytes) -> list[Sample]:
    """Every SAMPLE element of a browser-API SAMPLE_SET, attributes in document order."""
    parsed = []
    for node in ET.fromstring(body).iter("SAMPLE"):
        sample = Sample(
            accession=_text(node, "IDENTIFIERS/PRIMARY_ID") or node.get("accession", ""),
            secondary=_text(node, "IDENTIFIERS/SECONDARY_ID"),
            alias=_text(node, "IDENTIFIERS/SUBMITTER_ID") or node.get("alias", ""),
            title=_text(node, "TITLE"),
            scientific_name=_text(node, "SAMPLE_NAME/SCIENTIFIC_NAME"),
            center_name=node.get("center_name", ""),
            broker_name=node.get("broker_name", ""),
            description=_text(node, "DESCRIPTION"),
        )
        for attribute in node.iter("SAMPLE_ATTRIBUTE"):
            # Verbatim, unstripped: a value is the submitter's, whitespace included.
            tag, value = attribute.findtext("TAG") or "", attribute.findtext("VALUE") or ""
            sample.attributes.append((tag, value))
        parsed.append(sample)
    return parsed


def fetch_study(study: str, opener: Opener, counter: Requests) -> list[Sample]:
    """Every sample of ``study``, each with its runs and its full record.

    Raises:
        EnaError: The study does not resolve, or a sample's record did not come back.
    """
    if not STUDY.match(study):
        raise EnaError(f"{study}: not an ENA study or project accession (PRJ*, SRP, ERP, DRP)")
    url = f"{PORTAL}?accession={study}&result=read_run&fields={RUN_FIELDS}&format=tsv"
    samples = parse_runs(study, fetch(url, study, opener, counter).decode())
    accessions = list(samples)
    returned: set[str] = set()
    for start in range(0, len(accessions), BATCH):
        batch = accessions[start : start + BATCH]
        what = f"{study} samples {batch[0]}..{batch[-1]}"
        body = fetch(BROWSER_XML + ",".join(batch), what, opener, counter)
        for record in parse_sample_xml(body):
            if record.accession in samples:
                record.runs = samples[record.accession].runs
                record.layouts = samples[record.accession].layouts
                samples[record.accession] = record
                returned.add(record.accession)
    missing = [a for a in accessions if a not in returned]
    if missing:
        raise EnaError(f"{study}: ENA returned no record for {len(missing)} samples: {missing[:5]}")
    return list(samples.values())


def tag_order(samples: Iterable[Sample]) -> list[str]:
    """Every tag the samples carry, in first-seen order."""
    return list(dict.fromkeys(tag for s in samples for tag, _ in s.attributes))


def wide_rows(samples: list[Sample]) -> tuple[list[str], list[list[str]]]:
    """The header and one row per sample: fixed columns, then one per tag.

    A tag named as a fixed column (ArrayExpress samples carry ``scientific_name``) heads its
    column as ``attribute_<tag>``, so no header repeats; the long table keeps the raw tag.
    """
    tags = tag_order(samples)
    rows = []
    for s in samples:
        values: dict[str, list[str]] = {}
        for tag, value in s.attributes:
            values.setdefault(tag, []).append(value)
        fixed = [s.accession, s.secondary, s.alias, s.title, s.scientific_name, s.center_name,
                 s.broker_name, ";".join(s.runs), ";".join(s.layouts), s.description]  # fmt: skip
        rows.append(fixed + ["; ".join(values.get(tag, [])) for tag in tags])
    return [*FIXED, *(f"attribute_{t}" if t in FIXED else t for t in tags)], rows


def long_rows(samples: list[Sample]) -> list[list[str]]:
    """One row per sample, tag and value, in document order."""
    return [[s.accession, tag, value] for s in samples for tag, value in s.attributes]


def tsv_bytes(header: Iterable[str], rows: Iterable[list[str]]) -> bytes:
    """A tab-separated table, quoted only where csv's minimal quoting must."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return buffer.getvalue().encode()


def write_atomic(files: dict[Path, bytes]) -> None:
    """Write every file under a temporary name, then rename each; none is left half-written."""
    temporary = {path: path.with_name(f".{path.name}.part") for path in files}
    try:
        for path, body in files.items():
            temporary[path].write_bytes(body)
        for path, part in temporary.items():
            os.replace(part, path)
    finally:
        for part in temporary.values():
            part.unlink(missing_ok=True)


def output_paths(
    study: str, out_dir: Path, version: int | None, label: str | None = None
) -> tuple[Path, Path]:
    """The wide and long table paths, deliverable-named when a version is given.

    A label goes before the accession; a deliverable name holds the accession lowercased.
    """
    named = f"{label}_{study}" if label else study
    if version is None:
        return out_dir / f"{named}_samples.tsv", out_dir / f"{named}_attributes_long.tsv"
    stem = f"{date.today().isoformat()}_ena_samples_{named.lower()}"
    return out_dir / f"{stem}_v{version}.tsv", out_dir / f"{stem}_attributes_long_v{version}.tsv"


def summary(study: str, samples: list[Sample], counter: Requests) -> str:
    """One line: samples, runs, requests, and each tag with how many samples carry it."""
    carried = Counter(tag for s in samples for tag in dict.fromkeys(t for t, _ in s.attributes))
    tags = ", ".join(f"{tag} {carried[tag]}" for tag in tag_order(samples))
    return (f"{study}: {len(samples)} samples, {sum(len(s.runs) for s in samples)} runs, "
            f"{counter.made} requests ({counter.retried} retries); tags: {tags}")  # fmt: skip


def parse_labels(pairs: Iterable[str], studies: Iterable[str]) -> dict[str, str]:
    """``<study>=<label>`` pairs as a mapping, each label checked against the naming rule.

    Raises:
        EnaError: A pair lacks ``=``, names a study not asked for, repeats one, or its label
            is not lowercase words joined by underscores.
    """
    asked, labels = set(studies), {}
    for pair in pairs:
        study, sep, label = pair.partition("=")
        if not sep or study not in asked:
            raise EnaError(f"--label {pair!r}: expected <study>=<label> for a study given")
        if study in labels:
            raise EnaError(f"--label {pair!r}: {study} is labelled twice")
        if not LABEL.match(label):
            raise EnaError(f"--label {pair!r}: a label is lowercase words joined by underscores")
        labels[study] = label
    return labels


def run(
    study: str,
    out_dir: Path,
    version: int | None,
    opener: Opener = urlopen,
    label: str | None = None,
) -> str:
    """Fetch one study and write its two tables; returns its summary line."""
    counter = Requests()
    samples = fetch_study(study, opener, counter)
    header, rows = wide_rows(samples)
    wide, long = output_paths(study, out_dir, version, label)
    write_atomic({wide: tsv_bytes(header, rows), long: tsv_bytes(LONG, long_rows(samples))})
    return summary(study, samples, counter)


def main(argv: list[str] | None = None) -> int:
    """Command line: one or more study accessions, an output folder, an optional version."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("studies", nargs="+", help="PRJNA/PRJEB/SRP/ERP accessions")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--version", type=int, help="write deliverable names, _v<N>")
    parser.add_argument(
        "--label", action="append", default=[], metavar="STUDY=LABEL",
        help="put LABEL in STUDY's file names, e.g. PRJNA000001=lung_tissue",
    )  # fmt: skip
    args = parser.parse_args(argv)
    try:
        labels = parse_labels(args.label, args.studies)
        args.out_dir.mkdir(parents=True, exist_ok=True)
        for study in args.studies:
            print(run(study, args.out_dir, args.version, label=labels.get(study)), flush=True)
    except EnaError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
