"""A gene's structure: MANE Select exons, Ensembl REST's MANE record, the qPCR exon pair.

MANE exons come from the GENCODE 50 GTF, the qPCR pair from GTEx's exon shares.

The GTF the dataset file registers is the counts' own annotation (Ensembl 116); REST is asked
for the same release, and adds the RefSeq transcript a primer tool takes. The two must name
the same MANE transcript and exon count, or the card says they differ.
"""

import gzip
import json
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from ctx.genecard import GenecardError

ENSEMBL_RELEASE = "116"
REST = "https://rest.ensembl.org"
PRIMER_BLAST = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/"
_GENE = re.compile(r'gene_id "ENSG([0-9]{11})\.[0-9]+"')
_TRANSCRIPT = re.compile(r'transcript_id "([^"]+)"')
_EXON_NUMBER = re.compile(r"exon_number ([0-9]+)")
MANE_TAG = 'tag "MANE_Select"'


def lookup_url(gene: str) -> str:
    """Ensembl REST's lookup of one gene, its transcripts and MANE records expanded."""
    return f"{REST}/lookup/id/{gene}?expand=1;mane=1;content-type=application/json"


def release_url() -> str:
    """Ensembl REST's own release record."""
    return f"{REST}/info/data?content-type=application/json"


def rest_release(body: bytes) -> str:
    """The release ``info/data`` reports, as text."""
    releases = json.loads(body).get("releases") or []
    return ",".join(str(r) for r in releases)


@dataclass(frozen=True)
class Mane:
    """A gene's MANE Select transcript and its exons, in transcription order.

    Attributes:
        transcript: The Ensembl transcript ID with its version.
        refseq: The matching RefSeq transcript, or "" where the source gives none.
        chrom: Chromosome.
        strand: "+" or "-".
        exons: (start, end) of each exon, 1-based and inclusive, in transcription order.
    """

    transcript: str
    refseq: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]

    def introns(self) -> tuple[tuple[int, int], ...]:
        """(start, end) of each intron between consecutive exons, 1-based and inclusive."""
        ordered = sorted(self.exons)
        return tuple((a[1] + 1, b[0] - 1) for a, b in zip(ordered, ordered[1:], strict=False))


def gtf_mane(path: Path, numbers: Iterable[int]) -> dict[int, Mane]:
    """The MANE Select exons of the wanted genes, streamed from a gzipped GTF.

    Raises:
        GenecardError: If a gene has exons of two MANE Select transcripts.
    """
    wanted = set(numbers)
    found: dict[int, tuple[str, str, str, list[tuple[int, int, int]]]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if MANE_TAG not in line or "\texon\t" not in line:
                continue
            fields = line.split("\t", 8)
            match = _GENE.search(fields[8])
            if match is None or (n := int(match[1])) not in wanted:
                continue
            transcript = _TRANSCRIPT.search(fields[8])[1]  # type: ignore[index]
            number = int(_EXON_NUMBER.search(fields[8])[1])  # type: ignore[index]
            held = found.setdefault(n, (transcript, fields[0], fields[6], []))
            if held[0] != transcript:
                raise GenecardError(f"ENSG{n:011d}: two MANE Select transcripts in {path.name}")
            held[3].append((number, int(fields[3]), int(fields[4])))
    return {n: Mane(t, "", c, s, tuple((a, b) for _, a, b in sorted(e)))
            for n, (t, c, s, e) in found.items()}  # fmt: skip


def rest_mane(body: bytes) -> Mane | None:
    """The MANE Select transcript of an Ensembl REST lookup, or None where it has none."""
    gene = json.loads(body)
    for transcript in gene.get("Transcript") or []:
        for record in transcript.get("MANE") or []:
            if record.get("type") == "MANE_Select":
                exons = sorted((e["start"], e["end"]) for e in transcript["Exon"])
                if gene.get("strand") == -1:
                    exons.reverse()
                return Mane(f"{record['id']}.{record['version']}", record.get("refseq_match", ""),
                            str(gene.get("seq_region_name", "")),
                            "+" if gene.get("strand") == 1 else "-", tuple(exons))  # fmt: skip
    return None


@dataclass(frozen=True)
class Pair:
    """The adjacent exon pair whose larger CV across tissues is smallest."""

    first: str
    second: str
    cv: float


def exon_cv(shares: pl.DataFrame) -> dict[str, float]:
    """Each exon's coefficient of variation of median share across tissues (SD over mean).

    An exon with no reads in any tissue has no CV.
    """
    out = {}
    for (exon,), part in shares.group_by("exon", maintain_order=True):
        values = part["share"].to_numpy()
        mean = values.mean()
        out[str(exon)] = float(values.std(ddof=1) / mean) if mean > 0 else float("nan")
    return out


def qpcr_pair(cv: dict[str, float]) -> Pair | None:
    """The preview's rule on GTEx's collapsed exons: neighbours, smallest larger CV.

    Exons are taken in their GTEx numbering, so neighbours in the collapsed model; a pair
    with an exon of no reads is never chosen.
    """
    names = sorted(cv, key=lambda e: int(e[1:]))
    best = None
    for a, b in zip(names, names[1:], strict=False):
        # Both checked: max(0.05, nan) is 0.05 in Python, so a NaN partner would pass.
        if not (np.isfinite(cv[a]) and np.isfinite(cv[b])):
            continue
        worst = max(cv[a], cv[b])
        if best is None or worst < best.cv:
            best = Pair(a, b, worst)
    return best
