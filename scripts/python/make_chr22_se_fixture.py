#!/usr/bin/env python3
# /// codebase-meta
# exempt: instructor-side builder of tests/fixtures/chr22_se; standard library only, offline
# ///
"""Build tests/fixtures/chr22_se: single-end reads of known origin from the chr22 fixture.

Step-1 is single-end, 75-101 nt; the chr22 fixture's reads are paired-end dUTP. This writes
five single-end samples simulated from its committed transcripts.fa, checked first against the
SHA-256 in that fixture's manifest: four reverse-stranded (Salmon SR, featureCounts -s 2, STAR
GeneCounts column 4) and one unstranded, two expression profiles crossed with two read
lengths, and two planted problems for ctx qc: an adapter-heavy sample and a degraded one.
Offline and deterministic: the transcript parser, the seed rule and the helpers are the chr22
generator's own, imported.

Usage:
    pixi run python scripts/python/make_chr22_se_fixture.py --out tests/fixtures/chr22_se
"""

from __future__ import annotations

import argparse
import importlib.util
import random
import sys
from bisect import bisect_right
from itertools import accumulate
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "tests" / "fixtures" / "chr22"
_spec = importlib.util.spec_from_file_location("make_chr22_fixture", Path(__file__).with_name(
    "make_chr22_fixture.py"))  # fmt: skip
chr22 = importlib.util.module_from_spec(_spec)  # type: ignore[arg-type]
_spec.loader.exec_module(chr22)  # type: ignore[union-attr]

SEED = 42
N_READS = 1000
# The Illumina TruSeq read 1 adapter as a read runs into it, then poly-A where the read is
# longer still, as a sequencer reads past the library.
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
QUALITY = "I"
# Phred scores of the degraded sample's last bases: falling from 30 to 2.
TAIL, TAIL_FROM, TAIL_TO = 30, 30, 2
MANIFEST = "MANIFEST.sha256"


class Sample(NamedTuple):
    """One simulated sample: how its reads are made, and the problem planted in it."""

    name: str
    read_length: int
    fragments: tuple[int, int]
    profile: int
    stranded: bool
    planted: str


SAMPLES = (
    Sample("se_a", 76, (150, 350), 1, True, "none"),
    Sample("se_b", 76, (150, 350), 2, True, "none"),
    Sample("se_c", 101, (60, 200), 1, True, "adapter: fragments shorter than the read"),
    Sample("se_d", 101, (150, 350), 2, True, "degraded: last 30 qualities fall from 30 to 2"),
    Sample("se_e", 76, (150, 350), 1, False, "none; unstranded"),
)


def quality(sample: Sample) -> str:
    """The quality string every read of a sample carries."""
    if not sample.planted.startswith("degraded"):
        return QUALITY * sample.read_length
    steps = [TAIL_FROM - (TAIL_FROM - TAIL_TO) * i // (TAIL - 1) for i in range(TAIL)]
    return QUALITY * (sample.read_length - TAIL) + "".join(chr(33 + q) for q in steps)


def shares(profile: int, n: int) -> list[float]:
    """Expression shares of one profile, 1 to 64, from a squared uniform as the chr22 reads."""
    rng = random.Random(SEED + profile)
    return [(1 + 7 * u) * (1 + 7 * u) for u in (rng.random() for _ in range(n))]


def simulate(sample: Sample, index: int, records: list, out: Path) -> dict[str, list[int]]:
    """Write one sample's FASTQ; return [reads, reads into adapter] per transcript."""
    low, high = sample.fragments
    eligible = [r for r in records if len(r.sequence) >= high]
    cumulative = list(accumulate(shares(sample.profile, len(eligible))))
    rng = random.Random(SEED + 10 + index)
    counts = {r.transcript_id: [0, 0] for r in records}
    qualities = quality(sample)
    with (out / f"{sample.name}.fq").open("w", newline="\n") as handle:
        for n in range(1, N_READS + 1):
            origin = eligible[bisect_right(cumulative, rng.random() * cumulative[-1])]
            length = low + int(rng.random() * (high - low + 1))
            start = int(rng.random() * (len(origin.sequence) - length + 1))
            fragment = origin.sequence[start : start + length]
            # Reverse-stranded: the read is antisense, from the fragment's 3' end. Unstranded:
            # either strand, half the time each.
            antisense = sample.stranded or rng.random() < 0.5  # noqa: PLR2004 - a fair coin
            insert = chr22.revcomp(fragment) if antisense else fragment
            into_adapter = length < sample.read_length
            read = (insert + ADAPTER + "A" * sample.read_length)[: sample.read_length]
            strand = "as" if antisense else "s"
            handle.write(f"@{n:06d}:{origin.transcript_id}:{start}:{length}:{strand}\n"
                         f"{read}\n+\n{qualities}\n")  # fmt: skip
            counts[origin.transcript_id][0] += 1
            counts[origin.transcript_id][1] += into_adapter
    return counts


def write(out: Path) -> None:
    """Check the source, simulate every sample, and write truth, sample sheet and manifest."""
    chr22.require_empty(out)
    source = SOURCE / "transcripts.fa"
    pinned = dict(
        reversed(line.split("  ")) for line in (SOURCE / MANIFEST).read_text().splitlines()
    )
    if chr22.digest(source, "sha256") != pinned["transcripts.fa"]:
        raise ValueError(f"{source} does not match its manifest; rebuild tests/fixtures/chr22")
    records = chr22.read_transcripts(source)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, sample in enumerate(SAMPLES):
        counts = simulate(sample, index, records, out)
        rows += [
            f"{sample.name}\t{r.transcript_id}\t{r.gene_id}\t{counts[r.transcript_id][0]}\t"
            f"{counts[r.transcript_id][1]}\n"
            for r in records
        ]
    (out / "truth.tsv").write_text(
        "sample\ttranscript_id\tgene_id\tn_reads\tn_adapter\n" + "".join(rows), newline="\n"
    )
    library = {True: "SR", False: "U"}
    (out / "samples.tsv").write_text(
        "sample\tread_length\tfragment_min\tfragment_max\tprofile\tlibrary_type\tplanted\n"
        + "".join(
            f"{s.name}\t{s.read_length}\t{s.fragments[0]}\t{s.fragments[1]}\t{s.profile}\t"
            f"{library[s.stranded]}\t{s.planted}\n"
            for s in SAMPLES
        ),
        newline="\n",
    )
    readme = README.format(n=N_READS, seed=SEED, adapter=ADAPTER)
    (out / "README.md").write_text(readme, newline="\n")
    manifest = "".join(
        f"{chr22.digest(p, 'sha256')}  {p.name}\n"
        for p in sorted(out.iterdir())
        if p.name != MANIFEST
    )
    (out / MANIFEST).write_text(manifest, newline="\n")


README = """# chr22_se fixture: single-end reads of known origin

Written by `scripts/python/make_chr22_se_fixture.py`; never edit it by hand. To rebuild, run
`pixi run python scripts/python/make_chr22_se_fixture.py --out <empty directory>` and compare
`MANIFEST.sha256`. The reads are simulated from `tests/fixtures/chr22/transcripts.fa`, checked
against that fixture's manifest first: no sample or study is involved.

Step-1 is single-end, 75-101 nt; the chr22 fixture's reads are paired-end. `samples.tsv` says how
each of the five samples was made: {n} reads each, seed {seed}. Four are reverse-stranded
(`library_type` SR: Salmon SR, featureCounts -s 2, STAR GeneCounts column 4) and one is
unstranded (U). Two expression profiles are crossed with two read lengths, so samples cluster
by profile, not by read length.

Planted for `ctx qc`: `se_c`'s fragments are 60-200 nt against 101 nt reads, so reads from shorter
fragments run into the adapter `{adapter}` and then poly-A; `se_d`'s last 30 qualities fall from
Phred 30 to 2, its bases unchanged, so the truth still holds.

Read names are `<n>:<transcript>:<0-based start>:<fragment length>:<as|s>`, the strand of the
read against the transcript. `truth.tsv` counts reads per sample and transcript, and how many ran
into the adapter.
"""


def main(argv: list[str] | None = None) -> int:
    """Build the fixture into an empty directory.

    Args:
        argv: Argument list without the program name; ``sys.argv[1:]`` when None.

    Returns:
        0 on success.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True, type=Path, help="Empty fixture directory.")
    write(parser.parse_args(argv).out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
