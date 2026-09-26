"""The committed chr22 fixture is what its manifest says, self-consistent, and reproducible.

The parsers and the splice here are written independently of the builder script, so a bug
there cannot be copied into the check.
"""

import hashlib
import importlib.util
import random
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "chr22"
SCRIPT = ROOT / "scripts" / "python" / "make_chr22_fixture.py"
BUDGET = 2_000_000
READS = ("reads_R1.fq", "reads_R2.fq", "truth.tsv")
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
ATTR = re.compile(r'(\S+) "([^"]*)"')


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("make_chr22_fixture", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest():
    lines = (FIXTURE / "MANIFEST.sha256").read_text().splitlines()
    return {name: digest for digest, name in (line.split("  ", 1) for line in lines)}


def revcomp(sequence):
    return "".join(COMPLEMENT[base] for base in reversed(sequence))


def fasta(path):
    records, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        else:
            records[name].append(line)
    return {name: "".join(lines) for name, lines in records.items()}


def parsed(rows):
    return [(r[0], r[2], int(r[3]), int(r[4]), r[6], dict(ATTR.findall(r[8]))) for r in rows]


def gtf(path):
    lines = path.read_text().splitlines()
    return parsed(line.split("\t") for line in lines if not line.startswith("#"))


def independent_splice(sequence, rows, kind="exon"):
    """Concatenate one feature type in genomic order, reverse-complementing minus-strand joins."""
    exons = {}
    for _, feature, low, high, strand, attrs in rows:
        if feature == kind:
            exons.setdefault(attrs["transcript_id"], []).append((low, high, strand))
    out = {}
    for transcript, parts in exons.items():
        joined = "".join(sequence[low - 1 : high] for low, high, _ in sorted(parts))
        out[transcript] = revcomp(joined) if parts[0][2] == "-" else joined
    return out


def fastq(path):
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0
    return [
        (lines[i][1:], lines[i + 1], lines[i + 2], lines[i + 3]) for i in range(0, len(lines), 4)
    ]


def test_manifest():
    entries = manifest()
    on_disk = sorted(p.name for p in FIXTURE.iterdir() if p.name != "MANIFEST.sha256")
    assert len(entries) > 0
    assert sorted(entries) == on_disk
    assert {name: sha256(FIXTURE / name) for name in entries} == entries
    total = sum(p.stat().st_size for p in FIXTURE.iterdir())
    assert 0 < total <= BUDGET
    table = re.findall(
        r"^\| (\S+) \| (\d+) \| ([0-9a-f]{64}) \|$", (FIXTURE / "README.md").read_text(), re.M
    )
    assert {name: digest for name, _, digest in table} == {
        n: d for n, d in entries.items() if n != "README.md"
    }
    assert {name: int(size) for name, size, _ in table} == {
        n: (FIXTURE / n).stat().st_size for n, *_ in table
    }


def test_tracked_not_ignored():
    files = sorted(p.relative_to(ROOT).as_posix() for p in FIXTURE.iterdir())
    assert len(files) > 1
    tracked = subprocess.run(
        ["git", "ls-files", "--", *files], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.split()
    assert sorted(tracked) == files
    ignored = subprocess.run(
        ["git", "check-ignore", "--no-index", "--", *files],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert ignored.returncode == 1, f"ignored: {ignored.stdout.split()}"


def test_gtf_fasta_consistent(script):
    contigs = fasta(FIXTURE / "slice.fa")
    assert list(contigs) == [script.CONTIG]
    sequence = contigs[script.CONTIG]
    assert len(sequence) == script.WINDOW_END - script.WINDOW_START + 1
    assert set(sequence) <= set("ACGT")
    rows = gtf(FIXTURE / "slice.gtf")
    assert rows
    assert all(
        contig == script.CONTIG and 1 <= low <= high <= len(sequence)
        for contig, _, low, high, _, _ in rows
    )

    transcripts = fasta(FIXTURE / "transcripts.fa")
    assert independent_splice(sequence, rows) == transcripts
    # Anchor to biology, not to the builder: a wrong rebase offset keeps the three files
    # consistent with each other but moves every codon off its ATG and stop.
    starts = independent_splice(sequence, rows, "start_codon")
    stops = independent_splice(sequence, rows, "stop_codon")
    assert starts and stops
    assert set(starts.values()) == {"ATG"}
    assert set(stops.values()) <= {"TAA", "TAG", "TGA"}
    exon_attrs = {r[5]["transcript_id"]: r[5] for r in rows if r[1] == "exon"}
    for line in (FIXTURE / "transcripts.fa").read_text().splitlines():
        if line.startswith(">"):
            name, *fields = line[1:].split()
            header = dict(f.split("=", 1) for f in fields)
            source = exon_attrs[name]
            assert (header["gene_id"], header["gene_name"]) == (
                source["gene_id"],
                source["gene_name"],
            )
    exons_per = Counter((r[5]["transcript_id"], r[4]) for r in rows if r[1] == "exon")
    for strand in "+-":
        assert any(n > 1 for (_, s), n in exons_per.items() if s == strand), (
            f"no multi-exon {strand}"
        )

    # Order: rebase first, then extract from the slice. That must equal extracting from the
    # whole chromosome at genomic coordinates; extracting before rebasing must fail loudly.
    rng = random.Random(7)
    chromosome = "".join("ACGT"[int(rng.random() * 4)] for _ in range(3000))
    attrs = 'gene_id "{g}"; transcript_id "{t}"; gene_name "{g}";'
    genes = [
        ("g1", "+", [(1001, 1100), (1301, 1450), (1701, 1800)]),
        ("g2", "-", [(1200, 1260), (1500, 1650), (1850, 1900)]),
    ]
    genomic = []
    for gene, strand, exons in genes:
        genomic.append(row("gene", exons[0][0], exons[-1][1], strand, attrs.format(g=gene, t=gene)))
        genomic += [
            row("exon", lo, hi, strand, attrs.format(g=gene, t=f"{gene}.1")) for lo, hi in exons
        ]
    start, end = 901, 2000
    window = script.window_features(genomic, start, end)
    assert window == genomic
    sliced = chromosome[start - 1 : end]
    rebased = script.rebase(window, start - 1, "chrS_slice")
    truth = independent_splice(chromosome, parsed(genomic))
    built = {t.transcript_id: t.sequence for t in script.splice(sliced, rebased)}
    assert built == truth
    assert len(truth["g2.1"]) == 61 + 151 + 51
    with pytest.raises(ValueError, match="outside"):
        script.splice(sliced, window)
    crossing = row("gene", 850, 950, "+", 'gene_id "g0";')
    with pytest.raises(ValueError, match="crosses"):
        script.window_features([crossing, *genomic], start, end)


def row(feature, low, high, strand, attrs):
    return ["chrS", "t", feature, str(low), str(high), ".", strand, ".", attrs]


def test_reads_known_truth(script):
    transcripts = fasta(FIXTURE / "transcripts.fa")
    read1, read2 = fastq(FIXTURE / "reads_R1.fq"), fastq(FIXTURE / "reads_R2.fq")
    assert len(read1) == len(read2) == script.N_PAIRS
    counts = Counter()
    for (name1, seq1, plus1, qual1), (name2, seq2, plus2, qual2) in zip(read1, read2, strict=True):
        assert name1.endswith("/1") and name2.endswith("/2") and name1[:-2] == name2[:-2]
        assert plus1 == plus2 == "+"
        assert len(seq1) == len(seq2) == script.READ_LENGTH
        assert qual1 == qual2 == script.QUALITY * script.READ_LENGTH
        _, transcript, start, length = name1[:-2].split(":")
        start, length = int(start), int(length)
        assert script.FRAGMENT_MIN <= length <= script.FRAGMENT_MAX
        fragment = transcripts[transcript][start : start + length]
        assert len(fragment) == length
        # dUTP: read 1 antisense from the fragment's 3' end, read 2 sense from its 5' end.
        assert seq1 == revcomp(fragment[-script.READ_LENGTH :])
        assert seq2 == fragment[: script.READ_LENGTH]
        counts[transcript] += 1
    assert len({name[:-2] for name, *_ in read1}) == script.N_PAIRS

    header = {}
    for line in (FIXTURE / "transcripts.fa").read_text().splitlines():
        if line.startswith(">"):
            name, *fields = line[1:].split()
            header[name] = dict(f.split("=", 1) for f in fields)["gene_id"]
    lines = (FIXTURE / "truth.tsv").read_text().splitlines()
    assert lines[0] == "transcript_id\tgene_id\tn_pairs"
    truth = {t: (g, int(n)) for t, g, n in (line.split("\t") for line in lines[1:])}
    assert {t: g for t, (g, _) in truth.items()} == header
    assert {t: n for t, (_, n) in truth.items()} == {t: counts[t] for t in transcripts}
    assert sum(n for _, n in truth.values()) == script.N_PAIRS


def test_regenerate_byte_identical(script, tmp_path):
    source = tmp_path / "input" / "transcripts.fa"
    source.parent.mkdir()
    shutil.copyfile(FIXTURE / "transcripts.fa", source)
    out = tmp_path / "out"
    out.mkdir()
    script.reads(source, script.N_PAIRS, out_dir=out)  # the default seed is the committed one
    produced = {p.name: sha256(p) for p in out.iterdir()}
    assert sorted(produced) == sorted(READS)
    assert produced == {name: manifest()[name] for name in READS}
    other = tmp_path / "other_seed"
    script.reads(source, script.N_PAIRS, seed=script.SEED + 1, out_dir=other)
    assert sha256(other / "reads_R1.fq") != produced["reads_R1.fq"]


def test_refuses_nonempty_dir(script, tmp_path):
    (tmp_path / "keep.txt").write_text("x")
    with pytest.raises(FileExistsError):
        script.reads(FIXTURE / "transcripts.fa", 10, out_dir=tmp_path)
    cache = tmp_path.parent / f"{tmp_path.name}_cache"
    with pytest.raises(FileExistsError):
        script.main(["--out", str(tmp_path), "--cache", str(cache)])
    assert [p.name for p in tmp_path.iterdir()] == ["keep.txt"]
    assert not cache.exists(), "refusal must come before any download"


@pytest.mark.network
def test_slice_rebuilds(script, tmp_path):
    out, cache = tmp_path / "chr22", tmp_path / "cache"
    assert script.main(["--out", str(out), "--cache", str(cache)]) == 0
    assert (out / "MANIFEST.sha256").read_bytes() == (FIXTURE / "MANIFEST.sha256").read_bytes()
    # Independent of the rebase: splice the whole chromosome at genomic coordinates.
    chromosome = script.read_chromosome(cache / script.FASTA.url.rsplit("/", 1)[1])
    rows = script.gtf_rows(cache / script.GTF.url.rsplit("/", 1)[1], script.CHROM)
    genomic = script.window_features(rows, script.WINDOW_START, script.WINDOW_END)
    assert independent_splice(chromosome, parsed(genomic)) == fasta(FIXTURE / "transcripts.fa")
