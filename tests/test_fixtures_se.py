"""tests/fixtures/chr22_se: single-end reads of known origin, checked read by read."""

import hashlib
import importlib.util
import shutil
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "chr22_se"
SCRIPT = ROOT / "scripts" / "python" / "make_chr22_se_fixture.py"
BUDGET = 1_500_000
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
# What a read runs into past a short fragment: the TruSeq adapter, then poly-A.
TAIL = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" + "A" * 101


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("make_chr22_se_fixture", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def transcripts() -> dict[str, str]:
    """transcripts.fa of the chr22 fixture, parsed here without the generator's parser."""
    records: dict[str, list[str]] = {}
    name = ""
    source = ROOT / "tests" / "fixtures" / "chr22" / "transcripts.fa"
    for line in source.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        else:
            records[name].append(line)
    return {k: "".join(v) for k, v in records.items()}


def reads(sample: str):
    lines = (FIXTURE / f"{sample}.fq").read_text().splitlines()
    for i in range(0, len(lines), 4):
        yield lines[i][1:], lines[i + 1], lines[i + 3]


def samples() -> list[dict[str, str]]:
    rows = (FIXTURE / "samples.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    return [dict(zip(header, r.split("\t"), strict=True)) for r in rows[1:]]


def test_manifest_and_budget():
    listed = dict(
        reversed(line.split("  "))
        for line in (FIXTURE / "MANIFEST.sha256").read_text().splitlines()
    )
    on_disk = {p.name for p in FIXTURE.iterdir()} - {"MANIFEST.sha256"}
    assert set(listed) == on_disk
    for name, digest in listed.items():
        assert hashlib.sha256((FIXTURE / name).read_bytes()).hexdigest() == digest, name
    assert sum(p.stat().st_size for p in FIXTURE.iterdir()) < BUDGET


def test_every_read_is_where_its_name_says():
    """Orientation, adapter read-through and truth counts, from the names alone."""
    source = transcripts()
    truth = Counter()
    adapters = Counter()
    for sample in samples():
        length = int(sample["read_length"])
        strands = Counter()
        for name, read, _ in reads(sample["sample"]):
            _, tx, start, size, strand = name.split(":")
            fragment = source[tx][int(start) : int(start) + int(size)]
            insert = (
                "".join(COMPLEMENT[b] for b in reversed(fragment)) if strand == "as" else fragment
            )
            assert read.startswith(insert[:length])
            assert (int(size) < length) == (len(insert) < length)
            if int(size) < length:
                assert read[int(size) :] == TAIL[: length - int(size)]
                adapters[(sample["sample"], tx)] += 1
            truth[(sample["sample"], tx)] += 1
            strands[strand] += 1
        assert set(strands) == ({"as"} if sample["library_type"] == "SR" else {"as", "s"})
    rows = (FIXTURE / "truth.tsv").read_text().splitlines()[1:]
    for row in rows:
        sample, tx, _, n_reads, n_adapter = row.split("\t")
        assert truth[(sample, tx)] == int(n_reads)
        assert adapters[(sample, tx)] == int(n_adapter)


def test_planted_problems_are_there():
    # Reads from fragments shorter than the read run into the adapter: 276 of 1,000 in se_c.
    into = sum(1 for name, _, _ in reads("se_c") if int(name.split(":")[3]) < 101)  # noqa: PLR2004
    assert into == 276  # noqa: PLR2004 - measured when the fixture was built
    for sample in ("se_a", "se_b", "se_d", "se_e"):
        assert not any("AGATCGGAAGAGC" in read for _, read, _ in reads(sample))
    (_, _, quality) = next(reads("se_d"))
    phred = [ord(c) - 33 for c in quality[-30:]]
    assert phred[0] == 30 and phred[-1] == 2 and phred == sorted(phred, reverse=True)  # noqa: PLR2004
    assert set(next(reads("se_a"))[2]) == {"I"}


def test_rebuilds_byte_for_byte(script, tmp_path):
    script.main(["--out", str(tmp_path / "out")])
    for path in FIXTURE.iterdir():
        assert (tmp_path / "out" / path.name).read_bytes() == path.read_bytes(), path.name


def test_refuses_a_changed_source(script, tmp_path, monkeypatch):
    copy = tmp_path / "chr22"
    shutil.copytree(ROOT / "tests" / "fixtures" / "chr22", copy)
    (copy / "transcripts.fa").write_text((copy / "transcripts.fa").read_text() + "\n")
    monkeypatch.setattr(script, "SOURCE", copy)
    with pytest.raises(ValueError, match="does not match its manifest"):
        script.main(["--out", str(tmp_path / "out")])


def test_refuses_a_non_empty_directory(script, tmp_path):
    (tmp_path / "old.fq").write_text("x")
    with pytest.raises(FileExistsError):
        script.main(["--out", str(tmp_path)])
