"""ctx.process.reference on the chr22 fixture slice: genome, GTF subset, start codons, indices."""

import gzip
import re
import shutil
from pathlib import Path

import pytest

from ctx.process import STRAND, ProcessError, Strandedness, reference
from ctx.qc import registry

ROOT = Path(__file__).resolve().parents[1]
SLICE = ROOT / "tests" / "fixtures" / "chr22" / "slice.fa"
GTF = ROOT / "tests" / "fixtures" / "chr22" / "slice.gtf"
REG = registry.load(ROOT / "conf" / "qc_metrics.toml")


def gene_names(gtf: Path) -> set[str]:
    """The gene_name of every gene record, read here without the module."""
    return {re.search(r'gene_name "([^"]+)"', line)[1]
            for line in gtf.read_text().splitlines() if "\tgene\t" in line}  # fmt: skip


def test_one_strand_map():
    """STAR column and featureCounts -s per strand (STAR manual section 8; Subread guide 6.2.10).

    Mutation caught: columns 3 and 4 swapped."""
    assert STRAND == {Strandedness.UNSTRANDED: (2, 0), Strandedness.FORWARD: (3, 1),
                      Strandedness.REVERSE: (4, 2)}  # fmt: skip
    assert set(Strandedness) - set(STRAND) == {Strandedness.AUTO}


@pytest.mark.parametrize(
    ("length", "bases"), [(101_000, 7), (50_818_468, 11), (39_159_782, 11), (2**32, 14)]
)
def test_sa_index_bases(length, bases):
    """min(14, floor(log2(L)/2 - 1)): the slice, chr22 with and without its N gaps, the cap.

    Mutations caught: round (12 for chr22); ceil (8 for the slice); no cap (15 at 2^32)."""
    assert reference.sa_index_bases(length) == bases


def test_the_subset_matches_contigs_exactly(tmp_path):
    """Mutation caught: substring matching (a "chr2" genome would take chr22_slice's lines)."""
    with pytest.raises(ProcessError, match=r"share no contig.*\['chr2'\].*\['chr22_slice'\]"):
        reference.annotation_subset(GTF, {"chr2"}, tmp_path / "a.gtf")


def test_the_slice_subset_keeps_its_header_and_genes(tmp_path):
    """Mutation caught: the ## header lines dropped."""
    subset = reference.annotation_subset(GTF, {"chr22_slice"}, tmp_path / "a.gtf")
    written = (tmp_path / "a.gtf").read_text().splitlines()
    source = GTF.read_text().splitlines()
    assert written[:2] == source[:2] and all(line.startswith("##") for line in written[:2])
    assert (subset.genes, subset.lines) == (10, len(source) - 2)
    assert subset.release == "GENCODE 50 (Ensembl 116)"
    assert subset.sha256 == reference.sha256(tmp_path / "a.gtf")


def test_start_codons_read_atg_on_both_strands(tmp_path):
    """31 of 31 on the slice, minus-strand ones among them; one base of shift reads none.

    Mutation caught: no reverse complement on the minus strand."""
    minus = [line for line in GTF.read_text().splitlines()
             if "\tstart_codon\t" in line and line.split("\t")[6] == "-"]  # fmt: skip
    assert minus
    codons = reference.start_codon_share(SLICE, GTF)
    assert (codons.records, codons.atg, codons.split) == (31, 31, 0)
    header, *lines = SLICE.read_text().splitlines()
    (tmp_path / "shifted.fa").write_text(f"{header}\nN{''.join(lines)}\n")
    shifted = reference.start_codon_share(tmp_path / "shifted.fa", GTF)
    assert (shifted.records, shifted.atg) == (31, 0)
    metric = REG.metrics["reference_start_codons"]
    assert (metric.band(codons.share), metric.band(shifted.share)) == ("good", "bad")


def test_soft_masked_and_split_codons(tmp_path):
    """A soft-masked codon still reads ATG; a split record is counted, never read.

    Mutations caught: no upper(); split records read as whole ones."""
    (tmp_path / "g.fa").write_text(">c1\nggatgcccatgg\n")
    rows = [(3, 5, "+"), (8, 10, "-"), (11, 12, "+")]
    lines = [f'c1\tt\tstart_codon\t{a}\t{b}\t.\t{s}\t0\tgene_id "g";\n' for a, b, s in rows]
    (tmp_path / "g.gtf").write_text("".join(lines))
    codons = reference.start_codon_share(tmp_path / "g.fa", tmp_path / "g.gtf")
    assert codons == reference.StartCodons(2, 2, 1, {"ATG": 2})


def test_a_repeated_contig_name_is_refused(tmp_path):
    """Two records named alike would hand STAR one contig's annotation for the other."""
    (tmp_path / "g.fa").write_text(">c1\nACGT\n>c1 again\nACGT\n")
    with pytest.raises(ProcessError, match="repeated"):
        list(reference.records(tmp_path / "g.fa"))


def test_a_gzipped_genome_is_copied_uncompressed(tmp_path):
    """STAR refuses a gzipped FASTA; the hash recorded is of the file as given.

    Mutation caught: hashing the uncompressed text instead of the file."""
    given = tmp_path / "slice.fa.gz"
    with SLICE.open("rb") as src, gzip.open(given, "wb") as dst:
        shutil.copyfileobj(src, dst)
    ref = reference.prepare(given, GTF, tmp_path / "index")
    assert ref.source.sha256 == reference.sha256(given) != reference.sha256(SLICE)
    assert ref.source.bytes == given.stat().st_size
    assert list(reference.records(ref.fasta)) == list(reference.records(SLICE))
    assert ref.fasta.read_bytes()[:1] == b">"
    assert (ref.genome.lengths, ref.genome.n_bases) == ((101_000,), 0)


def test_the_cut_keeps_a_gene_by_its_end(tmp_path):
    """Half the slice, 50,500 bp, keeps exactly SF3A1 and CCDC157, every line of them.

    Mutation caught: keeping a gene by its start (RNF215 starts at 38,811 and ends at 91,771)."""
    ref = reference.prepare(SLICE, GTF, tmp_path, half=True)
    assert ref.genome.lengths == (50_500,)
    assert gene_names(ref.gtf) == {"SF3A1", "CCDC157"}
    kept = set(re.findall(r'gene_id "([^"]+)"', ref.gtf.read_text()))
    assert len(kept) == ref.annotation.genes == 2
    assert ref.start_codons.share == 1.0


def test_a_feature_before_its_gene_is_refused_when_cutting(tmp_path):
    """The cut decides by the gene record, so a transcript before it cannot be placed.

    Mutation caught: such a line kept by default."""
    lines = GTF.read_text().splitlines(keepends=True)
    first_gene = next(i for i, line in enumerate(lines) if "\tgene\t" in line)
    lines.insert(first_gene + 2, lines.pop(first_gene))
    (tmp_path / "late.gtf").write_text("".join(lines))
    with pytest.raises(ProcessError, match="before its gene record"):
        reference.annotation_subset(tmp_path / "late.gtf", {"chr22_slice"}, tmp_path / "a.gtf",
                                    cut={"chr22_slice": 50_500})  # fmt: skip


def test_each_input_keys_its_own_folders(tmp_path):
    """A half reference, a changed annotation and a changed genome each get their own folders;
    HISAT2's plain index has no annotation, so a changed one reuses it.

    Mutations caught: half left out of prepare's key (the half call gets the whole genome); the
    annotation left out of STAR's key (old junctions reused); the genome left out of HISAT2's."""
    full = reference.prepare(SLICE, GTF, tmp_path)
    half = reference.prepare(SLICE, GTF, tmp_path, half=True)
    assert half.folder != full.folder
    assert (full.genome.lengths, half.genome.lengths) == ((101_000,), (50_500,))
    star = reference.star_index(full, tmp_path, 2)[0]
    hisat2 = reference.hisat2_index(full, tmp_path, 2)[0]
    renamed = tmp_path / "renamed.gtf"
    renamed.write_text(GTF.read_text().replace('gene_name "CCDC157"', 'gene_name "CCDC157x"'))
    annotated = reference.prepare(SLICE, renamed, tmp_path)
    assert annotated.genome.sha256 == full.genome.sha256
    assert annotated.annotation.sha256 != full.annotation.sha256
    assert reference.star_index(annotated, tmp_path, 2)[0] != star
    assert reference.hisat2_index(annotated, tmp_path, 2) == (hisat2, False)
    header, *lines = SLICE.read_text().splitlines()
    assert lines[500][30] == "C"
    lines[500] = lines[500][:30] + "A" + lines[500][31:]
    (tmp_path / "edited.fa").write_text("\n".join([header, *lines]) + "\n")
    edited = reference.prepare(tmp_path / "edited.fa", GTF, tmp_path)
    assert edited.genome.sha256 != full.genome.sha256
    assert edited.annotation.sha256 == full.annotation.sha256
    assert reference.star_index(edited, tmp_path, 2)[0] != star
    assert reference.hisat2_index(edited, tmp_path, 2)[0] != hisat2


def test_relative_folders_still_reach_the_tools(tmp_path, monkeypatch):
    """Each index builds inside its own folder, where a path relative to the caller misses.

    Mutation caught: keyed's folder left as given (STAR exits 104, hisat2-build 1)."""
    monkeypatch.chdir(tmp_path)
    ref = reference.prepare(SLICE, GTF, Path("index"))
    star, _ = reference.star_index(ref, Path("index"), 2)
    hisat2, _ = reference.hisat2_index(ref, Path("index"), 2)
    assert {ref.folder.parent, star.parent, hisat2.parent} == {tmp_path / "index"}
    assert (star / "SA").is_file() and (hisat2 / "genome.1.ht2").is_file()
    assert not list((tmp_path / "index").glob("*.stale*"))


def test_both_indices_build_and_are_reused(tmp_path):
    """The pinned STAR builds indices on this Mac; it is mapping that fails there (D-036).

    Mutation caught: --genomeSAindexNbases not passed (STAR's default 14 on a 101 kb genome)."""
    ref = reference.prepare(SLICE, GTF, tmp_path)
    star, built = reference.star_index(ref, tmp_path, 2)
    assert built
    assert (star / "chrName.txt").read_text().split() == ["chr22_slice"]
    assert (star / "chrLength.txt").read_text().split() == ["101000"]
    assert (star / "geneInfo.tab").read_text().split("\n", 1)[0] == "10"
    params = (star / "genomeParameters.txt").read_text()
    assert re.search(r"^genomeSAindexNbases\t7$", params, re.M)
    assert re.search(r"^sjdbOverhang\t100$", params, re.M)
    hisat2, built = reference.hisat2_index(ref, tmp_path, 2)
    assert built
    assert sorted(p.name for p in hisat2.glob("genome.*.ht2")) == [
        f"genome.{n}.ht2" for n in range(1, 9)
    ]  # fmt: skip
    assert reference.prepare(SLICE, GTF, tmp_path).built is False
    assert reference.star_index(ref, tmp_path, 2) == (star, False)
    assert reference.hisat2_index(ref, tmp_path, 2) == (hisat2, False)
