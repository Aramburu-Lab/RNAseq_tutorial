# chr22 fixture: a GRCh38 slice and simulated reads of known origin

Written by `scripts/python/make_chr22_fixture.py`; never edit it by hand. To rebuild, run
`pixi run python scripts/python/make_chr22_fixture.py --out <empty directory>` and compare
`MANIFEST.sha256`, or check this copy with `shasum -a 256 -c MANIFEST.sha256`. The sequence is
the reference genome and the reads are simulated: no sample or study is involved.

## Source

| | |
|---|---|
| Annotation | GENCODE 50 basic, GRCh38.p14 |
| URL | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz |
| Publisher checksum | MD5 `11e77cf1d4b78b43f57d608e3f7f2202`, from `MD5SUMS` in the same directory |
| SHA-256 | `5b9ad9aa71f6ab5d55123c51f74723644a895bc0cbaa9cee8fbc11bb9dec0bf7` (first download, 2026-09-26; GENCODE publishes none) |
| Sequence | Ensembl 116, chromosome 22, unmasked (`dna`) |
| URL | https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz |
| Publisher checksum | BSD `sum` `31309 11123`, from `CHECKSUMS` in the same directory |
| SHA-256 | `974f97ac8ef7ffae971b63b47608feda327403be40c27e391ee4a1a78b800df5` (first download, 2026-09-26; Ensembl publishes none) |

GENCODE 50 is Ensembl 116
(https://www.gencodegenes.org/human/releases.html; the GTF header reads
`version 50 (Ensembl 116)`). Ensembl names the chromosome `22` and
GENCODE `chr22`; the coordinates are the same.

## Window

`chr22:30330001-30431000`, 1-based and closed, 101000 bp, written as the contig
`chr22_slice`. Offset 30330000: chr22 position = chr22_slice position + 30330000. Both edges fall
between genes, and every gene lying wholly inside is kept, whatever its biotype. It was picked
by a scan of chr22 for the narrowest run of whole genes, 2 kb clear of other genes at both
ends, holding at least two multi-exon protein-coding genes on each strand in at most 45
transcripts; the edges were then rounded.

| gene | gene_id | strand | type | slice start | slice end | transcripts | pairs |
|---|---|---|---|---|---|---|---|
| SF3A1 | ENSG00000099995.21 | - | protein_coding | 1988 | 26919 | 3 | 131 |
| CCDC157 | ENSG00000187860.13 | + | protein_coding | 26631 | 48673 | 8 | 506 |
| RNF215 | ENSG00000099999.17 | - | protein_coding | 38811 | 91771 | 6 | 311 |
| SEC14L2 | ENSG00000100003.18 | + | protein_coding | 66941 | 95303 | 4 | 385 |
| ENSG00000249590 | ENSG00000249590.7 | + | protein_coding | 79255 | 98990 | 1 | 114 |
| ENSG00000273350 | ENSG00000273350.1 | + | lncRNA | 90512 | 90912 | 1 | 40 |
| RNU6-564P | ENSG00000222915.1 | + | snRNA | 91045 | 91148 | 1 | 0 |
| ENSG00000272689 | ENSG00000272689.1 | - | lncRNA | 91206 | 91536 | 1 | 0 |
| ENSG00000303189 | ENSG00000303189.1 | - | lncRNA | 93830 | 94211 | 1 | 13 |
| MTFP1 | ENSG00000242114.8 | + | protein_coding | 95623 | 99058 | 10 | 1000 |

`transcripts.fa` holds the 36 transcripts spliced from `slice.fa` at the
coordinates in `slice.gtf`, minus-strand ones reverse-complemented.

Exonic bases shared by two genes, counted from `slice.gtf`. On opposite strands a stranded
count separates the genes and an unstranded one cannot; on the same strand neither can.

- SF3A1 (-) and CCDC157 (+): 165 bp
- CCDC157 (+) and RNF215 (-): 4 bp
- SEC14L2 (+) and ENSG00000249590 (+): 583 bp
- SEC14L2 (+) and ENSG00000303189 (-): 382 bp
- ENSG00000249590 (+) and MTFP1 (+): 887 bp

## Reads

- Seed 42, 2500 pairs, 2 x 100 nt, fragment length uniform on
  150-350 nt.
- Reverse-stranded (dUTP): read 1 is antisense, from the fragment's 3' end; read 2 is sense.
  Salmon library type `ISR`; featureCounts `-s 2`.
- No sequencing errors and no adapter (every fragment is longer than a read); quality `I`
  (Phred 40) on every base.
- Read name `<pair>:<transcript_id>:<start>:<length>/1` (or `/2`); start is 0-based on the
  transcript, 5' to 3'.
- Each transcript of at least 350 nt gets a fragment share from 1 to 64 (a squared
  uniform draw); the 2 shorter ones get no reads and stay in `truth.tsv` at 0.
- `truth.tsv`: `transcript_id`, `gene_id`, `n_pairs`.

## What it can test

- Parsing and pairing of FASTQ, FASTA and GTF, and their agreement with each other.
- Spliced alignment on both strands: every read's transcript and offset are in its name.
- Strandedness: the wrong strand setting loses reads, and the opposite-strand pairs above
  show what an unstranded count cannot separate.
- Transcript and gene quantification against exact counts, including zeros.
- Byte-identical regeneration from the seed.

## What it cannot test

- Sequencing errors, adapters, low-quality tails or trimming: none are planted.
- Intronic, intergenic, pre-mRNA, rRNA or genomic DNA reads: every read is from a mature
  transcript.
- Multi-mapping to paralogues or pseudogenes elsewhere: the reference is 101000 bp.
- GC, positional or 3' bias, PCR duplicates, a realistic fragment-length distribution.
- Anything between samples, such as differential expression: there is one sample.
- Novel transcripts, or transcripts shorter than 350 nt, which get no reads.

## Files

| file | bytes | SHA-256 |
|---|---|---|
| reads_R1.fq | 600319 | c0f346af62648d378afd4c877d138ec4395a9d5e7833dfe539926a67d462c751 |
| reads_R2.fq | 600319 | 970c5d203e224d4f467d00db2c411c90b52f72e97598300c39cc037e0da7bf2f |
| slice.fa | 102721 | c15340b9dc82d16e440441936a12c852e48117e629f2d75cff7aa96e149187ff |
| slice.gtf | 297410 | 05cb4efac5d81ff83204a9123a1a8ba56d40a149f4b9004062e25be15740adea |
| transcripts.fa | 74717 | a88233b41cc9c70d29709ddff9b8857e81680ebca62d54e66c11caf019c0422d |
| truth.tsv | 1460 | 867589549e3f3d3930adcbd79bce6fced7d106ca6856333119bfc95cfa5461e0 |
