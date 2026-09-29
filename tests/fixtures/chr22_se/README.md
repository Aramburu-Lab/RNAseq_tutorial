# chr22_se fixture: single-end reads of known origin

Written by `scripts/python/make_chr22_se_fixture.py`; never edit it by hand. To rebuild, run
`pixi run python scripts/python/make_chr22_se_fixture.py --out <empty directory>` and compare
`MANIFEST.sha256`. The reads are simulated from `tests/fixtures/chr22/transcripts.fa`, checked
against that fixture's manifest first: no sample or study is involved.

Step-1 is single-end, 75-101 nt; the chr22 fixture's reads are paired-end. `samples.tsv` says how
each of the five samples was made: 1000 reads each, seed 42. Four are reverse-stranded
(`library_type` SR: Salmon SR, featureCounts -s 2, STAR GeneCounts column 4) and one is
unstranded (U). Two expression profiles are crossed with two read lengths, so samples cluster
by profile, not by read length.

Planted for `ctx qc`: `se_c`'s fragments are 60-200 nt against 101 nt reads, so reads from shorter
fragments run into the adapter `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` and then poly-A; `se_d`'s last 30 qualities fall from
Phred 30 to 2, its bases unchanged, so the truth still holds.

Read names are `<n>:<transcript>:<0-based start>:<fragment length>:<as|s>`, the strand of the
read against the transcript. `truth.tsv` counts reads per sample and transcript, and how many ran
into the adapter.
