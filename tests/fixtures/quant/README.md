# quant fixture: Salmon and R tximport on chr22_se

Written by `scripts/python/make_quant_fixture.py`; never edit it by hand. Salmon 2.7.0
(the pinned `reads` environment) quantified the five `tests/fixtures/chr22_se` samples with
`-l A` against `tests/fixtures/chr22/transcripts.fa`, its `slice.fa` as decoys, k 25. R
tximport 1.30.0 in `community.wave.seqera.io/library/bioconductor-tximeta_jq:b3a71c30effcaa9d` wrote the merged tables with nf-core/rnaseq's columns.

- `rnaseq/` is an nf-core/rnaseq outdir as `ctx quant` reads it. `pipeline_info/` is the
  yeast fixture's (`tests/fixtures/rnaseq`), so the release check passes; its params name
  the yeast genome, not this one. `salmon/` holds what copy Q1 holds per sample: `quant.sf`,
  `cmd_info.json`, `lib_format_counts.json`, `aux_info/meta_info.json`,
  `aux_info/ambig_info.tsv` and `libParams/flenDist.txt`, and the merged tables.
- `planted/` holds four made-up `quant.sf` files, their `tx2gene.tsv` and R's gene tables
  (`r.gene_*.tsv`): effective lengths differ by sample, G1 is zero in s3, G2 is zero
  everywhere and T9 is in no gene.

A rebuild repeats every count; `meta_info.json` records run times, so its bytes differ.
