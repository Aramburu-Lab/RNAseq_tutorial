# nf-core/rnaseq fixture: STAR and Salmon

A real nf-core/rnaseq run on its test profile with `aligner: star_salmon`, copied as written, never edited. Built by `scripts/python/make_rnaseq_fixture.py --aligned`.

| Field | Value |
| --- | --- |
| Command | `nextflow run nf-core/rnaseq -r 3.27.0 -profile test,docker,arm64 -params-file params.yml -c site.config, then again with -resume after skip_preseq: false was added to params.yml (.nextflow.log, 2026-10-23 11:15 and 11:41); site.config is publish.config` |
| Publish override | `publish.config`: TRIMGALORE also publishes `*_trimming_report.json` (D-025) |
| nf-core/rnaseq | v3.27.0-ga1fcddd |
| Nextflow | 26.04.6 |
| star | 2.7.11b |
| salmon | 2.7.0 |
| samtools | 1.21, 1.24 |
| rseqc | 5.0.4 |
| qualimap | 2.3 |
| preseq | 3.2.0 |
| bioconductor-dupradar | 1.38.0 |
| subread | 2.1.1 |
| trimgalore | 2.3.0 |
| fastqc | 0.12.1 |
| MultiQC (multiqc.log) | 1.35 |
| Profile | test,docker,arm64 |
| Host | Darwin 27.0.0 arm64 |
| Copied (UTC) | 2026-10-23 |

## Launches, as their params files record them

| Params file | Keys |
| --- | --- |
| `params_2026-10-23_11-15-51.json` | `aligner` = "star_salmon", `pseudo_aligner` = "salmon", `skip_alignment` = false, `skip_preseq` = true, `publish_dir_mode` = "copy" |
| `params_2026-10-23_11-42-10.json` | `aligner` = "star_salmon", `pseudo_aligner` = "salmon", `skip_alignment` = false, `skip_preseq` = false, `publish_dir_mode` = "copy" |

## Left out, per folder

| Folder | Files | Bytes |
| --- | --- | --- |
| `bbsplit/` | 5 | 1108 |
| `fastqc/` | 48 | 21096444 |
| `fq_lint/` | 15 | 19248 |
| `multiqc/` | 356 | 27862787 |
| `pipeline_info/` | 6 | 5239369 |
| `salmon/` | 77 | 223466 |
| `star_salmon/` | 524 | 34709523 |
| `trimgalore/` | 8 | 21933 |

## Scope

Can test: the files M2's read stage and M3 read from an aligned outdir (nf-core's own MultiQC tables, both Salmon matrices, params and traces), two launches of one outdir, a sample with two runs merged, single- and paired-end samples, auto and reverse strandedness.

Cannot test: GENCODE and `--gencode`; human reads, genes and biotypes; an unstranded or a failed strand check; the cluster's Nextflow with Apptainer.
