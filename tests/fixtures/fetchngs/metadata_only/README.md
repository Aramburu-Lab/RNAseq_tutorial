# nf-core/fetchngs fixture: metadata_only

A real fetchngs run copied as written, never edited. Built by `scripts/python/make_fetchngs_fixture.py`.

| Field | Value |
| --- | --- |
| Command | `NXF_SYNTAX_PARSER=v1 nextflow run nf-core/fetchngs -r 1.12.0 -profile docker,arm --input out/ids.csv --outdir out --nf_core_pipeline rnaseq -work-dir <scratch>/work --max_cpus 4 --max_memory 8.GB --skip_fastq_download` (run in a scratch directory) |
| nf-core/fetchngs | v1.12.0-g8ec2d93 |
| Nextflow | 26.04.6 |
| Profile | docker,arm |
| Host | Darwin 24.6.0 arm64 |
| Date (UTC) | 2026-09-26 |
| Input list | `ids.csv`, from https://raw.githubusercontent.com/nf-core/test-datasets/2732b911c57e607fa7aea5ba0c3d91b25bafb662/testdata/v1.12.0/sra_ids_test.csv, sha256 `bf6f1b2fc582c63f55914f27856f301e6d08c9f039cad5d326b9b810ce2c5fbf` |

## Left out

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `pipeline_info/execution_report_2026-09-26_22-57-58.html` | 1665434 | `d2cb1378a2f0504f7662b37432768d18aafef4a3bc9c7506b77b3ef68dde9454` |
| `pipeline_info/execution_timeline_2026-09-26_22-57-58.html` | 260368 | `9b0805edbf4d99be903566ccda96af8f268f799e9402dc671677afb08b135259` |
| `pipeline_info/pipeline_dag_2026-09-26_22-57-58.html` | 2440 | `7a720f608f13565a10bc6706217a80842f99a7d73833cf5b63dc07d0608340dd` |

## Scope

Can test: a metadata-only outdir (not_downloaded), a multi-run experiment (SRX9504942, 4 runs), GEO and DDBJ input IDs.

Cannot test: runs with more than two ENA files, Aspera, a -resume into the same outdir, any cohort's sample titles, and any fetchngs release but the pinned one.
