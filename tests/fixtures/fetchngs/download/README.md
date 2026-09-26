# nf-core/fetchngs fixture: download

A real fetchngs run copied as written, never edited. Built by `scripts/python/make_fetchngs_fixture.py`.

| Field | Value |
| --- | --- |
| Command | `NXF_SYNTAX_PARSER=v1 nextflow run nf-core/fetchngs -r 1.12.0 -profile docker,arm --input out/ids.csv --outdir out --nf_core_pipeline rnaseq -work-dir <scratch>/work --max_cpus 4 --max_memory 8.GB` (run in a scratch directory) |
| nf-core/fetchngs | v1.12.0-g8ec2d93 |
| Nextflow | 26.04.6 |
| Profile | docker,arm |
| Host | Darwin 24.6.0 arm64 |
| Date (UTC) | 2026-09-26 |
| Input list | `ids.csv`, from the four IDs in DOWNLOAD_IDS, sha256 `7957fa42ec64e35aaa4bed43ead6a45f8ba88c8faaa8a2a1e742796f28526ec5` |

## Left out

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `pipeline_info/execution_report_2026-09-26_22-56-31.html` | 1663953 | `ca675fe3aa11fd8e17aa186715f3fda59b436772979b664abb6437c5f5e7b35e` |
| `pipeline_info/execution_timeline_2026-09-26_22-56-31.html` | 256684 | `cebd835bb809dab7cc65d6eb01345520050eee0d4fe4b36d629099223260e691` |
| `pipeline_info/pipeline_dag_2026-09-26_22-56-31.html` | 3249 | `74e6a94ba3c9ceff674f4a35b2d9df16791b77397618989f5a603e89b7b20221` |

## Scope

Can test: verified and no_md5 files, paired and single-end, FTP and the sra-tools fallback, sheet paths that point nowhere after the move.

Cannot test: runs with more than two ENA files, Aspera, a -resume into the same outdir, any cohort's sample titles, and any fetchngs release but the pinned one.
