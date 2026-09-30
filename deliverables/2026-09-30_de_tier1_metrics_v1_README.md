# 2026-09-30_de_tier1_metrics_v1

Written by `ctx`. Do not edit; rerun the command instead.

| Field | Value |
| --- | --- |
| Created (UTC) | 2026-09-30T09:45:40Z |
| Command | `scripts/python/record_sanity.py --root /Users/rytis/ctx_tier1 --dataset tier1 --dataset-file conf/datasets/tier1.toml` |
| ctx | 0.1.0, git d3d4f3973e0e0cdbf4ca3d1cf6a9c1e89d4b09b2 |
| Directory | /private/tmp/claude-505/-Users-rytis-Scripts-Pipeline-RNAseq-tutorial/bd98644a-984e-40e7-be65-bc8c990c4492/scratchpad/record_wt/deliverables |
| Rows | 92 |
| Columns | 4 |
| Estimated in-memory size (bytes) | 7139 |

## Files

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `2026-09-30_de_tier1_metrics_v1.csv` | 7561 | `641121260c9a0663a77d7228743800a5a248b2173f8749b99045c538d9a90453` |

## Data releases

| Resource | Release |
| --- | --- |
| hub verdicts | 2026-09-30_report_tier1_sanity_v1.csv |
| sha256:2026-09-30_report_tier1_sanity_v1.csv | ebb82fd7aa98aca35581cf77304b7a90ea026a24945bd024d423916c4a59812e |
| sha256:2026-09-30_de_tier1_metrics_v1.parquet | a057153c772392679b4daef38d5a183914aa853bc1fd2c1948a95cb602bb83ef |
| sha256:2026-09-30_function_tier1_metrics_v1.parquet | 80217d022ca7776ce998e4eac3b0bbab8015ffaafc35d9abc7997d6a035b6f47 |
| sha256:tier1.toml | d3c37ccc7080cdf82a7975d7f08e2e8d3b7ae3b213fe2c6971c3b181a902941a |
| rule commit | 3620021a0a51b8efbd257e12d9be75790f2ed526 |

## Software versions

| Package | Version |
| --- | --- |
| Jinja2 | 3.1.6 |
| MarkupSafe | 3.0.3 |
| charset-normalizer | 3.5.1 |
| ctx | 0.1.0 |
| numpy | 2.5.3 |
| polars | 1.44.2 |
| polars-runtime-32 | 1.44.2 |
| python | 3.12.14 |
| scipy | 1.18.1 |
| setuptools | 84.0.0 |
| typing_extensions | 4.16.0 |

## Data dictionary

| Column | Unit | Meaning | Type (csv) |
| --- | --- | --- | --- |
| `module` | none | The module whose metrics table holds the metric: de (M5), function (M6). | `String` |
| `metric` | none | Name the sanity rules read, e.g. permanova.<pca>.<factor>.r2. | `String` |
| `value` | varies | The measured value; the unit follows from the metric's name. | `Float64` |
| `detail` | none | What the value was computed on, or why it is missing. | `String` |
