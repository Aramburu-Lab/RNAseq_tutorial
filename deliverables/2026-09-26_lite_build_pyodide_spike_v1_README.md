# 2026-09-26_lite_build_pyodide_spike_v1

Written by `ctx`. Do not edit; rerun the command instead.

| Field | Value |
| --- | --- |
| Created (UTC) | 2026-09-26T18:52:15Z |
| Command | `scripts/python/spike_pyodide.py --node-dir /private/tmp/claude-505/p15_spike/node --workdir /private/tmp/claude-505/p15_spike/run2` |
| ctx | 0.1.0, git 4a54ff07227e798cbe93fa4cd01738823b1be1b5 |
| Rows | 378 |
| Columns | 4 |
| Estimated in-memory size (bytes) | 20261 |

## Files

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `2026-09-26_lite_build_pyodide_spike_v1.csv` | 21809 | `e838f9b49d46ef800496215fbb340209c1081cf20a33ed2ead230c4fb7473322` |

## Data releases

| Resource | Release |
| --- | --- |
| Pyodide | 314.0.6, cdn.jsdelivr.net/pyodide/v314.0.6/full/ |
| PyPI | the wheels listed in the table, downloaded 2026-09-26 |

## Software versions

| Package | Version |
| --- | --- |
| ctx | 0.1.0 |
| numpy | 2.5.3 |
| polars | 1.44.2 |
| polars-runtime-32 | 1.44.2 |
| python | 3.12.14 |
| setuptools | 84.0.0 |
| typing_extensions | 4.16.0 |

## Data dictionary

| Column | Unit | Meaning | Type (csv) |
| --- | --- | --- | --- |
| `stage` | none | Source of the value: meta (inputs, tools, rule hashes), cpython (the reference in the pixi environment), stage_a (Pyodide 314.0.6 in Node), stage_b (Pyodide in a browser), decision (decide() applied to every row above it). | `String` |
| `component` | none | What was measured: a package, a grid size (fit_<genes>x<samples> for pydeseq2 as the CLI calls it, fit_shim_<size> for pydeseq2 under the in-process joblib backend, which decide() does not read, scipy_tests_<size> for the fallback), an input (input_<size>) or a downloaded file (<group>:<file name>). | `String` |
| `field` | none | The quantity. Its suffix is the unit of value: _s seconds, _bytes bytes, _n count, _log2 log2 fold change; any other field is text or a dimensionless number. | `String` |
| `value` | given by the field suffix | The recorded value, as text. | `String` |
