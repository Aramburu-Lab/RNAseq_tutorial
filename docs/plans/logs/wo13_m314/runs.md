# M3.14: the suite's cost, measured before and after (2026-10-21)

All runs on this Mac (14 cores, 24 GB), `-n auto` = 12 workers, from a clean clone of the
branch in the session's scratch folder (`suite.sh`, the probe `m314_probe.py`, the tables
`m314_table.py`). The machine was shared: other agents' suites and Microsoft Defender's
real-time scan (`wdavdaemon_unprivileged` at 304 % CPU at 16:16) ran throughout, so every
wall time here is inflated, and CPU time is the steadier measure. The full tables are
[before.md](before.md), [after.md](after.md) and [final.md](final.md), the run of the
committed code; every session run in it was built once and every module fixture once.

## The four full runs

| Run | Code | `--dist` | Result | Wall (s) | CPU user + sys (s) | Load, 1 min, during |
| --- | --- | --- | --- | --- | --- | --- |
| before | 9ac06d1 | loadscope | 2,118 passed, 10 skipped | 3,099 | 4,673 | 19.6 to 55.8 |
| v1 | `versions()` reads the installed map once | loadscope | 1 failed, 2,118 passed | 1,796 | 3,553 | 15.5 to 46.8 |
| v2 | v1, interrogate fixed, shared fixtures through `once`, loadgroup | loadgroup | 2,119 passed | 2,371 | 5,196 | 20.0 to 160.5 |
| after | v2, the grouping fixed, MultiQC and redraws 4 at a time | loadgroup | 2,122 passed, 10 skipped | 1,193 | 3,372 | 21.3 to 104.8 |
| final | after, two docstrings and these logs: the commit | loadgroup | 2,122 passed, 10 skipped | 923 | 3,182 | 22.0 to 54.2 |

v1's failure was interrogate: the two new runtime helpers had no docstring. v2 grouped
nothing: the hook read `--dist`, which an xdist worker sees as `no`, so every module was
spread and each browser module launched Chromium on every worker that ran one of its tests
(`test_export_browser.py` 51 s CPU in v1, 644 s in v2). `test_session_runs.py` now fails on
that. The after run's workers were busy 942 to 1,171 s each, the final run's 765 to 857 s,
against the before run's 1,840 to 3,031 s.

## The ten slowest files, each alone without xdist

The ten files of most summed time in the before run, each in its own `pytest -p no:xdist`
process (`serial.sh`), the ten of each tree and both trees at once: 20 processes, load 18.9 to
39.6. Each builds its own session runs.

| File | Before wall (s) | Before CPU (s) | After wall (s) | After CPU (s) | Result after |
| --- | --- | --- | --- | --- | --- |
| `test_figures_m2.py` | 1,535 | 414 | 875 | 357 | 6 passed |
| `test_cli.py` | 1,494 | 424 | 1,050 | 348 | 41 passed |
| `test_qc_run.py` | 1,136 | 300 | 749 | 227 | 16 passed |
| `test_figures_m7.py` | 1,198 | 220 | 684 | 189 | 9 passed |
| `test_outputs_contract.py` | 1,274 | 340 | 831 | 286 | 27 passed |
| `test_figures.py` | 913 | 260 | 701 | 251 | 11 passed |
| `test_process_run.py` | 844 | 151 | 458 | 142 | 20 passed, 2 skipped |
| `test_pipeline.py` | 839 | 185 | 659 | 142 | 15 passed |
| `test_figures_m8.py` | 836 | 147 | 491 | 124 | 8 passed |
| `test_figures_m9.py` | 843 | 165 | 393 | 117 | 6 passed |
| Sum | 10,911 | 2,606 | 6,891 | 2,183 | |

## What cost what, and what changed

| Cost found | Measured | Change |
| --- | --- | --- |
| `runtime.versions()` read every installed distribution's files for each README | 1.9 to 2.8 s a call at a load of 20; 2,374 READMEs left in the before run's basetemp; 45 in one fixture `ctx process` | the module-to-distribution map and each version read once per process; the modules listed stay live |
| `--dist loadscope` held a module's tests on one worker | `test_figures_m2.py` 1,813 s on one worker; busiest worker 3,031 s of 3,099 | `--dist loadgroup`: only the tests that share a module fixture stay together (`conftest.module_group`) |
| Module fixtures built once per requesting file | `test_qc_views` page 3 times; `with_reads`, `with_reports`, `clinical_root` twice | session fixtures through `once` |
| The same run under two fixture names | M7's cross-dataset root 3 times, across-samples, tumour-normal and MaxQuant study twice; `ctx ncrna` on each made-up root twice; `ctx clinical` twice | one `once` run each, in conftest |
| 11 MultiQC runs one after another in `ctx process` | 247 of the 324 s fixture run in v1, 238 s of the rerun | runs of one stage `threads` at a time (`multiqc_io.at_once`): the build 139 s and the rerun 131 s in the after run |
| Byte-for-byte redraws one `ctx figure` at a time | 13.7 to 15.6 s CPU for each of 4 M3 redraws timed alone | `REDRAWS = 4` at a time (`test_figures.redrawn_each`) |

## Not done, and the next largest costs

- Under 5 minutes is not shown: the final run took 923 s at a load of 22 to 54. Its CPU,
  3,182 s over 12 workers, is 265 s if nothing else ran; the measurement on a quiet machine
  is still to be made.
- `test_pipeline.py::test_run_matches_the_module_commands_through_de_and_function`, 397 s:
  it runs de and function twice, through `ctx run` and module by module, which is what it
  tests.
- The DE fits: `pipeline` 207 s and `handoff` 209 s, each de and function on the same
  synthetic design, one with made-up Ensembl IDs and Hallmark set names. Sharing one fit
  would need the pipeline's tests to read those IDs.
- A `ctx figure` redraw costs 13.7 to 15.6 s CPU; a stack sampler on one M3 redraw put 67 %
  of its wall time in imports and 19 % in plotly's validated `update_layout`. The after run's
  basetemp holds 158 redraw folders, counted by their `params.toml`.
- CI still runs serially on 2 cores (D-060); the `once` sharing, the MultiQC and the redraw
  concurrency apply there too, but no CI time was measured.
