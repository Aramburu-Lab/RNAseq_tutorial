# Before: 9ac06d1, loadscope

Command: `pytest -q -n auto --dist loadscope --durations=80`

| Measure | Value |
| --- | --- |
| Result | 2118 passed, 10 skipped, 3 warnings in 3089.93s (0:51:29) |
| Wall time (`/usr/bin/time -l`) | 3,099 s |
| CPU time, user + sys | 3,431 + 1,242 = 4,673 s |
| Largest resident set | 1,133,035,520 B |
| Load average, 1/5/15 min, at start | 42.45 24.14 17.35 |
| Load average, 1 min, during (every 30 s) | 19.6 to 55.8 |
| Load average, 1/5/15 min, at end | 30.99 26.77 26.48 |
| Worker busy time, least to most | 1,840 to 3,031 s over 12 workers |

## Session runs built through `conftest.once`

Builds counts the run's `session_runs.log`, which `built_twice` reads; wall is the build on the worker that took the lock; the others waited on it.

| Name | Builds | Built on | Build wall (s) | Workers that asked | Longest wait (s) |
| --- | --- | --- | --- | --- | --- |
| process | 1 | gw1 | 569.5 | 5 | 570.0 |
| handoff | 1 | gw8 | 514.4 | 1 | 515.6 |
| process_again | 1 | gw7 | 347.7 | 3 | 350.1 |
| pipeline | 1 | gw6 | 329.2 | 6 | 329.3 |
| qc_cli | 1 | gw7 | 191.3 | 4 | 195.9 |
| clinical_cli | 1 | gw7 | 142.4 | 3 | 144.4 |
| quant | 1 | gw11 | 137.8 | 5 | 138.2 |
| proteo_cli | 1 | gw11 | 66.3 | 3 | 69.1 |
| ncrna_cli | 1 | gw8 | 58.1 | 4 | 61.0 |
| genecard_cli | 1 | gw11 | 32.2 | 3 | 33.1 |
| gather_root | 1 | gw11 | 21.1 | 3 | 21.5 |
| root | 1 | gw11 | 18.1 | 4 | 18.4 |
| wheel | 1 | gw6 | 6.9 | 6 | 7.1 |
| guide_cli | 1 | gw11 | 5.6 | 3 | 6.7 |

## The 40 slowest setups and calls

A setup's time includes building or waiting for the session runs it is the first on its worker to ask for. CPU is the whole test's, self and reaped children, where measured.

| Wall (s) | Phase | Worker | Test | CPU of the test (s) | Shared fixtures it reaches |
| --- | --- | --- | --- | --- | --- |
| 794.1 | call | gw7 | `tests/test_cli.py::test_parameter_cards_match_cli[process]` |  | process_run (by request) |
| 595.2 | setup | gw8 | `tests/test_handoff_browser.py::test_a_clicked_link_shows_the_genes_on_m6` |  | none; module: test_handoff_browser.ran |
| 570.0 | setup | gw1 | `tests/test_process_run.py::test_every_output_is_named_described_and_offline` |  | process_root |
| 561.3 | call | gw0 | `tests/test_pipeline.py::test_run_matches_the_module_commands_through_de_and_function` |  | none |
| 474.9 | call | gw1 | `tests/test_figures_m2.py::test_report_figures_redraw_byte_for_byte` |  | none; module: test_qc_run.with_reports |
| 460.5 | call | gw6 | `tests/test_figures.py::test_every_spec_redraws_its_data_byte_for_byte` |  | pipeline, wheel_site |
| 445.7 | call | gw11 | `tests/test_outputs_contract.py::test_every_output_named_and_described[process_root]` |  | process_root (by request) |
| 410.4 | call | gw3 | `tests/test_figures_m3.py::test_every_figure_redraws_byte_for_byte` |  | process_root |
| 385.0 | call | gw10 | `tests/test_page_cap.py::test_m9_decides_its_cards_on_the_page_with_its_panels` |  | none |
| 370.5 | call | gw2 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[made]` |  | made (by request) |
| 366.3 | call | gw9 | `tests/test_gather_counts.py::test_gather_qc_de_run_on_a_matrix` |  | none |
| 348.6 | setup | gw1 | `tests/test_process_run.py::test_a_moved_rerun_reuses_every_tool_folder` |  | process_rerun, process_root |
| 336.4 | setup | gw6 | `tests/test_figures.py::test_every_spec_redraws_its_data_byte_for_byte` |  | pipeline, wheel_site |
| 335.5 | setup | gw1 | `tests/test_figures_m2.py::test_fastqc_plots_redraw_from_their_tables` |  | none; module: test_figures_m2.with_fastqc |
| 322.6 | setup | gw3 | `tests/test_de_views_browser.py::test_every_grid_setting_draws_what_classify_counts[wide_light]` |  | pipeline, wheel_site; module: test_de_views_browser.m5 |
| 311.4 | call | gw4 | `tests/test_figures_m9.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte` |  | none; module: test_figures_m9.built |
| 283.5 | setup | gw10 | `tests/test_qc_run.py::test_reads_are_judged_run_by_run` |  | none; module: test_qc_run.with_reads |
| 271.7 | call | gw4 | `tests/test_figures_m9.py::test_a_second_build_in_one_process_carries_no_panel_of_the_first` |  | none |
| 269.1 | call | gw11 | `tests/test_figures_m10.py::test_every_figure_redraws_byte_for_byte` |  | clinical_run |
| 266.9 | setup | gw10 | `tests/test_page_cap.py::test_m2_drops_its_flag_view_under_a_cap_one_byte_below_its_page` |  | none; module: test_qc_views.page |
| 264.5 | call | gw2 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[across_root]` |  | across_root (by request) |
| 262.5 | call | gw8 | `tests/test_figures_m8.py::test_gdc_figures_redraw_byte_for_byte` |  | none; module: test_figures_m8.tcga |
| 254.6 | call | gw1 | `tests/test_figures_m2.py::test_every_count_figure_redraws_byte_for_byte` |  | qc_run |
| 253.9 | setup | gw3 | `tests/test_proteo_views_browser.py::test_each_fold_change_regroups_the_genes_as_classify_does[wide_light]` |  | none; module: test_proteo_views_browser.m7 |
| 243.8 | call | gw11 | `tests/test_outputs_contract.py::test_every_output_named_and_described[de_run]` |  | de_run (by request) |
| 232.2 | setup | gw5 | `tests/test_qc_views.py::test_the_sex_view_compares_each_sample_with_its_own_record` |  | none; module: test_qc_views.page |
| 230.7 | call | gw0 | `tests/test_pipeline.py::test_run_matches_the_module_commands` |  | none |
| 225.2 | setup | gw9 | `tests/test_figure_panel.py::test_every_spec_has_its_panel_on_the_module_page` |  | pipeline, wheel_site |
| 221.1 | call | gw2 | `tests/test_figures_m7.py::test_a_run_with_other_options_redraws_with_them` |  | none; module: test_figures_m7.across_root |
| 220.2 | setup | gw4 | `tests/test_function_views_browser.py::test_every_set_draws_lites_running_sum_and_its_own_words[wide_light]` |  | pipeline, wheel_site; module: test_function_views_browser.m6 |
| 214.9 | setup | gw1 | `tests/test_figures_m2.py::test_report_figures_redraw_byte_for_byte` |  | none; module: test_qc_run.with_reports |
| 214.0 | call | gw3 | `tests/test_genecard.py::test_cli_offline_reads_only_the_cache` |  | none |
| 212.3 | call | gw2 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[cross_root]` |  | cross_root (by request) |
| 210.1 | call | gw2 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[study_root]` |  | study_root (by request) |
| 209.4 | setup | gw2 | `tests/test_process_count.py::test_the_strand_is_called_from_the_forward_share` |  | none; module: test_process_count.run |
| 204.0 | setup | gw10 | `tests/test_qc_run.py::test_the_planted_problems_are_flagged` |  | none; module: test_qc_run.ran |
| 201.2 | call | gw7 | `tests/test_genecard_absent.py::test_offline_rerun_prints_the_online_pages_no_entry` |  | none |
| 199.5 | call | gw0 | `tests/test_genecard_pack.py::test_cards_after_the_first_travel_packed_and_read_back_whole` |  | none |
| 196.0 | call | gw7 | `tests/test_cli.py::test_parameter_cards_match_cli[qc]` |  | qc_run (by request) |
| 196.0 | call | gw11 | `tests/test_outputs_contract.py::test_every_output_named_and_described[qc_run]` |  | qc_run (by request) |

## Module fixtures built more than once

A module fixture is built once per worker and per file that requests it.

| Fixture | Builds | Files that request it, with the workers |
| --- | --- | --- |
| test_qc_views.page | 3 | test_page_cap.py (gw10); test_qc_judge.py (gw9); test_qc_views.py (gw5) |
| test_clinical_page.clinical_root | 2 | test_clinical_live.py (gw5); test_clinical_page.py (gw5) |
| test_qc_run.with_reads | 2 | test_figures_m2.py (gw1); test_qc_run.py (gw10) |
| test_qc_run.with_reports | 2 | test_figures_m2.py (gw1); test_qc_run.py (gw10) |

Module fixtures built once: 83.

## Files by summed setup and call time

| File | Summed wall (s) | CPU (s) |
| --- | --- | --- |
| `tests/test_figures_m2.py` | 1813.0 | |
| `tests/test_cli.py` | 1365.4 | |
| `tests/test_qc_run.py` | 1318.0 | |
| `tests/test_figures_m7.py` | 1310.2 | |
| `tests/test_outputs_contract.py` | 1268.1 | |
| `tests/test_figures.py` | 1041.1 | |
| `tests/test_process_run.py` | 940.8 | |
| `tests/test_pipeline.py` | 883.9 | |
| `tests/test_figures_m8.py` | 873.1 | |
| `tests/test_figures_m9.py` | 850.8 | |
| `tests/test_page_cap.py` | 836.6 | |
| `tests/test_handoff_browser.py` | 595.6 | |
| `tests/test_figures_m10.py` | 547.3 | |
| `tests/test_gather_counts.py` | 508.7 | |
| `tests/test_genecard.py` | 456.3 | |
| `tests/test_figures_m3.py` | 410.5 | |
| `tests/test_genecard_pack.py` | 376.2 | |
| `tests/test_de_views_browser.py` | 349.7 | |
| `tests/test_process_count.py` | 314.4 | |
| `tests/test_ncrna_tcga.py` | 296.0 | |
