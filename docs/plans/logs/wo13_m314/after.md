# After: this change, loadgroup

Command: `pytest -q -n auto --dist loadgroup --durations=80`

| Measure | Value |
| --- | --- |
| Result | 2122 passed, 10 skipped, 2 warnings in 1192.30s (0:19:52) |
| Wall time (`/usr/bin/time -l`) | 1,193 s |
| CPU time, user + sys | 2,640 + 732 = 3,372 s |
| Largest resident set | 1,119,256,576 B |
| Load average, 1/5/15 min, at start | 21.29 47.02 58.45 |
| Load average, 1 min, during (every 30 s) | 21.3 to 104.8 |
| Load average, 1/5/15 min, at end | 68.61 54.84 53.75 |
| Worker busy time, least to most | 942 to 1,171 s over 12 workers |

## Session runs built through `conftest.once`

Builds counts the run's `session_runs.log`, which `built_twice` reads; wall is the build on the worker that took the lock; the others waited on it.

| Name | Builds | Built on | Build wall (s) | Workers that asked | Longest wait (s) |
| --- | --- | --- | --- | --- | --- |
| handoff | 1 | gw5 | 208.6 | 1 | 208.6 |
| pipeline | 1 | gw11 | 207.2 | 12 | 207.3 |
| qc_fastqc | 1 | gw3 | 205.8 | 1 | 206.2 |
| qc_reads | 1 | gw1 | 196.7 | 4 | 197.3 |
| qc_views | 1 | gw4 | 146.2 | 5 | 146.3 |
| process | 1 | gw1 | 139.1 | 10 | 139.6 |
| process_again | 1 | gw9 | 130.9 | 3 | 133.9 |
| qc_planted | 1 | gw5 | 83.4 | 3 | 86.6 |
| qc_reports | 1 | gw1 | 83.2 | 2 | 83.3 |
| proteo_made | 1 | gw10 | 82.9 | 3 | 83.0 |
| proteo_study_maxquant | 1 | gw8 | 77.4 | 2 | 77.7 |
| ncrna_tcga | 1 | gw11 | 72.6 | 5 | 72.8 |
| proteo_across | 1 | gw2 | 53.3 | 7 | 53.7 |
| quant | 1 | gw8 | 51.0 | 6 | 51.5 |
| m10_options | 1 | gw8 | 44.8 | 1 | 45.8 |
| ncrna_options | 1 | gw2 | 44.2 | 1 | 44.4 |
| proteo_study_diann | 1 | gw11 | 40.0 | 1 | 40.3 |
| qc_cli | 1 | gw10 | 38.4 | 6 | 40.5 |
| clinical_cli | 1 | gw6 | 37.0 | 8 | 40.4 |
| m9_figures | 1 | gw9 | 35.4 | 3 | 36.2 |
| ncrna_cli | 1 | gw4 | 25.0 | 5 | 26.4 |
| proteo_cli | 1 | gw4 | 20.1 | 5 | 24.2 |
| gather_root | 1 | gw10 | 13.2 | 4 | 13.5 |
| genecard_cli | 1 | gw7 | 12.7 | 3 | 12.8 |
| guide_cli | 1 | gw0 | 7.3 | 5 | 7.7 |
| root | 1 | gw10 | 3.3 | 5 | 3.6 |
| wheel | 1 | gw11 | 2.4 | 12 | 2.5 |

## The 40 slowest setups and calls

A setup's time includes building or waiting for the session runs it is the first on its worker to ask for. CPU is the whole test's, self and reaped children, where measured.

| Wall (s) | Phase | Worker | Test | CPU of the test (s) | Shared fixtures it reaches |
| --- | --- | --- | --- | --- | --- |
| 397.4 | call | gw3 | `tests/test_pipeline.py::test_run_matches_the_module_commands_through_de_and_function` | 110.9 | none |
| 308.0 | setup | gw11 | `tests/test_de_views_browser.py::test_every_grid_setting_draws_what_classify_counts[wide_light]` | 150.3 | pipeline, wheel_site; module: test_de_views_browser.m5 |
| 306.4 | setup | gw5 | `tests/test_handoff_browser.py::test_a_clicked_link_shows_the_genes_on_m6` | 166.4 | none; module: test_handoff_browser.ran |
| 284.4 | setup | gw3 | `tests/test_function_views_browser.py::test_every_set_draws_lites_running_sum_and_its_own_words[wide_light]` | 15.1 | pipeline, wheel_site; module: test_function_views_browser.m6 |
| 250.5 | call | gw1 | `tests/test_figures_m2.py::test_report_figures_redraw_byte_for_byte` | 124.0 | with_reports |
| 207.5 | setup | gw1 | `tests/test_process_views_browser.py::test_the_scatter_draws_each_type_and_sample[wide_light]` | 82.5 | process_root; module: test_process_views_browser.m3 |
| 206.2 | setup | gw3 | `tests/test_figures_m2.py::test_fastqc_plots_redraw_from_their_tables` | 80.1 | with_fastqc |
| 201.1 | call | gw7 | `tests/test_gather_counts.py::test_gather_qc_de_run_on_a_matrix` | 46.3 | none |
| 197.6 | setup | gw1 | `tests/test_figures_m2.py::test_reads_figures_redraw_and_shared_stems_have_no_spec` | 65.5 | with_reads |
| 181.9 | call | gw5 | `tests/test_figures_m3.py::test_every_figure_redraws_byte_for_byte` | 154.6 | process_root |
| 173.6 | call | gw8 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[study_root]` | 33.3 | study_root (by request) |
| 169.6 | call | gw2 | `tests/test_figures_m7.py::test_a_run_with_other_options_redraws_with_them` | 56.3 | across_root, proteo_across |
| 169.2 | call | gw10 | `tests/test_figures_m10.py::test_every_figure_redraws_byte_for_byte` | 45.5 | clinical_run |
| 162.6 | call | gw10 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[made]` | 42.5 | made (by request) |
| 157.2 | call | gw4 | `tests/test_figures_m9.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte` | 47.8 | built |
| 152.1 | call | gw4 | `tests/test_page_cap.py::test_m2_drops_its_flag_view_under_a_cap_one_byte_below_its_page` | 72.3 | page |
| 146.3 | setup | gw4 | `tests/test_page_cap.py::test_m2_drops_its_flag_view_under_a_cap_one_byte_below_its_page` | 72.3 | page |
| 144.8 | call | gw0 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[across_root]` | 29.1 | across_root (by request) |
| 142.4 | call | gw10 | `tests/test_pipeline.py::test_run_matches_the_module_commands` | 15.6 | none |
| 133.9 | call | gw9 | `tests/test_cli.py::test_parameter_cards_match_cli[process]` | 59.0 | process_run (by request) |
| 132.7 | setup | gw0 | `tests/test_proteo_views_browser.py::test_each_fold_change_regroups_the_genes_as_classify_does[wide_light]` | 31.0 | none; module: test_proteo_views_browser.m7 |
| 131.1 | call | gw11 | `tests/test_figures_m8.py::test_gdc_figures_redraw_byte_for_byte` | 62.4 | ncrna_tcga_run, tcga |
| 125.1 | setup | gw7 | `tests/test_process_count.py::test_the_strand_is_called_from_the_forward_share` | 26.0 | none; module: test_process_count.run |
| 114.9 | call | gw11 | `tests/test_qc_run.py::test_fetchngs_reads_alone_within_the_track_a_budget` | 31.6 | none |
| 112.8 | call | gw8 | `tests/test_page_width_browser.py::test_no_fixture_page_is_wider_than_a_phone` | 25.9 | none |
| 105.7 | setup | gw8 | `tests/test_view_browser.py::test_a_view_swaps_the_figure_and_its_words_and_comes_back[m5]` | 5.0 | pipeline, wheel_site; module: test_view_browser.ran |
| 103.8 | call | gw7 | `tests/test_figures_m9.py::test_a_page_short_of_cards_redraws_the_whole_list` | 15.7 | none |
| 102.1 | call | gw6 | `tests/test_figures.py::test_gather_figures_redraw_from_the_sheet` | 16.5 | none |
| 102.0 | setup | gw2 | `tests/test_ncrna_views_browser.py::test_the_volcano_calls_each_biotype_at_the_readers_cutoffs[wide_light]` | 36.2 | none; module: test_ncrna_views_browser.m8 |
| 98.7 | call | gw10 | `tests/test_figures.py::test_every_spec_redraws_its_data_byte_for_byte` | 80.7 | pipeline, wheel_site |
| 96.8 | setup | gw10 | `tests/test_pick_browser.py::test_a_real_drag_picks_the_genes_inside_it[False-1280-box]` | 21.8 | none; module: test_pick_browser.ran |
| 96.8 | call | gw0 | `tests/test_page_cap.py::test_m9_decides_its_cards_on_the_page_with_its_panels` | 41.4 | none |
| 93.5 | call | gw9 | `tests/test_own_fastq.py::test_run_takes_the_route` | 9.4 | none |
| 92.4 | call | gw8 | `tests/test_gates.py::test_zero_collected_fails` | 6.4 | none |
| 88.5 | call | gw3 | `tests/test_figures_m8.py::test_the_mirna_sheet_passes_the_gdc_door` | 8.8 | ncrna_tcga_run, tcga |
| 86.6 | setup | gw5 | `tests/test_qc_run.py::test_the_planted_problems_are_flagged` | 25.8 | ran |
| 85.6 | setup | gw7 | `tests/test_qc_run.py::test_size_against_genes_marks_the_samples_flagged_on_its_axes` | 0.0 | ran |
| 84.6 | call | gw6 | `tests/test_figures_m7.py::test_a_cross_study_n_counts_each_level_of_the_contrast` | 5.8 | cross_root, proteo_cross, proteo_run |
| 83.9 | setup | gw9 | `tests/test_qc_run.py::test_a_second_run_writes_every_table_byte_for_byte` | 28.1 | ran |
| 83.3 | setup | gw1 | `tests/test_figures_m2.py::test_report_figures_redraw_byte_for_byte` | 124.0 | with_reports |

## Module fixtures built more than once

A module fixture is built once per worker and per file that requests it.

| Fixture | Builds | Files that request it, with the workers |
| --- | --- | --- |

Module fixtures built once: 70.

## Files by summed setup and call time

| File | Summed wall (s) | CPU (s) |
| --- | --- | --- |
| `tests/test_figures_m2.py` | 1028.9 | 344.5 |
| `tests/test_figures_m7.py` | 820.4 | 185.0 |
| `tests/test_pipeline.py` | 592.2 | 132.3 |
| `tests/test_qc_run.py` | 576.0 | 133.4 |
| `tests/test_figures_m8.py` | 514.9 | 114.1 |
| `tests/test_figures_m9.py` | 433.1 | 119.5 |
| `tests/test_page_cap.py` | 395.4 | 113.8 |
| `tests/test_de_views_browser.py` | 330.6 | 157.7 |
| `tests/test_gather_counts.py` | 317.4 | 58.2 |
| `tests/test_cli.py` | 309.1 | 113.9 |
| `tests/test_handoff_browser.py` | 306.8 | 166.5 |
| `tests/test_function_views_browser.py` | 302.5 | 19.9 |
| `tests/test_figures_m10.py` | 267.8 | 82.9 |
| `tests/test_kit.py` | 258.2 | 17.9 |
| `tests/test_figures.py` | 255.1 | 109.9 |
| `tests/test_process_views_browser.py` | 207.8 | 82.8 |
| `tests/test_gates.py` | 195.3 | 18.1 |
| `tests/test_process_align.py` | 189.4 | 36.3 |
| `tests/test_figures_m3.py` | 182.4 | 154.6 |
| `tests/test_process_count.py` | 178.4 | 35.8 |
