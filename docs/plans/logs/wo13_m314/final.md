# Final: the committed code, loadgroup

Command: `pytest -q -n auto --dist loadgroup --durations=80`

| Measure | Value |
| --- | --- |
| Result | 2122 passed, 10 skipped, 2 warnings in 919.41s (0:15:19) |
| Wall time (`/usr/bin/time -l`) | 923 s |
| CPU time, user + sys | 2,519 + 662 = 3,182 s |
| Largest resident set | 1,119,813,632 B |
| Load average, 1/5/15 min, at start | 22.66 21.86 26.02 |
| Load average, 1 min, during (every 30 s) | 22.0 to 54.2 |
| Load average, 1/5/15 min, at end | 23.24 29.31 28.99 |
| Worker busy time, least to most | 765 to 857 s over 12 workers |

## Session runs built through `conftest.once`

Builds counts the run's `session_runs.log`, which `built_twice` reads; wall is the build on the worker that took the lock; the others waited on it.

| Name | Builds | Built on | Build wall (s) | Workers that asked | Longest wait (s) |
| --- | --- | --- | --- | --- | --- |
| handoff | 1 | gw5 | 294.1 | 1 | 294.2 |
| pipeline | 1 | gw11 | 287.7 | 10 | 287.8 |
| process_again | 1 | gw7 | 178.5 | 3 | 181.5 |
| process | 1 | gw2 | 161.1 | 7 | 162.2 |
| qc_fastqc | 1 | gw4 | 139.9 | 1 | 140.0 |
| qc_reports | 1 | gw5 | 104.7 | 2 | 105.0 |
| qc_reads | 1 | gw5 | 95.2 | 2 | 95.5 |
| qc_cli | 1 | gw1 | 78.6 | 5 | 81.2 |
| quant | 1 | gw8 | 71.3 | 7 | 71.4 |
| proteo_across | 1 | gw3 | 66.5 | 4 | 66.9 |
| ncrna_tcga | 1 | gw11 | 64.8 | 3 | 65.0 |
| qc_views | 1 | gw11 | 64.2 | 6 | 65.4 |
| m10_options | 1 | gw6 | 58.8 | 1 | 59.6 |
| clinical_cli | 1 | gw6 | 54.5 | 6 | 56.4 |
| proteo_made | 1 | gw4 | 46.6 | 4 | 46.9 |
| qc_planted | 1 | gw8 | 43.1 | 4 | 44.5 |
| m9_figures | 1 | gw7 | 40.2 | 3 | 40.9 |
| ncrna_options | 1 | gw9 | 35.8 | 1 | 35.9 |
| ncrna_cli | 1 | gw2 | 33.6 | 5 | 35.0 |
| proteo_cli | 1 | gw2 | 31.8 | 6 | 35.1 |
| proteo_study_maxquant | 1 | gw9 | 28.6 | 2 | 28.9 |
| genecard_cli | 1 | gw2 | 28.2 | 3 | 28.7 |
| proteo_study_diann | 1 | gw6 | 22.0 | 1 | 22.1 |
| wheel | 1 | gw11 | 11.1 | 10 | 11.2 |
| gather_root | 1 | gw1 | 9.8 | 3 | 10.0 |
| root | 1 | gw1 | 7.9 | 6 | 8.5 |
| guide_cli | 1 | gw0 | 6.5 | 4 | 7.7 |

## The 40 slowest setups and calls

A setup's time includes building or waiting for the session runs it is the first on its worker to ask for. CPU is the whole test's, self and reaped children, where measured.

| Wall (s) | Phase | Worker | Test | CPU of the test (s) | Shared fixtures it reaches |
| --- | --- | --- | --- | --- | --- |
| 377.2 | setup | gw5 | `tests/test_handoff_browser.py::test_a_clicked_link_shows_the_genes_on_m6` | 160.0 | none; module: test_handoff_browser.ran |
| 339.2 | setup | gw11 | `tests/test_de_views_browser.py::test_every_grid_setting_draws_what_classify_counts[wide_light]` | 144.0 | pipeline, wheel_site; module: test_de_views_browser.m5 |
| 308.3 | setup | gw3 | `tests/test_function_views_browser.py::test_every_set_draws_lites_running_sum_and_its_own_words[wide_light]` | 14.6 | pipeline, wheel_site; module: test_function_views_browser.m6 |
| 201.1 | setup | gw2 | `tests/test_process_views_browser.py::test_the_scatter_draws_each_type_and_sample[wide_light]` | 79.2 | process_root; module: test_process_views_browser.m3 |
| 181.5 | call | gw7 | `tests/test_cli.py::test_parameter_cards_match_cli[process]` | 60.7 | process_run (by request) |
| 174.2 | setup | gw7 | `tests/test_process_count.py::test_the_strand_is_called_from_the_forward_share` | 25.0 | none; module: test_process_count.run |
| 147.1 | call | gw5 | `tests/test_figures_m2.py::test_report_figures_redraw_byte_for_byte` | 113.4 | with_reports |
| 145.9 | setup | gw1 | `tests/test_proteo_views_browser.py::test_each_fold_change_regroups_the_genes_as_classify_does[wide_light]` | 29.3 | none; module: test_proteo_views_browser.m7 |
| 140.0 | setup | gw4 | `tests/test_figures_m2.py::test_fastqc_plots_redraw_from_their_tables` | 72.5 | with_fastqc |
| 131.8 | call | gw0 | `tests/test_pipeline.py::test_run_matches_the_module_commands_through_de_and_function` | 96.5 | none |
| 116.6 | setup | gw10 | `tests/test_ncrna_views_browser.py::test_the_volcano_calls_each_biotype_at_the_readers_cutoffs[wide_light]` | 33.4 | none; module: test_ncrna_views_browser.m8 |
| 115.9 | call | gw10 | `tests/test_figures.py::test_every_spec_redraws_its_data_byte_for_byte` | 81.5 | pipeline, wheel_site |
| 106.9 | call | gw8 | `tests/test_figures_m3.py::test_every_figure_redraws_byte_for_byte` | 136.3 | process_root |
| 106.1 | call | gw7 | `tests/test_gather_counts.py::test_gather_qc_de_run_on_a_matrix` | 42.6 | none |
| 105.0 | setup | gw5 | `tests/test_figures_m2.py::test_report_figures_redraw_byte_for_byte` | 113.4 | with_reports |
| 104.1 | call | gw9 | `tests/test_page_width_browser.py::test_no_fixture_page_is_wider_than_a_phone` | 25.1 | none |
| 101.2 | call | gw1 | `tests/test_figures.py::test_gather_figures_redraw_from_the_sheet` | 16.4 | none |
| 95.5 | setup | gw5 | `tests/test_figures_m2.py::test_reads_figures_redraw_and_shared_stems_have_no_spec` | 61.8 | with_reads |
| 93.5 | call | gw11 | `tests/test_page_cap.py::test_m2_drops_its_flag_view_under_a_cap_one_byte_below_its_page` | 67.1 | page |
| 90.0 | setup | gw8 | `tests/test_quant_views_browser.py::test_each_threshold_draws_its_own_count_of_the_same_draws[wide_light]` | 19.5 | quant_built; module: test_quant_views_browser.m4 |
| 89.3 | setup | gw10 | `tests/test_pick_browser.py::test_a_real_drag_picks_the_genes_inside_it[False-1280-box]` | 20.6 | none; module: test_pick_browser.ran |
| 81.2 | call | gw1 | `tests/test_cli.py::test_parameter_cards_match_cli[qc]` | 22.0 | qc_run (by request) |
| 81.1 | setup | gw4 | `tests/test_view_browser.py::test_a_view_swaps_the_figure_and_its_words_and_comes_back[m5]` | 4.7 | pipeline, wheel_site; module: test_view_browser.ran |
| 81.1 | call | gw3 | `tests/test_figures_m7.py::test_a_run_with_other_options_redraws_with_them` | 50.5 | across_root, proteo_across |
| 81.0 | call | gw1 | `tests/test_page_cap.py::test_m9_decides_its_cards_on_the_page_with_its_panels` | 39.7 | none |
| 80.3 | setup | gw6 | `tests/test_clinical_live.py::test_a_list_from_a_link_brings_the_survival_figure_into_view` | 19.4 | clinical_root, clinical_run; module: test_clinical_live.ran |
| 80.2 | call | gw4 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[made]` | 39.2 | made (by request) |
| 78.3 | call | gw10 | `tests/test_figures_m10.py::test_every_figure_redraws_byte_for_byte` | 38.1 | clinical_run |
| 76.0 | call | gw0 | `tests/test_figures_m9.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte` | 40.9 | built |
| 75.8 | call | gw9 | `tests/test_figures_m7.py::test_every_figure_has_a_spec_and_redraws_byte_for_byte[study_root]` | 30.2 | study_root (by request) |
| 74.6 | call | gw0 | `tests/test_ncrna_tcga.py::test_m5_takes_a_gdc_sheet_through_chosen` | 19.8 | none; module: test_ncrna_tcga.root |
| 71.5 | call | gw7 | `tests/test_process_count.py::test_a_moved_root_and_a_reordered_sheet_rerun_clean` | 9.5 | none; module: test_process_count.run |
| 71.2 | call | gw6 | `tests/test_qc_run.py::test_fetchngs_reads_alone_within_the_track_a_budget` | 29.0 | none |
| 70.2 | setup | gw0 | `tests/test_export_browser.py::test_buttons_give_the_files_at_the_defaults_and_say_when_they_draw[linked]` | 48.6 | none; module: test_export_browser.built, test_export_browser.ran |
| 69.5 | call | gw9 | `tests/test_genecard_pack.py::test_a_card_past_the_cap_is_named_with_its_command_and_writes_no_files` | 21.7 | none |
| 67.0 | setup | gw8 | `tests/test_quant_salmon.py::test_laptop_page_and_reuse` | 12.6 | none; module: test_quant_salmon.laptop |
| 66.9 | call | gw10 | `tests/test_figures_m2.py::test_a_run_with_its_own_registry_names_that_registry` | 31.8 | none |
| 66.9 | setup | gw3 | `tests/test_figures_m7.py::test_a_run_with_other_options_redraws_with_them` | 50.5 | across_root, proteo_across |
| 66.9 | setup | gw8 | `tests/test_box_stack_browser.py::test_the_page_draws_every_name_and_shown_label_apart[1280-11]` | 20.2 | none; module: test_box_stack_browser.ran |
| 65.4 | setup | gw11 | `tests/test_page_cap.py::test_m2_drops_its_flag_view_under_a_cap_one_byte_below_its_page` | 67.1 | page |

## Module fixtures built more than once

A module fixture is built once per worker and per file that requests it.

| Fixture | Builds | Files that request it, with the workers |
| --- | --- | --- |

Module fixtures built once: 70.

## Files by summed setup and call time

| File | Summed wall (s) | CPU (s) |
| --- | --- | --- |
| `tests/test_figures_m2.py` | 698.6 | 316.9 |
| `tests/test_cli.py` | 530.6 | 116.1 |
| `tests/test_figures_m7.py` | 384.9 | 168.0 |
| `tests/test_handoff_browser.py` | 377.5 | 160.1 |
| `tests/test_de_views_browser.py` | 371.0 | 151.3 |
| `tests/test_figures_m8.py` | 355.2 | 105.2 |
| `tests/test_function_views_browser.py` | 333.3 | 19.7 |
| `tests/test_qc_run.py` | 282.0 | 118.7 |
| `tests/test_figures.py` | 265.9 | 111.2 |
| `tests/test_figures_m9.py` | 250.9 | 108.6 |
| `tests/test_process_count.py` | 248.7 | 34.9 |
| `tests/test_page_cap.py` | 240.0 | 106.9 |
| `tests/test_pipeline.py` | 209.3 | 116.1 |
| `tests/test_process_views_browser.py` | 201.5 | 79.6 |
| `tests/test_process_align.py` | 200.9 | 33.8 |
| `tests/test_gather_counts.py` | 185.9 | 53.8 |
| `tests/test_figures_m10.py` | 181.3 | 75.1 |
| `tests/test_kit.py` | 147.3 | 15.7 |
| `tests/test_ncrna_views_browser.py` | 146.7 | 38.8 |
| `tests/test_proteo_views_browser.py` | 146.0 | 29.3 |
