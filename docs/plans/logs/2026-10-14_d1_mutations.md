## CAUGHT: src/ctx/report/templates/_page.js.j2
- planted: '      var seen = heightsSet.some(function (h) { return Math.abs(h - height) <= 1; });' -> '      var seen = heightsSet.length > 0 && Math.abs(heightsSet[heightsSet.length - 1] - height) <= 1;'
- ran: tests/test_heatmap_cycle.py
- exit 1
```
FAILED tests/test_heatmap_cycle.py::test_a_height_cycle_settles_square_at_every_width_in_both_themes
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 71.68s (0:01:11)
```

## CAUGHT: src/ctx/viz/heatmap.py
- planted: '        return unlabelled or self._pitch_mm() >= HEAT["row_mm"] * self._row_lines()' -> '        return unlabelled'
- ran: tests/test_quant_figures.py
- exit 1
```
FAILED tests/test_quant_figures.py::test_labelled_rows_are_square_where_a_column_holds_their_label
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed, 3 passed in 2.00s
```

## CAUGHT: src/ctx/viz/heatmap.py
- planted: '        return unlabelled or self._pitch_mm() >= HEAT["row_mm"] * self._row_lines()' -> '        return True'
- ran: tests/test_quant_figures.py
- exit 1
```
FAILED tests/test_quant_figures.py::test_labelled_rows_are_square_where_a_column_holds_their_label
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed, 3 passed in 1.82s
```

## CAUGHT: src/ctx/de/deseq.py
- planted: '    return library_replicates(libraries, dict(order)) == 1' -> '    return rank_replicates(libraries, REPLICATE_KEYS, dict(order)) == 1'
- ran: tests/test_de.py::test_a_librarys_lanes_are_one_sample_and_a_second_library_is_dropped
- exit 1
```
FAILED tests/test_de.py::test_a_librarys_lanes_are_one_sample_and_a_second_library_is_dropped
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 4.98s
```

## CAUGHT: src/ctx/de/run.py
- planted: '        column = "log2FoldChange" if what == "lfc" else "padj"' -> '        column = "padj"'
- ran: tests/test_de.py::test_a_rule_names_a_gene_and_m5_records_its_fold_change_and_padj
- exit 1
```
FAILED tests/test_de.py::test_a_rule_names_a_gene_and_m5_records_its_fold_change_and_padj
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 14.06s
```

## CAUGHT: src/ctx/de/run.py
- planted: 'once.unique("gene_name", keep="first", maintain_order=True)' -> 'once.unique("gene_name", keep="last", maintain_order=True)'
- ran: tests/test_de.py::test_a_rule_names_a_gene_and_m5_records_its_fold_change_and_padj
- exit 1
```
FAILED tests/test_de.py::test_a_rule_names_a_gene_and_m5_records_its_fold_change_and_padj
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 13.23s
```

## CAUGHT: src/ctx/de/page.py
- planted: '    if words is None:\n        return level.replace("_", " ")' -> '    if True:\n        return level'
- ran: tests/test_de_words.py::test_every_step1_contrast_reads_in_plain_words
- exit 1
```
FAILED tests/test_de_words.py::test_every_step1_contrast_reads_in_plain_words
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 12.03s
```

## CAUGHT: src/ctx/de/page.py
- planted: '        if seen[found[c.name].label] > 1 and not c.label else found[c.name]' -> '        if False else found[c.name]'
- ran: tests/test_de_words.py::test_contrasts_with_no_label_that_read_alike_are_told_apart
- exit 1
```
FAILED tests/test_de_words.py::test_contrasts_with_no_label_that_read_alike_are_told_apart
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 1.15s
```

## CAUGHT: src/ctx/de/run.py
- planted: 'trend = word(PAGE, f"figures.pca.trend.{result.trend}")' -> 'trend = word(PAGE, "figures.pca.trend.parametric")'
- ran: tests/test_de_words.py::test_the_pca_caption_names_the_dispersion_trend_the_vst_used
- exit 1
```
FAILED tests/test_de_words.py::test_the_pca_caption_names_the_dispersion_trend_the_vst_used
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 31.58s
```

## CAUGHT: src/ctx/function/run.py
- planted: '        *_acknowledged(spec),\n' -> ''
- ran: tests/test_de_words.py::test_m6_of_gdc_data_acknowledges_the_project_and_the_gdc
- exit 1
```
FAILED tests/test_de_words.py::test_m6_of_gdc_data_acknowledges_the_project_and_the_gdc
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 25.30s
```

## CAUGHT: src/ctx/report/builder.py
- planted: '        required = ("title", *(s for s, _ in SECTIONS))' -> '        required = ("title", *(s for s, _ in SECTIONS if s != "input"))'
- ran: tests/test_report.py::test_skeleton_order
- exit 1
```
FAILED tests/test_report.py::test_skeleton_order - Failed: DID NOT RAISE Valu...
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 3.13s
```

## NOT CAUGHT: src/ctx/report/prose.py
- planted: '        items = list(found.values()) if isinstance(found, dict) else found' -> '        items = list(found) if isinstance(found, dict) else found'
- ran: tests/test_prose.py
- exit 0
```
.....................................                                    [100%]
37 passed in 1.71s
```

## CAUGHT: src/ctx/report/prose.py
- planted: '        items = list(found.values()) if isinstance(found, dict) else found' -> '        items = list(found) if isinstance(found, dict) else found'
- ran: tests/test_prose.py
- exit 1
```
FAILED tests/test_prose.py::test_a_table_of_question_variants_holds_each_prompt_to_its_budget
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed, 9 passed in 0.33s
```

## CAUGHT: src/ctx/de/deseq.py
- planted: '    return library_replicates(libraries, dict(order)) == 1' -> '    return rank_replicates(libraries, REPLICATE_KEYS, dict(order)) == 1'
- ran: tests/test_de.py::test_a_librarys_lanes_are_one_sample_and_a_second_library_is_dropped
- exit 1
```
FAILED tests/test_de.py::test_a_librarys_lanes_are_one_sample_and_a_second_library_is_dropped
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 6.11s
```

## CAUGHT: src/ctx/viz/strip.py
- planted: '        refit_y_title(fig, ax, self._value_title())' -> '        pass'
- ran: tests/test_title_room.py
- exit 1
```
FAILED tests/test_title_room.py::test_a_vertical_title_fits_its_own_drawn_axis[strip-long title]
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed, 21 passed in 5.39s
```

## CAUGHT: src/ctx/viz/scatter.py
- planted: '        refit_y_title(fig, ax, self._axis_title(self.y_label, self.y))' -> '        pass'
- ran: tests/test_title_room.py
- exit 1
```
FAILED tests/test_title_room.py::test_a_vertical_title_fits_its_own_drawn_axis[volcano-long title]
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed, 1 passed in 2.04s
```

## CAUGHT: src/ctx/viz/axes.py
- planted: '    return "<br>".join(wrap_title(title, room).split("\\n"))' -> '    return "<br>".join(fit_title(title, room / 72 / 0.55).split("\\n"))'
- ran: tests/test_page_axes.py
- exit 1
```
FAILED tests/test_page_axes.py::test_the_page_y_title_is_the_files[volcano-as built]
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed in 9.64s
```

Note: the axes.py mutation above failed by a NameError (fit_title no longer imported there), so it proves nothing; the scatter one below replaces it.
## CAUGHT: src/ctx/viz/scatter.py
- planted: 'page_y_title(y_title, files_room(self.static))' -> 'page_y_title(y_title, self.inches()[1] * 72 * 0.55)'
- ran: tests/test_page_axes.py
- exit 1
```
FAILED tests/test_page_axes.py::test_the_page_y_title_is_the_files[volcano-long title]
!!!!!!!!!!!!!!!!!!!!!!!!!! stopping after 1 failures !!!!!!!!!!!!!!!!!!!!!!!!!!!
1 failed, 1 passed in 2.36s
```

