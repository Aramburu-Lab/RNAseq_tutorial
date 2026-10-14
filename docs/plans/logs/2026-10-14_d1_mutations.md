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

