# D1, 2026-10-14: the page script's height cycle and square labelled heatmaps

Measured in headless Chromium through `tests/browser.mjs`, on the M4 page of the laptop path
(`ctx quant --dataset tier1_teach`, 8 teaching runs, built into a scratch root from
`~/ctx_tier1`'s gather and process outputs; quant wall 79.30 s cold, 30.49 s warm by
`/usr/bin/time -l`).

## Before: `_square` true for labelled rows, the old guard (`heightSet !== height`)

A probe in `apply` recorded [held height, drawn height, margin t, margin b, plot w, div width]
at 1,280 px; the page's load event never fired (driver timed out at 60 s):

| call | held | drawn | t | b | w |
|---|---|---|---|---|---|
| 1 | 938.1 | 1101.8 | 192 | 159 | 427 |
| 2 | 946.5 | 938.1 | 221 | 159 | 412 |
| 3 | 938.1 | 946.5 | 192 | 159 | 427 |
| 4 | 946.5 | 938.1 | 221 | 159 | 412 |

and so on to the probe's limit of 30. Cause: plotly's top margin and plot width are a step
function of the height: at 946.5 px they hold the cells at 938.1 px, at 938.1 px at 946.5 px.
No height between is a fixed point, so no guard on a repeat can end it.

## After: `settle` stops on any height set before at this width, takes the cycle's tallest once

| page | window | theme | height (px) | heights drawn in 3 s after settling | cell height / width |
|---|---|---|---|---|---|
| M4 overview, real | 1280 | white, then black | 946.5 | 0, then 1 (the toggle's redraw) | 1.014 |
| M4 overview, real | 375 | black, then white | 787.6 | 0, then 1 | 1.017 |
| figure opened alone, real | 1280 | white, black | 1891.5 | 0, 1 | 1.000 |
| figure opened alone, real | 375 | black, white | 787.6 | 0, 1 | 1.017 |
| fixture twin, old script | 1280 | | page never loaded | | |
| fixture twin, new script | 1280 / 375 | both | 946.5 / 787.6 | 0 | 1.014 / 1.017 |

The fixture (`tests/fixtures/heatmap_cycle.csv`) is the real table's rho, cell types and
batches with the 8 samples renamed XRX0000001 to XRX0000008; it reproduces the real page's
heights and aspects to the tenth of a pixel. `tests/test_heatmap_cycle.py` holds it at both
widths, both themes, and after the window narrows and widens.
