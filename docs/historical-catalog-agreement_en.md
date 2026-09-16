# Historical catalog agreement — 48-day 2005 audit

[Japanese](historical-catalog-agreement_jp.md)

Most strong detections have close numerical counterparts in the supplied
catalog, but the catalogs are not identical. This comparison uses historical
`../analysis/2004_2048_0.099609-0.250000.dat` and the current **original-formula**
outputs, before the candidate self-term correction. The historical catalog
selects dates; the [matrix correction comparison](matrix-ratio-threshold-2005_en.md)
itself uses identical current waveforms and fitted events for both formulas.

## Counts and matching

On the 48 calculated days, historical/current counts are 666/677 overall.
After R/T max/MAD > 7 and U > 35, counts are R 167/171, T 37/36, U 149/151.

Matching requires equal start time, end time and detected component. Within
each group, a one-to-one assignment first maximizes the number of matches
within 0.005 s/km Euclidean distance in (px, py), then minimizes total slowness
distance. This tolerance is one search-grid step. Power values are not used
in matching. It is an operational correspondence, not proof of event identity.

There are 604 matched pairs, 62 unmatched historical rows and 73 unmatched
current rows. Of the 358 current selected rows, 348 (97.2%) match. The unmatched
sets include changed detections and candidates outside the matching tolerance;
they cannot all be labeled new or missing physical events.

For matched current selected events, the table gives median absolute differences.
Relative differences use the current original-formula value as denominator.
`MS_max` is the scalar fitted-component power in column 33.

| Component | Matched/current selected | Slowness distance, s/km | max/MAD relative difference | MS_max relative difference |
| --- | ---: | ---: | ---: | ---: |
| R | 164/171 | 3.09e-06 | 0.134% | 0.0061% |
| T | 35/36 | 3.29e-06 | 0.141% | 0.0093% |
| U | 149/151 | 1.78e-06 | 0.092% | 0.0036% |

Typical slowness differences are approximately 1/1500–1/2800 of a grid step.
Some pairs are less close: maximum MS_max differences in selected R/T/U pairs
are 13.71/11.46/6.41%, and source-distance differences can also be large.
Thus small medians do not imply that every fitted solution agrees.

Tolerance sensitivity: at 0.0001 s/km (1/50 grid step), 325/358 selected rows
match (R 150, T 32, U 143); at 0.01 s/km, 351/358 match (R 166, T 36, U 149).
This checks how much the matching rate depends on the tolerance. No sign reversal
of slowness or waveform-orbit equivalence is silently applied.

## Input identity and other columns

All 604 matched rows have identical printed array-center coordinates.
Among the 348 matched selected rows, accepted-window counts are identical in
318; the historical count is exactly one larger in the other 30. The subsequent
[consistency investigation](catalog-consistency-investigation_en.md) identifies
restart-dependent QC history as a cause and verifies it by a preceding-day rerun.
Current input SHA-256 values remain unchanged before and after each run; the
historical waveform hashes and full run configuration are unavailable, so
identity of the historical input files cannot be established.

The 38-column layout matches, but not every field reproduces. For all 348
matched selected historical rows, column 9 equals hypot(px,py) to 1e-6 and
column 10 equals minus source distance. Current code fills these columns with
slowness and distance changes from the initial fit parameters (`prm.dp` and
`-prm.dD` in `output_result`). The user identifies the public CPU source as the historical basis. Public code
had an uninitialized `prm_init`; commit `393a93e` fixed it. This known diagnostics
fix is consistent with the observed column differences; see the investigation. Diagnostics, iteration counts
and bootstrap outputs should not be assumed identical merely from the layout.

- [Main comparison, pair IDs, unmatched rows, metrics and hashes](benchmarks/matrix-ratio-threshold-2005-20260916/historical-comparison.json)
- [Tight tolerance](benchmarks/matrix-ratio-threshold-2005-20260916/historical-comparison-tight.json) and [wide tolerance](benchmarks/matrix-ratio-threshold-2005-20260916/historical-comparison-wide.json)
- [Estimate distributions and signed implementation shifts](matrix-estimate-distributions_en.md)

Reproduce with `tests/compare_historical_catalog.py HISTORICAL EXPANDED_JSON
--output REPORT_JSON`; the default matching distance is 0.005 s/km.
