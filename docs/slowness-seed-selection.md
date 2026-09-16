# Slowness range and initial seeds — 2026-09-16

Static inspection of `src/cal_ccf.cpp` confirms that increasing
`AUTOFOCUSING_SLOWNESS_MAX` expands both the slant-stack grid and the radial mask
used by `search_max`. There is no separate 0.165 s/km cap on the selected seeds.
For horizontal data, R and T peaks are selected independently, with at most 30
seeds per component. Each accepted `(px, py)` supplies
`prm.p = sqrt(px*px + py*py)` and `prm.θ = atan2(py, px)` to distance-grid search
and then gradient fitting. Distance-grid search retains the seed slowness while
searching distance and slowness gradient. Subsequent fitting can change slowness.

Candidate acceptance also requires contrast relative to MAD, a strict local
maximum in a neighborhood of approximately 0.01 s/km, and separation from
previously selected peaks. Candidates are visited in descending stack power,
not by slowness magnitude. Expanding the domain can change MAD and candidate
competition, so it does not guarantee selection of every large-slowness peak.

## Issues identified by inspection

- In `search_max`, the local background ring indexes `ssRTU` at offsets of up to
  `imask` without checking the array bounds. The array itself spans only
  `[-ipmax, ipmax]` in both directions. The circular candidate mask does not
  exclude all points whose neighborhoods cross a square-grid edge, particularly
  near the coordinate axes. Such candidates can cause out-of-bounds reads.
- If a candidate fails the strict local-maximum count, `break` terminates the
  entire candidate scan rather than skipping that candidate. This can prevent
  later candidates, including large-slowness candidates, from being examined.

These behaviors predate the date-range/workspace changes. No peak-selection
code was modified and no synthetic peak test or large-slowness waveform
validation was performed for this inspection. Before claiming reliable seed
selection across the full expanded domain, verify and address edge handling and
early termination with synthetic interior, boundary and competing peaks.
