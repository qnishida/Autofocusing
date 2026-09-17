# Peak-search handoff and validation (2026-09-17)

The agreed initial-seed selection changes are implemented on
`test/peak-search-boundaries`, based on main `485d761`. The implementation,
regression tests and these validation notes are included in the test-branch
change set; no installed binary, main branch or remote was updated. The earlier bias fix
is already on main. No covariance or convergence change is included here.

This report records the peak-only comparison. The subsequent
[covariance precision correction](parameter-covariance-precision-20260917_en.md)
addresses the NaN found below; the original failed-run records are retained.

## Preserved behavior and correction

Candidate centers remain inside the configured search circle. Background and
local comparisons use the computed square grid, including grid points outside
that circle. The empirical peak-minus-background threshold of 7 × MAD, the
approximately 0.01 s/km neighborhood, candidate limit, and the 0.1 s/km exclusion
for vertical-derived horizontal seeds are preserved. Frequency-dependent tuning
of these widths is outside this change.

- Rejecting a non-strict local maximum now skips that candidate rather than
  ending the entire search. The valid global threshold early exit remains.
- Background reads check both coordinates; the average divides by the actual
  number of contributing references. Interior accumulation order and duplicate
  ring endpoint weights are retained.
- Local-maximum comparisons count the valid neighbors actually visited,
  excluding the center. Missing neighbors do not disqualify edge seeds.
- The former second background loop had no iterations because
  `dpy_ary[imask] == 0`; removing it does not change interior sampling.

## Synthetic and regression checks

25 synthetic checks pass, covering central peaks in all three components,
capacity limits, plateaus followed by a separate peak, four edges and nearby
positions, exclusion of outside-circle centers while still using their grid
values as neighbors, actual boundary background counts, threshold cases,
nearby sidelobe suppression, empty detections, ±0.4 range and a finer dp.
Boost array bounds checks are enabled in this test translation unit.

The old production function fails the plateau test (misses the isolated peak)
and aborts on the boundary-only case at a Boost bounds assertion. The corrected
code passes. All 11 selected CTests pass after rebuilding affected targets.

## Real-data comparison

CPU, three components, March 15–17, 2005; 0.1–0.25 Hz requested, ±0.165 s/km,
dp=0.005 s/km, eight threads. Paired runs use identical start/end dates, waveform
and source-catalog hashes, accepted windows, fixed bootstrap seeds and FFTW
ESTIMATE | UNALIGNED. Baseline output reproduces the prior 60-event CPU catalog
byte-for-byte. Candidates are paired by time window, component and original
slowness seed, not by approximate final location.

| Measure | Baseline | Corrected |
|---|---:|---:|
| Initial candidates | 60 | 147 |
| Converged / emitted rows | 60 | 143 |
| Iteration-limit failures | 0 | 3 |
| Line-search failures | 0 | 1 |
| Selected R / T / U | 28 / 10 / 15 | 31 / 10 / 34 |

All 60 existing seeds and events are retained, with identical final parameters,
iteration counts and all 38 output columns, including matrices and energy ratios.
There are 87 additional seeds: 83 converge, three reach the iteration limit and
one fails line search. Of the additions, 3 R and 19 U pass the existing downstream
max/MAD thresholds (>7 for R/T, >35 for U). These are threshold-selected rows,
not independently confirmed new physical events.

None of the added seeds has an incomplete grid neighborhood in this sample
(maximum absolute coordinate 0.144 s/km). Thus these real-data additions exercise
continued search and propagation of vertical seeds; the edge policy is qualified
by the synthetic tests. Increasing the candidate count alone is not a failure.

## Failed finite-output check: downstream covariance

The real-data audit is **not an unconditional pass**. One added T event has NaN
in all ten parameter-covariance columns (21–30). Its seed is (0.063,-0.054) s/km
in the window beginning 2005-03-16 11:51:28. Its max/MAD is -2.70888, so it does
not pass the existing >7 selection. All 75 threshold-selected rows are finite.

A targeted two-day replay preserving QC history reproduces that row exactly.
The Hessian/sigma input is finite, but the single-precision determinant and
inverse are nonfinite. The double-precision determinant is approximately
5.20497e39, exceeding float's maximum approximately 3.40282e38; the double inverse
is finite. This identifies a numerical range problem in the shared CPU
`est_dist_boot()` inverse, separate from the previously clarified covariance
output-scale convention. It was exposed by a newly admitted seed; the peak fix
does not modify covariance arithmetic. A finite double inverse alone does not
establish statistical validity of the resulting uncertainty.

The runner's all-finite check remains failed and is not relaxed. The comparison
tool can diagnose this saved failure without reporting it as a pass. Before
production adoption, handle or explicitly decide how to treat nonfinite
covariance for weak candidates in a separate change. Keep the agreed seed
selection criteria rather than tightening them to hide this example.

## Reproduction and local records

```bash
cmake -S . -B build-clang -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-clang --target test_peak_search --parallel 4
ctest --test-dir build-clang -R '^peak_search_boundaries$' --output-on-failure
```

`tests/build_peak_search_probe.py` reuses the local deterministic audit controls;
see its prerequisites. After the covariance change, pass `--source` with the saved
`docs/benchmarks/peak-search-boundaries-20260917/peak-only-production.cpp` to
reproduce this historical comparison. Its hash matches the original provenance.
`tests/compare_peak_search_runs.py` compares complete
baseline/fixed runs and records added/lost seeds and failed audit checks.
Detailed logs, catalogs, source/build hashes, before-test failures and covariance
replay are retained locally in `docs/benchmarks/peak-search-boundaries-20260917/`,
which is ignored by Git. They are not intended for publication.
