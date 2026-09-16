# Three-component spectral-matrix validation, 2005 data — 2026-09-16

Detailed run logs, output catalogs and per-event comparison records are retained
locally and are not distributed in this repository. References to these records
below describe local artifacts.

[Japanese](spectral-matrix-3c-2005-verification_jp.md)

For the later R/T > 7 and U > 35 max/MAD selection, January extension, and
catalog-guided March dates, see the [thresholded ratio follow-up](matrix-ratio-threshold-2005_en.md).

The existing self-term indexing error affects matrix output materially on
approximately 685-station real data. A candidate correction on
`test/spectral-matrix-3c-2005` passes independent matrix tests and preserves
all non-matrix results in the three-day comparison.

## Conditions

- Reference source: `9ca6a93`; corrected source differs only in `cal_S_matrix`
  numerical operations. Each complex same-station term is subtracted from its
  matching matrix entry, retaining its imaginary part.
- Input: `/Volumes/Seismic_Data/hdf5/Hi-net/2005`, January 1–3, 2005, three
  components. The daily loader reports 687, 685, and 684 stations. Per-window
  QC can reduce these counts. Selected files were available while the larger
  archive was still being copied; their SHA-256 hashes matched before and
  after both runs.
- Requested frequency band: 0.1–0.25 Hz. Existing floor-to-bin conversion
  gives bins 102–256, or 0.099609375–0.25 Hz, with 155 bins.
- Slowness: px/py ±0.165 s/km, step 0.005 s/km, 67-by-67 grid; existing radial
  candidate mask retained.
- CPU backend, GPU power off, four OpenMP threads, dynamic teams disabled.
  No performance claim: other analysis and data copying may be running.
- Test-only drivers fix FFTW to `ESTIMATE | UNALIGNED` and bootstrap seeds to
  `1837 + 104729 * replicate`. Production FFTW/seed behavior is unchanged.
- One bounded run per version; each uses the same catalog and three input
  files. Catalog hashes also match before and after the comparison.

## Independent matrix check

`tests/spectral_matrix.cpp` explicitly sums only pairs of distinct stations,
independently of the production beam-sum-minus-self-term algorithm. It checks
all nine complex entries, Hermitian symmetry, and agreement of each diagonal
with the scalar `cal_S(..., flag_red=1)` calculation. Cases use nonuniform
weights, component masks, two windows, and correlated complex components:
7 stations with all 155 default-band bins, and 700 stations with five bins.
Both three-component and horizontal modes are tested, including missing-U NaNs.
Zero delays isolate matrix accumulation from travel-time modeling.

The unmodified code fails the first case: RR is 0.461623 instead of 1.04686.
The corrected code passes all cases. All 10 selected CTests pass, including
frequency, CPU fitting, I/O, horizontal components, slant stack, workspace,
launcher, and date-range regressions.

## Real-data comparison

Both versions finish normally with 45 events (18 R, 10 T, 17 U), 12 segments,
and 450 accepted windows. No output values are non-finite. No loader errors
or missing input days were found.

All first 33 output columns match exactly, including location, slowness,
parameter covariance, and bootstrap power. Initial grid candidates and
accepted-window counts also match exactly. Thus this correction does not
change detection or fitted parameters in these data.

For each event's fitted component, compare the original diagonal with the
independently corrected scalar reference `MS_max/(155*df)`, using output
column 33 for `MS_max`. Absolute relative differences are:

| Fitted component | Events | Median | 95th percentile | Maximum |
| --- | ---: | ---: | ---: | ---: |
| R | 18 | 5.25% | 18.66% | 28.55% |
| T | 10 | 12.22% | 34.41% | 34.75% |
| U | 17 | 5.68% | 9.95% | 12.47% |

The corrected matrix diagonal matches this scalar reference with a maximum
scaled difference of `4.4421e-6`, within the six-significant-digit output
precision (check threshold `2e-5`).

These relative differences compare two calculations for the same event;
they do not measure event-to-event variation of SUU or statistical error bars.
Directly comparing SUU before/after for the 17 U-detected events, all corrected
SUU values are positive. The original value is higher in 15 cases and lower
in two. The signed relative mean is +5.46% and median is +5.68%, indicating
an upward shift in this sample. This does not establish a universal positive
bias: the extra real cross-component self terms can have either sign.

Across all 45 events, including components other than the fitted component,
median absolute relative changes are 6.82% for RR, 33.64% for TT, 8.93% for UU,
7.08% for Re(RU), and 137.33% for Im(RU). These percentages divide by the
corrected entry's magnitude; weak or near-zero entries can produce very large
ratios. They are not percentage errors in polarization angle or ellipticity.

For SUU across all 45 events, 38 differences are positive and seven negative;
eight corrected SUU values are negative. The relative-error definition uses
the absolute corrected entry in the denominator, so the all-event statistics
must not be interpreted as ratios of two positive powers in every case.

## Diagonal power ratios

For a ratio `a/b`, the reported implementation difference is
`100*((old_a/old_b)/(corrected_a/corrected_b)-1)` for the same event.
Only pairs with both diagonal values strictly positive in both versions are
included. No negative entries were clipped, converted to magnitudes, or
replaced by an arbitrary positive floor.

| Ratio | Included events | Median old ratio | Median corrected ratio | Median absolute relative difference |
| --- | ---: | ---: | ---: | ---: |
| RR/UU | 34 | 0.7491 | 0.7620 | 4.99% |
| TT/UU | 11 | 0.1103 | 0.07042 | 60.51% |
| TT/RR | 16 | 0.3503 | 0.2573 | 34.83% |

The medians of the ratios are separate distribution summaries; their ratio
is not the median event-wise relative difference. Signs vary across events:
the signed relative medians are -2.21%, -3.46%, and -5.43%, respectively.
The mixed signs and outliers do not support applying a single correction
factor to all events. Similar component power errors do not necessarily
cancel in a ratio, especially when one component is weak.

Only 10 events have all three powers positive in both versions. In this
common subset, the means of each event's normalized diagonal fractions
`S_aa/(S_RR+S_TT+S_UU)` change as follows:

| Component | Old mean share | Corrected mean share | Median absolute event-wise change |
| --- | ---: | ---: | ---: |
| R | 34.01% | 34.85% | 1.20 percentage points |
| T | 7.58% | 7.77% | 0.75 percentage points |
| U | 58.42% | 57.38% | 0.81 percentage points |

These are means of per-event fractions, not pooled energy fractions or
an estimate for the whole wavefield. Individual changes can be much larger:
the maximum absolute change in U's share is 18.69 percentage points.
Across all 45 events, corrected RR/TT/UU are nonpositive in 3/24/8 events.
Self-term-subtracted coherent beam powers are not guaranteed positive;
such entries cannot directly be interpreted as positive energy ratios.
Even a small positive entry can give an unstable ratio. This analysis does
not add a statistical significance threshold for component power.

Per-event ratios, inclusion rules, signs, and fraction statistics
are saved with the comparison artifacts.

## Reproduction and records

Use the existing `build-clang` setup with Homebrew LLVM C/C++ and Apple Clang
Objective-C++. Build targets before generating the corrected driver:

```bash
cmake -S . -B build-clang
cmake --build build-clang --target cal_ccf_clang test_spectral_matrix -j2
python3 tests/build_io_probe.py build-clang --reference 9ca6a93
python3 tests/build_io_probe.py build-clang
OMP_NUM_THREADS=4 ctest --test-dir build-clang -R spectral_matrix_station_pairs --output-on-failure
python3 -B tests/run_spectral_matrix_events.py \
  build-clang/io-fixed-reference/cal_ccf_io_probe \
  build-clang/io-fixed-current/cal_ccf_io_probe \
  /Volumes/Seismic_Data/hdf5/Hi-net ../moment_loc_76_24 \
  build-clang/spectral-matrix-3c-2005-new --threads 4 --start 2005-01-01 --days 3
```

Use a new output directory. The runner rejects changing inputs and compares
output fields without relaxing non-matrix equality requirements.

- Measurements, checks, commands, and input/binary hashes
- Reference events and corrected events
- Reference log and corrected log
- Original matrix-test failure, 10 passing CTests, and artifact checks
- Source/build provenance and source change

The evidence qualifies this matrix correction on these three days; it does
not validate all modeling assumptions, the full archive, GPU arithmetic, or
the separate slowness-selection issues. Correctly self-term-subtracted
matrices need not be positive semidefinite for finite data.
These runs were performed on the test branch without replacing the installed
binary, merging, or pushing. The adoption decision is recorded below.

## Interpretation of the real-only subtraction (2026-09-17 clarification)

The original author's intent is not established. A common instrument-noise
model with real channel couplings and no relative phase can have a real
cross-component noise covariance; discarding its imaginary part may therefore
have been intentional. Common origin alone does not imply zero relative phase.
The measured same-station product used here includes signal as well as noise,
so it is not independently an estimate of instrument noise.

The candidate correction implements removal of all same-station pairs,
consistent with the existing denominator sum(i != j) w_ai w_bj. Equal unit
weights give N(N-1) per window. Under that estimator, the full complex
same-station term must be removed. An instrument-noise-only estimator is a
different modeling choice and needs a noise estimate and matching normalization.
Dividing the original sum over component pairs by a component count alone does
not generally repair its indexing or unequal weights/couplings.

Once the destination index is corrected, retaining versus discarding the
imaginary part affects off-diagonal entries only: same-component products
are real. Thus the real-arithmetic change in SRR, STT and SUU is attributable
to the destination-index correction, not to retaining imaginary self terms.
This last separation follows algebraically; a separate index-only run was
not performed. No production estimator was changed in response to this
clarification.

## Adopted correction (2026-09-17)

Remove the full complex same-station contribution from its corresponding matrix
entry: `S_RTU(a,b) -= sum_i,k conj(X_ai) X_bi w_ai w_bi`. Both its real and
imaginary parts are subtracted, consistent with the existing distinct-station
pair normalization. The source and independent station-pair regression test
are included with this report. The previously tested source is unchanged.

The peak-search findings and convergence-condition change are separate work;
the covariance output convention is preserved. See the
[CPU/Metal comparison](cpu-gpu-3c-2005-verification_en.md) and the
commit-time validation record.
