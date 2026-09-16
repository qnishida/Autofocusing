# Spectral-matrix ratios after event selection — 2026-09-16

Detailed run logs, output catalogs and per-event comparison records are retained
locally and are not distributed in this repository. References to these records
below describe local artifacts.

[Japanese](matrix-ratio-threshold-2005_jp.md)

For historical catalog, current-before, and current-after consistency, see the [regression investigation](catalog-consistency-investigation_en.md).

For estimate distributions and signed shifts, see [the distribution follow-up](matrix-estimate-distributions_en.md).

This extends the [three-day matrix verification](spectral-matrix-3c-2005-verification_en.md).
The comparison is between the original same-station subtraction in `9ca6a93`
and the candidate correction on `test/spectral-matrix-3c-2005`, at the same
fitted events. It is not a CPU/GPU precision comparison.

## Selection and interpretation

Use component column 2 and `max/MAD` from columns 12/13 of the 38-column
catalog: strictly greater than 7 for R or T, and strictly greater than 35
for U (vertical, also called V). Selection applies to the detected component;
it does not establish significance of the other two powers in that event.

For each ratio, report `100 * (ratio_original / ratio_corrected - 1)`.
An absolute relative difference is the absolute value of this quantity.
Include an event only if both powers in that pair are positive in both
versions. Do not clip negative powers or replace small denominators.
The diagonal powers are columns 34–36: SRR, STT, SUU.
These differences are implementation effects, not statistical error bars
or event-to-event variability. A small positive component can cause a large
relative ratio difference. Normalized fractions use only the common subset
with all three powers positive in both versions.

## Catalog-guided dates

The supplied `../analysis/2004_2048_0.099609-0.250000.dat` actually contains
2004–2024. In its 2005 rows, 895 R, 37 T, and 892 U detections pass the
respective thresholds. Historical detections guide date selection only;
the new outputs are independently thresholded.

| Additional date | Historical T detections > 7 | Largest historical T max/MAD |
| --- | ---: | ---: |
| 2005-03-16 | 5 | 36.52 |
| 2005-03-17 | 4 | 46.01 |
| 2005-03-27 | 2 | 8.55 |
| 2005-03-28 | 1 | 9.77 |

The first pair supplies strong T events; the second supplies events closer
to the cutoff. All four daily waveform files exist. January 1–31 is retained
as a consecutive-period control, without duplicating dates. The selected
days do not provide an unbiased estimate of annual event populations or
independent earthquake counts. Catalog rows can describe related signals.

Catalog hash, selected dates, and all 37 historical T row references
preserve the selection without copying the full private catalog.

## Test method

The test-only audit driver evaluates both matrix formulas after fitting and
writes paired catalogs. On January 1–3, both catalogs match the independently
run original and corrected binaries byte for byte (45 events). Equality of
the first 33 columns in subsequent paired runs is by construction, not an
independent repeated fit. Production binaries are not replaced.

Settings: three components, CPU, GPU power off, requested 0.1–0.25 Hz
(bins 102–256, actual 0.099609375–0.25 Hz), px/py ±0.165 s/km, step
0.005 s/km, existing radial mask. Test-only deterministic FFTW plans and
bootstrap seeds follow the three-day verification. Each run records thread
count, station counts, source hashes, command, and checks. Input and source
catalog hashes are checked again after computation to detect copying changes.

Example, using fresh output directories and the existing documented build:

```bash
python3 tests/build_matrix_audit_probe.py build-clang
python3 -B tests/run_matrix_audit_events.py \
  build-clang/matrix-audit/cal_ccf_matrix_audit \
  /Volumes/Seismic_Data/hdf5/Hi-net ../moment_loc_76_24 \
  build-clang/matrix-audit-new --start 2005-03-16 --days 2 --threads 8
python3 -B tests/analyze_matrix_ratios.py build-clang/matrix-audit-new \
  --output build-clang/matrix-audit-new/ratios.json
```

The two `io-fixed-*` drivers must first be generated as described in the
three-day verification. `moment_loc_76_24` is the runtime source catalog;
the supplied spectral catalog above is used only to choose dates.

## Consecutive January control

January 1–31 produces 331 detections; 50 R, 4 T, and 49 U pass the thresholds.
Median absolute relative ratio differences are:

| Detected component | Ratio | Positive-pair events | Median absolute difference |
| --- | --- | ---: | ---: |
| R | R/U | 50 | 2.29% |
| U | R/U | 49 | 3.36% |
| T | T/R | 4 | 34.97% |
| T | T/U | 3 | 228.19% |

The four T detections occur on January 27–28. Their small sample motivates
the additional dates. In R/U detections, T may be weak or negative despite
passing the detected-component threshold; the threshold alone cannot make
all energy ratios reliable.

January per-event results and all ratio/fraction summaries
and test provenance
are retained with individual run logs and checks.

## Expanded sample: 48 days

After the initial 35-day evaluation, all remaining 2005 dates with historical
T max/MAD > 7 were added: February 12–13, March 15, April 2, August 9,
September 18, October 5–6, 9 and 28, November 4–5 and 25. These 13 days
add 276 detections, with 90 R, 20 T and 86 U passing selection. Each additional
run uses four OpenMP threads; at most two runs execute concurrently.

| Coverage | R selected | T selected | U selected |
| --- | ---: | ---: | ---: |
| Initial 35 days | 81 | 16 | 65 |
| Expanded 48 days | 171 | 36 | 151 |

The expanded sample has 677 total detections. Selected T detections occur on
19 distinct days, compared with six previously. All 19 historical T dates
are covered. The historical catalog has 37 selected T rows versus 36 in the
current outputs: January 28 has four historical versus three current rows;
all other per-day T counts match. This count comparison is not one-to-one
identification of fitted events, and exact reproduction of the historical
catalog is not established. The paired matrix formulas use identical fitted
events and selection fields.

Every run passes its checks; all 48 requested days load, daily station counts
range from 669 to 693, and dates do not overlap between runs. Archived event
hashes match the run reports. The largest scalar-diagonal scaled difference
is `6.453e-6`, below the `2e-5` output-rounding tolerance.

| Detected component | Ratio | Positive-pair events | Median absolute difference | 95th percentile | Maximum absolute difference |
| --- | --- | ---: | ---: | ---: | ---: |
| R | R/U | 171 | 2.35% | 16.75% | 71.85% |
| U | R/U | 150 | 3.35% | 19.03% | 54.33% |
| T | T/R | 36 | 15.00% | 44.30% | 297.70% |
| T | T/U | 32 | 22.99% | 152.49% | 2539.38% |
| T | R/U | 32 | 24.53% | 139.46% | 1543.25% |

The 95th percentile describes the event-wise difference distribution; it is
not a confidence limit. Only 150 of the 151 U detections enter R/U because
one has nonpositive R. Only 32 of the 36 T detections enter ratios involving
U because four have nonpositive U in at least one version.

The main conclusion persists with more events: selected T detections show
material T/R and T/U ratio changes. Their median absolute differences change
from 13.67% to 15.00% for T/R, and from 16.81% to 22.99% for T/U. Signed medians
are +8.73% and -19.44%, respectively; these do not support one universal
multiplicative correction. Large outliers remain, especially with weak U.
R/U in R or U detections has a smaller typical difference, but individual
outliers reach 71.85% and 54.33%.

For the 32 T detections with all powers positive in both versions, mean
per-event R/T/U shares change from 38.84/39.66/21.50% to 40.93/40.97/18.11%.
Median absolute per-event share changes are 3.47/2.69/3.38 percentage points;
maxima are 8.47/16.26/16.11 points. These are not pooled wavefield energies.
Corrected T is nonpositive in 120/171 selected R and 119/151 selected U events.
A threshold on the detected component does not establish significance of the
other powers. The sample is selected by historical T activity, and the 36 T
rows are not 36 independent earthquakes or an unbiased annual population.

- Expanded per-event results and full ratio/fraction statistics
- Additional 13 days only
- Added-date manifest and run directories
- Coverage, hashes, checks and historical-count comparison
- Initial 35-day results and initial targeted March results

The validation runs did not replace an installed binary, merge, or push.
The [adoption decision](spectral-matrix-3c-2005-verification_en.md#adopted-correction-2026-09-17)
records the selected complex self-term subtraction. A subsequent [CPU/Metal comparison](cpu-gpu-3c-2005-verification_en.md)
passed for three days and 60 real three-component events.
