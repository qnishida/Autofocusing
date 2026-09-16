# Estimate distributions and signed shifts — 2005 matrix audit

[Japanese](matrix-estimate-distributions_jp.md)

This uses the expanded 48-day sample from the [thresholded matrix audit](matrix-ratio-threshold-2005_en.md):
171 R, 36 T and 151 U/V selected detections. It describes estimated energy
ratios themselves, rather than absolute percentage differences.

The spread is **between catalog events**. It mixes real differences in the
wavefield with estimation variability. It is not the uncertainty obtained by
repeatedly estimating the same event. Original-minus-corrected paired differences
isolate the effect of the matrix formula at the same fitted points. The corrected
formula is not known physical truth: bias against truth, including fitting or
selection bias, cannot be identified from these catalogs alone.

## Ratio distributions

Entries are median [25th, 75th percentile], in dimensionless ratio units.
Only the same events with both powers positive in both formulas are compared.

| Detected component | Ratio | N | Original estimate | Corrected estimate |
| --- | --- | ---: | --- | --- |
| R | R/U | 171 | 0.580 [0.406, 1.011] | 0.576 [0.394, 1.041] |
| U | R/U | 150 | 0.188 [0.126, 0.462] | 0.201 [0.128, 0.456] |
| T | T/R | 36 | 0.983 [0.762, 1.249] | 1.038 [0.680, 1.505] |
| T | T/U | 32 | 1.976 [1.240, 3.285] | 2.648 [1.524, 4.723] |

For example, corrected T/R has median 1.038 and its middle half spans
0.680–1.505. Its mean and sample SD are 1.140 and 0.564. Corrected T/U has
median 2.648, middle half 1.524–4.723, mean 3.803 and SD 4.204. These SDs
describe event-to-event spread, not a confidence interval or standard error.

## Direction of the implementation shift

Define `delta = original_ratio - corrected_ratio` for each matched event.
A positive value means the original estimate was higher. These are differences
in dimensionless ratio units, not percent errors. The median paired difference
need not equal the difference between the two marginal medians above.

| Detected component | Ratio | Mean delta | Median delta | Original higher / lower |
| --- | --- | ---: | ---: | ---: |
| R | R/U | -0.1475 | -0.0031 | 64 / 107 |
| U | R/U | +0.0005 | +0.0001 | 75 / 75 |
| T | T/R | +0.1128 | +0.0624 | 19 / 17 |
| T | T/U | +3.4385 | -0.4314 | 4 / 28 |

T/R has mixed signs (19 higher, 17 lower), although its mean paired shift is
positive. T/U is lower originally in 28/32 events. Its positive mean delta is
dominated by an outlier: original T/U 147.35 versus corrected 5.58, a difference
of 141.77. Thus a positive mean alone would misrepresent the prevailing direction
of the T/U shifts. R/U in U detections is balanced in this sample (75/75).

## Three-component shares in T detections

For the 32 T detections with all three powers positive in both formulas,
normalize each diagonal by `SRR + STT + SUU`. Entries show mean ± sample SD
across events; mean shifts are percentage points (original minus corrected).

| Share | Original mean ± SD | Corrected mean ± SD | Mean shift, points | Original higher / lower |
| --- | --- | --- | ---: | ---: |
| R | 38.84 ± 9.40% | 40.93 ± 10.98% | -2.09 | 10 / 22 |
| T | 39.66 ± 12.37% | 40.97 ± 11.95% | -1.30 | 18 / 14 |
| U | 21.50 ± 13.16% | 18.11 ± 12.34% | +3.39 | 29 / 3 |

The original U share is higher in 29/32 events, with a mean shift of +3.39
points and median shift +3.00 points. This gives a clearer directional effect
than unstable ratios with small U denominators. It is conditional on these
selected positive-power events, and is not a universal correction or proof of
bias relative to physical truth. Neither the SD nor the number of rows assumes
independent earthquakes; selected T detections occupy 19 calendar dates.

For repeated-estimation uncertainty and statistical bias, the next distinct
experiment would need resampling of the same event, and known-input simulations
(or another specified truth reference) to evaluate bias. The current catalogs
alone do not provide those quantities for all matrix entries and ratios.

[Complete distributions, paired shifts, row IDs and input/script hashes](benchmarks/matrix-ratio-threshold-2005-20260916/estimate-distributions.json)
are reproducible with:

```bash
python3 -B tests/summarize_matrix_estimates.py   docs/benchmarks/matrix-ratio-threshold-2005-20260916/expanded-thresholds.json   --output build-clang/estimate-distributions.json
```

[Historical catalog agreement and input-identity limits](historical-catalog-agreement_en.md).
