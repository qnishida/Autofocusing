# Estimate distributions and paired matrix shifts

The spread of estimated energy ratios across events combines variation in the
wavefield with estimation variability. It is not repeated-estimation uncertainty
for a single event, and it does not identify bias relative to physical truth.

## Ratio distributions

Report medians and quartiles in dimensionless ratio units. Compare identical
paired events, requiring positive numerator and denominator in both formulas;
report the excluded counts. Small denominators can make means and standard
deviations unstable. Between-event standard deviation is not a confidence interval.

## Direction of the implementation shift

Define a paired difference explicitly, for example
`delta = original_ratio - corrected_ratio`. Report its sign, median, mean,
and counts of positive/negative differences. The median paired difference
need not equal the difference between marginal medians. Inspect outliers before
using a mean to describe the prevailing direction of an implementation change.

## Three-component shares

For events with all three powers positive, diagonal shares can be normalized by
`SRR + STT + SUU`. Report share differences in percentage points and state the
conditioning sample. These shares are not a universal correction or evidence
of accuracy relative to truth. Rows from related events need not be independent.

Repeated estimation of the same event requires a resampling experiment. Bias
requires a specified truth reference, such as known-input simulations.

## Reproduction

```bash
python3 tests/summarize_matrix_estimates.py --help
```

The tool consumes an external expanded matrix-audit JSON. Save distributions,
paired shifts, row identifiers and hashes outside the public source checkout.
Detailed catalog distributions are maintained in the analysis workspace.
See [matrix selection](matrix-ratio-threshold-2005_en.md) for selection conventions.
