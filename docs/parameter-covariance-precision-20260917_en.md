# Parameter covariance precision (2026-09-17)

This follows the [peak-search validation](peak-search-validation-20260917_en.md).
The changes are included on `test/peak-search-boundaries`, based on `485d761`. Main, the remote and the installed binary have not been updated.

This records the covariance-only precision change. The subsequent
[Newton precision change](newton-precision-20260917_en.md) is evaluated separately;
references below to float Newton arithmetic describe this earlier baseline.

## Cause and correction

One newly admitted T candidate produced NaNs in all ten parameter-covariance
columns. Its window starts at 2005-03-16 11:51:28 and its seed is
(0.063, -0.054) s/km. A deterministic replay reproduced the original row.
All entries of H/sigma were finite, but its determinant was approximately
5.205e39, beyond float's maximum of approximately 3.403e38. The single-precision
4×4 inverse failed; the double-precision inverse was finite.

`est_dist_boot()` now delegates covariance conversion to
`set_parameter_covariance()`, which forms and inverts H/sigma in double
precision. The output remains

```text
C_output(i,j) = -inverse(H/sigma)(i,j) * W(i) * W(j)
W = (0.06, pi/2, pi/2, 0.04/(30*111))
```

Even the existing float representations of W are retained. This preserves
the convention that downstream location-error propagation cancels the W
scales in its Jacobian. There is no 0.06 rescaling of published location errors.
Bootstrap sampling, sigma, the optimizer, convergence tests, peak selection
and spectral-matrix calculations are unchanged by this precision correction.

## Regression tests

`tests/parameter_covariance.cpp` calls the production covariance conversion.
It covers a coupled ordinary Hessian at two sigma values, the observed H/sigma,
and the reconstructed raw Hessian with its actual bootstrap sigma. The fixture
also verifies that the old float inverse is nonfinite despite finite inputs.
An independent, equilibrated long-double LU solve checks every output entry,
including cross terms and normalization. On this Apple arm64 host, long double
and double both have 53 mantissa bits; independence comes from equilibration
and the different solution method, not additional mantissa precision.
Maximum relative disagreement is
3.81e-12 for the observed fixture and below 7e-16 for the ordinary Hessian.
All 12 selected CTests pass after rebuilding the affected targets.

## Real-data validation

The March 15–17, 2005 CPU run uses 0.1–0.25 Hz, ±0.165 s/km, dp=0.005 s/km
and eight threads, with the same fixed bootstrap seeds, FFTW settings and QC
history as the peak-only run. All input and catalog hashes match.

- The same 147 initial candidates produce the same 143 events. All optimizer
  records, iteration counts and outcomes are identical.
- All 143 rows are now finite; all original audit checks pass without relaxation.
- All 28 columns outside covariance (21–30) match exactly as written. This
  includes locations, slowness, max/MAD, bootstrap power/sigma, SRR/STT/SUU and
  the cross terms. Energy ratios and downstream selection are unchanged.
- The same 75 rows pass selection: R=31, T=10, U=34.
- Among previously finite covariance entries, relative changes
  `abs(new-old)/abs(old)` have median 0, 95th percentile 5.90e-6 and maximum
  5.55e-5. For selected rows, the maximum is 3.80e-5 (0.00380%). These statistics
  use the catalog's six-significant-digit text, not full internal precision.
  No previously finite covariance entry changes sign.

The recovered T row has negative covariance diagonals and fails max/MAD > 7
(ratio approximately -2.709). It therefore remains unsuitable for uncertainty
interpretation despite being finite. Another unselected R row already had a
negative diagonal before this change and retains it. None of the selected
75 rows has a negative diagonal; this alone does not validate uncertainty
coverage or positive definiteness.

## Scope of the result

This corrects the observed numerical overflow, not the statistical model for
uncertainty. Double precision does not guarantee a usable covariance for every
singular or poorly constrained Hessian, or for zero bootstrap sigma. No
regularization, event rejection or replacement of nonfinite values is added.

## Reproduction

```bash
cmake -S . -B build-clang -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-clang --target test_parameter_covariance --parallel 4
ctest --test-dir build-clang -R '^parameter_covariance_range$' --output-on-failure
```

`tests/build_covariance_probe.py` takes the saved deterministic peak-search
probe and replaces its bootstrap function with the current production function
and covariance helper, preserving the fixed random seed. It checks that the
old function matches the baseline and that peak search is unchanged.
`tests/compare_covariance_runs.py` checks input hashes, QC windows, candidates,
optimizer records, selection, and every non-covariance catalog column; it
requires the corrected run to pass the existing audit checks.

Local commands after building the saved peak-search audit setup:

```bash
python3 -B tests/build_covariance_probe.py build-clang
python3 -B tests/run_matrix_audit_events.py \
  build-clang/covariance-double-probe/cal_ccf_covariance \
  /Volumes/Seismic_Data/hdf5/Hi-net ../moment_loc_76_24 \
  build-clang/covariance-double-mar15-17 --start 2005-03-15 --days 3 --threads 8
python3 -B tests/compare_covariance_runs.py \
  build-clang/peak-search-mar15-17-fixed build-clang/covariance-double-mar15-17 \
  build-clang/covariance-double-comparison.json
```

Raw catalogs, logs, generated probe source, build commands and hashes are
retained locally under `docs/benchmarks/parameter-covariance-precision-20260917/`,
which is ignored by Git. The original failed float run remains in the separate
peak-search records.

## Clarification: matrices, signs and CPU precision

The 3×3 spectral matrix S_RTU and the 4×4 parameter covariance are different
outputs. Removing station self terms from S_RTU does not preserve positive
semidefiniteness; a negative corrected power estimate can occur. The parameter
covariance instead comes from -sigma W H^-1 W, not inversion of S_RTU. For the
fixed spectra and weights used in one Hessian evaluation, the scalar self term
is independent of the four fitted phase parameters. Its second derivative is
zero. Subtraction can affect the bootstrap sigma, but a positive scalar sigma
cannot reverse covariance definiteness.

For the formerly NaN T example, a 60-digit Decimal LDL decomposition of the
recorded H/sigma gives pivots approximately (-2.01770756e8, +4.83715177e6,
-6.83665288e11, +7.80057330e12): two negative and two positive eigenvalues by
congruence. Thus the Hessian used for covariance is indefinite even before
inversion in float; recovering finite numbers does not make it a covariance
of a valid local maximum. Why this candidate passed the optimizer's earlier
curvature check is not yet isolated. The code checks curvature before the last
accepted step and updates the R/T rotation after optimization, before evaluating
the bootstrap Hessian; these are separate evaluation points/inputs to investigate.

CPU processing is mixed precision. Objective and Hessian accumulation use
double, while the Newton eigensolver and update vectors still use float.
The covariance inverse ran on the CPU in float before this fix and now uses
double on that same CPU path. The demonstrated NaN was from a CPU-only run;
this change is not a CPU/GPU comparison and does not qualify the remaining
float Newton step against an all-double implementation.
