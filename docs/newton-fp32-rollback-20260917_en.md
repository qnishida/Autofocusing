# Interim restoration of FP32 Newton fitting (2026-09-17)

After [investigating the FFT-dependent acceptance regression](pre-bootstrap-repeatability-20260917_en.md),
restore only the Newton 4×4 arithmetic and stopping behavior from
`485d761800ad6981d505fea4870cc2e85fa3f7d7`. The production `est_dist_grad()`
function body is identical to that revision. Remove the FP64 roundoff fallback.
This is an interim return to established behavior, not a general claim that
single precision is more accurate or guarantees repeatability.

Keep the FP64 parameter-covariance inverse, spectral self-term corrections,
peak-search bug fixes, bootstrap parallelization and `all` / `selected` support.
Objective and Hessian accumulation remain double precision. The production
FFT policy, bootstrap random seeds and original convergence thresholds are
unchanged. In particular, the original signed distance-step check is retained.

## Verification

Use March 15–17, 2005, three components, CPU with 12 threads, 0.1–0.25 Hz,
±0.165 s/km and dp=0.005 s/km. Preserve the full three-day QC history and verify
the hashes of all three waveform files and the reference catalog. All runs
retain 520 accepted input windows and the same 147 initial candidates.

The validation copies use fixed bootstrap seeds `1837 + 104729 * sample_index`
to separate fitting/FFT changes from bootstrap sampling variability.

| Run | Fitted candidates | Output rows | Selected R / T / U |
|---|---:|---:|---:|
| Fixed FFT, `all` | 147 | 143 | 31 / 10 / 34 |
| Normal measured FFT, `all`, three fresh processes | 147 each | 143 each | 31 / 10 / 34 each |
| Fixed FFT, `selected` | 75 | 75 | 31 / 10 / 34 |

Selection means max/MAD > 7 for R/T and > 35 for U. Three distinct FFT plan
descriptions occur across the four full runs. Every candidate has the same
acceptance and iteration count in all four runs; all four rejections are below
selection thresholds. Compared with the earlier fixed-plan FP64 result, the
two fewer outputs are also below these thresholds.

- The fixed-plan 143-row catalog matches the archived FP32 Newton / FP64
  covariance baseline byte for byte in all 38 columns.
- Across FFT plans, all catalog columns except the final relative power gain
  epsilon are identical. Epsilon differs by at most 1.53e−15 in absolute value.
  Location, slowness, covariance and energy outputs match at catalog precision.
  Internal values need not be bit-identical: the largest selected slowness
  magnitude difference is 5.26e−12 s/km and distance-angle difference is
  1.85e−8 rad (about 0.12 m of arc distance).
- Selected R/U ratios match for all 75 rows. Positive T/R and T/U ratios match
  for the 30 rows with positive terms; component signs and energy fractions
  also match for all selected rows. Negative bias-corrected powers are not
  silently discarded to claim positive-ratio agreement.
- `selected` output is byte-identical to filtering the fixed-plan `all` output,
  including all covariance and bootstrap columns. It skips the other 72 fits.
- The March 17 11:51:28 R seed (-0.01, 0.09) s/km, max/MAD ≈ 24.27, is accepted
  in four iterations in every run. Both saved inputs that previously gave
  opposite FP64 decisions also return four iterations and exactly reproduce
  the previously recorded FP32 parameter values.
- All output values are finite. Selected covariance diagonals are positive.
  The two known unselected rows with negative covariance diagonals remain;
  this rollback does not validate their uncertainty interpretation.

Replace the removed fallback-predicate test with `newton_wavefront_recovery`.
It fits a known synthetic curved wavefront from perturbations on either side,
checks parameter recovery and negative final curvature, and rejects zero power
and an invalid initial distance. Retain the independent covariance-overflow
regression test.

All 17 CTests pass after rebuilding. The timing-tool test requires an
unsandboxed retry because macOS `/usr/bin/time -l` is restricted in the sandbox;
the other 16 pass there directly. This test uses stub drivers and is not a GPU
hardware validation.

## Installed executable with production settings

Rebuild and install `bin/cal_ccf_clang`, then rerun all three days in `selected`
mode using its unmodified measured FFT policy and time-derived bootstrap seeds.
It again emits 75 finite rows, with the same acceptance and iteration counts
and positive covariance diagonals. The originally problematic R event still
converges in four iterations.

The initial strict comparison of non-bootstrap columns does **not** pass:
one different R candidate at March 17 11:51:28, seed (-0.01, 0.045) s/km,
max/MAD ≈ 36.76, has a small distance-update change and a change in the last
printed STT digit. Location and slowness columns still match. The difference
in column 10 corresponds to about 0.30 m of arc distance. STT and positive
T/R, T/U ratios differ by at most 3.861e−6 relative (0.0003861%); R/U is
unchanged. The maximum absolute energy-fraction difference is 3.86e−11.

Retain the interim restoration given unchanged acceptance and these small
physical differences; do not claim bitwise equivalence for arbitrary FFT
plans or builds. The exact FFT plan of this uninstrumented run was not logged,
so the specific cause of this small difference is not isolated here. Bootstrap
covariance, mean and sigma are expected to vary with production seeds and are
checked for finite values rather than equality with fixed-seed output. The
failed strict assertion and quantitative follow-up are both preserved locally;
no production tolerance was relaxed to force agreement.

## Reproduction and scope

```bash
cmake -S . -B build-clang -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-clang --parallel 4
ctest --test-dir build-clang --output-on-failure
```

Detailed build/run/summary scripts, input hashes, catalogs and logs remain local
under `docs/benchmarks/newton-fp32-rollback-20260917/`; generated validation code
and binaries are under `build-clang/newton-fp32-rollback/`. Both are Git-ignored.
The adopted source is checked against the tested proposal, apart from its
documentation comment.

This validates the tested CPU data and settings, not every date, frequency
band, platform or GPU input. Revisit FP64 Newton stopping separately, using
the saved failure cases and an independent data period before adopting it.
