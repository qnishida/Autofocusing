# Newton termination at objective roundoff (2026-09-17)

Historical validation: subsequent FFT-plan checks exposed a remaining acceptance
regression. The implementation and predicate-only test described here were
replaced by the [interim FP32 restoration](newton-fp32-rollback-20260917_en.md).
Commands and branch status below refer to the earlier implementation.

This follows the [Newton precision comparison](newton-precision-20260917_en.md).
Double Newton arithmetic exposed a pre-existing termination failure: a selected
R candidate reached the previous solution, but all trial powers were no larger
than the current power in floating-point arithmetic. Strict `S1 > S0` rejected
every trial before the ordinary convergence condition could be evaluated.

## Restricted fallback

Strict improvement and the existing convergence predicate remain the ordinary
path. Only after all backtracking trials fail, the current point may be retained
as converged if all of the following hold:

1. Current power is finite and positive; every trial power is finite.
2. The current four Hessian eigenvalues are finite and strictly negative.
3. The full Newton step, expressed in the existing characteristic parameter
   scales W, has infinity norm no greater than `sqrt(eps_double)` ≈ 1.49e-8.
   This is the full proposed step, not a small step produced by backtracking.
4. Its quadratic predicted gain `-0.5 * g_scaled.dot(step_scaled)` is finite,
   nonnegative, and at most `8 * eps_double` times the current power.
5. The best trial power differs from the current power by at most that same
   relative bound, approximately 1.78e-15.

The quadratic power change is second order in the Newton step, motivating
a square-root-epsilon step scale together with an epsilon-scale gain check.
The factor eight is a small numerical allowance, not a rigorous error bound
for the entire waveform calculation and not a threshold calibrated from event
max/MAD. It remains subject to validation beyond this sample.

A tie alone is insufficient: a broad flat region can have no measurable power
gain while still having a substantial Newton step. A small step alone is also
insufficient: the curvature, predicted gain and actual trial powers are checked.
No worsening trial is accepted. The current parameters and power are retained;
the recorded relative power change epsilon is zero and the stop is logged as
`#NewtonStop reason=roundoff`. Bootstrap scaling, selection thresholds, spectral
matrices, iteration limit and the existing ordinary signed distance-step
convergence predicate are unchanged.

Step-size convergence tests also appear in established solvers, for example
[PETSc's convergence reasons](https://petsc.gitlab.io/petsc/release/manualpages/SNES/SNESConvergedReason/)
and [SciPy BFGS xrtol](https://docs.scipy.org/doc/scipy-1.15.0/reference/optimize.minimize-bfgs.html).
Those references support the general use of step tests, not this application's
specific thresholds or statistical-error interpretation.

## Checks and targeted replay

`tests/newton_roundoff.cpp` exercises the production predicate in 24 cases:
a concave quadratic whose improvement rounds away, exact stationarity, tiny
negative rounding differences, the observed update/gain bounds, power-unit
rescaling, large steps on flat regions, resolvable predicted gains, genuine
power decreases, saddles, singular curvature, negative gain, nonpositive power,
and nonfinite inputs. All cases pass.

A two-day replay preserves the March 15→16 QC history and fits only the lost R
seed, (-0.01,0.09) s/km, in the window starting 2005-03-16 11:51:28. It retains
the exact current point from the prior failure diagnostic and now converges
at iteration 8. Current/best trial power is 1.1304596832815571e-6, full scaled
step infinity norm is 1.0195759592448435e-8, and predicted relative gain is
5.5164869323820266e-16. No trial point is substituted for the current solution.

## Full three-day validation

All 13 selected CTests pass. The March 15–17 CPU run and both comparisons pass
the recorded checks, with matching waveform/catalog hashes, QC windows,
original seeds and initial Newton parameters.

| Measure | Float Newton | Double, strict line search | Double with fallback |
|---|---:|---:|---:|
| Initial candidates | 147 | 147 | 147 |
| Emitted rows | 143 | 144 | 145 |
| Iteration-limit failures | 3 | 2 | 2 |
| Line-search failures | 1 | 1 | 0 |
| Selected R/T/U | 31/10/34 | 30/10/34 | 31/10/34 |

The fallback fires exactly once, for the diagnosed R candidate. All other
optimizer results and all accepted-step records are identical to the strict
double run. Its 144 previously emitted rows are unchanged in every catalog
column. The recovered row matches the targeted replay exactly. All 145 output
rows are finite, and the two remaining iteration-limit failures are not converted
to successes.

Compared with float Newton, all original 143 emitted candidates and all
75 selected events remain present. The two additions from double arithmetic
are both below the existing selection thresholds. Among the 75 selected rows,
slowness vectors and R/U ratios are identical at catalog precision; maximum
geographic-location difference is 0.01112 km. T/R and T/U have maximum relative
change 3.861e-6 (0.0003861%) for the 30 rows with positive terms in both runs;
the other 45 rows are excluded from those positive-ratio comparisons.

This resolves the demonstrated roundoff-induced loss within the tested CPU
configuration. It does not establish equivalence over all dates, frequency
bands, slowness grids or GPU inputs. As already recorded in the type-only
study, some weak unselected candidates converge to different solutions in
float versus double. The two known unselected rows with negative covariance
diagonals also remain; this fallback is not a correction for their post-fit
curvature or statistical uncertainty.

## Reproduction and status

```bash
cmake -S . -B build-clang -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-clang --target test_newton_roundoff --parallel 4
ctest --test-dir build-clang -R '^newton_roundoff_convergence$' --output-on-failure
```

The local probe builder and detailed records are retained under
`docs/benchmarks/newton-roundoff-20260917/`, ignored by Git. They require the
saved deterministic Newton-double audit setup. Before building, the builder
reverses the fallback changes and verifies the exact prior production-source
SHA-256; its full probe changes no other numerical behavior. Both full runs use
fixed bootstrap seeds, FFTW ESTIMATE | UNALIGNED, three components, 0.1–0.25 Hz,
±0.165 s/km, dp=0.005 s/km and eight CPU threads.

The implementation, regression tests and validation notes are included on
`test/peak-search-boundaries`. Main, remote and installed binaries are unchanged.
