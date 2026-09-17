# Newton iteration precision (2026-09-17)

This records the type-only experiment. The subsequent
[roundoff termination study](newton-roundoff-20260917_en.md) addresses the
selected-event loss identified here; the original comparison is retained.

The initial public release, `e551a69`, already uses `Matrix4f` and `Vector4f`
for Newton iteration, including its self-adjoint eigensolver. Its parameter
covariance inverse also uses float. These were original CPU computations,
not a later conversion introduced for GPU support.

On `test/peak-search-boundaries`, the Newton Hessian workspace, eigensolver,
eigenvalues/eigenvectors and update vectors in `src/cal_ccf.cpp` now use double.
The Newton scaling vector W also uses double with the same constants/formula.
The objective, gradient and Hessian accumulation already used double.
This follows the separate
[covariance inverse correction](parameter-covariance-precision-20260917_en.md).
The output-covariance normalization retains its existing float W values.

Line search, the convergence predicate, iteration limit, peak-selection rules,
GPU kernels and output format are unchanged. The Newton routine executes on
the CPU for both CPU and GPU-assisted runs. Validation here is CPU-only;
it is not a new GPU equivalence qualification or a controlled timing benchmark.

## Validation setup

Compare against the saved covariance-double/Newton-float run, with the same
147 original candidate seeds from March 15–17, 2005: three components,
0.1–0.25 Hz, ±0.165 s/km, dp=0.005 s/km, eight threads. Bootstrap seeds are
fixed; both runs use FFTW ESTIMATE | UNALIGNED and the same QC history.
The probe builder reverses only the Newton type replacements and verifies
that the resulting production-source SHA-256 equals the saved baseline hash.
It then applies those replacements to the already instrumented baseline.

All 12 selected CTests pass after rebuilding the affected targets. The
real-data comparison records candidates and failed fits as well as emitted
rows, and distinguishes selected events from weak candidates. Relative
catalog differences use the six-significant-digit text output. Energy ratios
are compared only when numerator and denominator are positive in both runs,
with exclusions counted explicitly.

## Real-data result: not yet suitable for unchanged production selection

| Measure | Newton float | Newton double |
|---|---:|---:|
| Initial candidates | 147 | 147 |
| Converged/emitted | 143 | 144 |
| Iteration-limit failures | 3 | 2 |
| Line-search failures | 1 | 1 |
| Selected R/T/U | 31/10/34 | 30/10/34 |

Both runs pass the existing audit checks, including all-finite output. Input
hashes, QC windows, original seeds and the parameters entering Newton iteration
are identical. Outcomes or iteration counts change for 15 candidates.
Two previously failed candidates converge, but both fail the existing event
selection thresholds. One previously selected R candidate instead fails line
search: window 2005-03-16 11:51:28–17:56:32, initial (px,py)=(-0.01,0.09) s/km,
max/MAD approximately 25.40. Consequently, passing the audit is not sufficient
to approve this type-only change for production.

For the 74 selected events emitted in both runs, slowness vectors and SRR/SUU
are identical at catalog precision. The largest geographic-location difference
is approximately 0.01112 km. R/U is identical in all 74; T/R and T/U differ by
at most 3.861e-6 relative (0.0003861%) among the 29 rows with positive terms in
both runs. The other 45 rows are excluded from those positive-ratio comparisons.

Weak, unselected common candidates can change substantially: five have
slowness-vector differences above 1e-6 s/km, with maximum 0.02697 s/km and
maximum geographic-location change approximately 11005 km. Precision changes
can lead these candidates to different solutions; the small differences among
selected common events do not characterize all candidates. The two previously
identified unselected rows with negative covariance diagonals retain that issue.

The type-only double Newton experiment was not adopted alone; the linked
roundoff-termination follow-up records the additional correction and validation.
Do not equate higher arithmetic precision with unchanged event acceptance.
The strict line-search/near-convergence behavior needs a separate qualification
before production adoption; the acceptance rule has not been loosened here.

### Targeted line-search replay

A two-day replay preserving the March 15→16 QC history, while fitting only the
lost R seed, reproduces all seven accepted steps and the iteration-8 failure.
The current Hessian has four negative eigenvalues. S0 and the best trial power
both equal 1.1304596832815571e-6 in double arithmetic. The quadratic prediction
of relative gain is 5.5165e-16, the W-scaled update norm is 1.0541e-8, and the
distance update is -9.1762e-7 degrees. The current slowness magnitude differs
from the float run's accepted solution by only 4.36e-11 s/km.

The optimizer tests strict `S1 > S0` before applying its usual convergence
condition. When every trial has no representable positive gain, it returns
failure even this close to the prior solution. This is an existing termination
fragility exposed by different arithmetic, not evidence of a grossly displaced
double solution. A small-step or roundoff-aware convergence policy would need
its own criteria and tests; none is implicitly enabled by this experiment.

## Local reproduction and records

The existing deterministic audit setup is required:

```bash
python3 -B docs/benchmarks/newton-precision-20260917/build_probe.py
python3 -B tests/run_matrix_audit_events.py \
  build-clang/newton-double-probe/cal_ccf_newton \
  /Volumes/Seismic_Data/hdf5/Hi-net ../moment_loc_76_24 \
  build-clang/newton-double-mar15-17 --start 2005-03-15 --days 3 --threads 8
python3 -B docs/benchmarks/newton-precision-20260917/compare_runs.py
```

Raw logs, catalogs, generated source, comparison details and hashes are kept
locally in `docs/benchmarks/newton-precision-20260917/`, ignored by Git. The combined implementation and validation notes are included on the test
branch; main, remote and installed binaries are not updated.
