# Pre-Bootstrap repeatability investigation (2026-09-17)

Follow-up: [FP32 Newton fitting has been restored and verified](newton-fp32-rollback-20260917_en.md).
The unchanged-source statements below refer to the investigation before that restoration.

The R event missing from one March 15–17 production run can be rejected without
calling Bootstrap. A diagnostic replay identifies FFT-plan-dependent roundoff
followed by a sensitive Newton stopping decision. This is not Bootstrap sampling
variability. The production sources, installed executable and settings were not
changed during this investigation.

## Reproduction

The target starts at `20050317T115128`, has initial px/py = −0.01/0.09 s/km,
and max/MAD ≈ 24.2705. Diagnostics retain the original three-day waveform/QC
sequence, frequency band and station geometry. They skip slant-stack calculation,
replay the already observed initial seed, and run only this event's CPU grid and
Newton fits. Bootstrap is never called. All 520 accepted windows match the
production run; the target uses 47 windows and a 698-entry station vector,
including entries without complete three-component waveforms (685 complete).

An excluded single-day pilot demonstrated why retaining the QC history matters:
`eval_deri` carries state between windows, and skipping the preceding windows
changed the target to 46 windows. That pilot is not evidence for the target event.

| FFT planning | Fresh processes | Newton return values |
|---|---:|---|
| Production `FFTW_MEASURE` | 5 | 4, 4, −1, 4, −1 |
| `FFTW_ESTIMATE \| FFTW_UNALIGNED` | 2 | 4, 4 |

Each process repeats the identical in-memory Newton input three times and gets
the same return value every time. A return of 4 is accepted; −1 is rejected.
Every replay reports zero Bootstrap calls. Geometry, weights and initial fitting
parameters match byte-for-byte across all seven processes. The two fixed-plan
processes also have byte-identical spectra.

The logged FFT plans differ. For example, the accepted reference uses an
`hf2_32` stage, while both rejected runs use `hf_32`; those two rejected runs have
identical saved spectra. The rejected/reference R-spectrum difference is
`||R_rejected − R_reference||₂ / ||R_reference||₂ = 2.6343e−16` at Newton entry.
This compares numerical spectra, not event-to-event scatter or estimation bias.
Runtime plan selection and small run-to-run numerical changes with `FFTW_MEASURE`
are documented in [FFTW FAQ 3.8](https://fftw.org/faq/section3.html).

## Where acceptance changes

At iteration 4, both paths are near the same local maximum, with four negative
Hessian eigenvalues. The accepted run finds a trial with relative gain
`+1.65954e−15`, passes the strict `S1 > S0` test, and satisfies normal convergence.
In the rejected run, every trial is no better than the current value:

| Rejected run at iteration 4 | Value |
|---|---:|
| Predicted relative Newton gain | 2.67668e−16 |
| Best observed relative trial gain | −3.68787e−16 |
| Roundoff power tolerance (`8 * epsilon`) | 1.77636e−15 |
| Full normalized Newton-step infinity norm | 1.64979e−8 |
| Allowed step (`sqrt(epsilon)`) | 1.49012e−8 |

The roundoff fallback passes the curvature and power checks, but its full-step
bound is exceeded by about 11%. The function therefore returns −1 despite an
improvement prediction already at floating-point resolution. This explains how
a minute FFT difference becomes a catalog membership difference.

The preceding iteration 3 has an actual relative power increase of
`8.816484504155519e−9`, which fails the normal stopping condition
`abs(epsilon) < 1e−9` and therefore leads to iteration 4. The distance-angle
quantity used in the stopping test, `dprm(2)`, is −0.007693775 degrees, so the
other existing condition passes. The gain
at iteration 3 measures the change across that update, not the error remaining
afterward. That update reaches the neighborhood of the maximum, where the next
predicted gain falls to `2.67668e−16`. The failure to establish strict improvement
then returns a line-search failure before normal convergence can be checked.
The previous FP32 Newton function compared below also converges at iteration 4,
not iteration 3, on this input.

A separate executable reads the saved spectra, weights, geometry and parameters,
without waveform I/O, FFT planning, grid search or Bootstrap. It reproduces the
accepted and rejected returns. Thus identical numerical inputs give identical
results; the FFT output difference and stopping sensitivity are both observable.
The original production run did not save FFT plans, so its exact plan cannot be
recovered. This targeted replay did not trace the other candidates; the
all-candidate scan below extends that investigation.
On failure, the returned parameter object still holds the caller's grid input;
the near-maximum point is recorded in `#DIAG_ITER iteration=4`, not in the failed
`#DIAG_RESULT` parameter fields.

## Comparison with the previous Newton function

Both `FFTW_MEASURE` and strict `S1 > S0` are present in the initial public commit
`e551a69`. The FP64 4x4 Newton operations and the additional stopping condition
using `sqrt(eps_double)` were introduced in `dc38340`. The current uncommitted
`selected` changes do not modify this function.

The same two saved inputs were replayed with the Newton function extracted from
the preceding revision `485d761`, that function with only `Matrix4f` / `Vector4f`
converted to double, and the current function. The previous function is otherwise
unchanged except for its name. Objective evaluation, Hessian assembly, initial
parameters and build settings are shared across the three variants. No waveform
reading, FFT or Bootstrap is performed. This isolates the Newton change; it is
not a complete historical pipeline replay.

| Newton implementation | Saved input `measure-12-1` | Saved input `measure-12-3` |
|---|---|---|
| Previous: FP32 4x4 operations, no fallback | Accepted (4 iterations) | Accepted (4 iterations) |
| FP64 4x4 operations only | Accepted (4 iterations) | Rejected (−1) |
| Current: FP64 plus fallback | Accepted (4 iterations) | Rejected (−1) |

For this event, switching to double changes the iteration path and final rounding,
and acceptance changes in combination with the existing strict improvement test.
The fallback does not introduce a new rejection; it fails to recover this
post-conversion failure. This should therefore be treated as an acceptance
regression associated with the change. The longstanding source of FFT roundoff
does not establish that this event was already unstable with the previous solver.

The original three-day validation retained all 75 selected events with fixed FFT
planning, but did not verify repeatability across multiple FFT plans. This frozen
comparison covers two inputs, not repeatability of all events in the older version.
`run_history.py`, `history-replay.json` and logs are saved in the Git-ignored
investigation directory below. Production sources and executable remain unchanged.

## Additional scan of all candidate stopping decisions

To investigate whether the issue is confined to the target, CPU slant stack,
grid search and Newton fitting were run for all 147 March 15–17 candidates,
recording every step, eigenvalue, predicted gain and trial power. Bootstrap and
catalog emission were disabled; Newton decisions were unchanged. One fresh
process used `ESTIMATE | UNALIGNED` and another used production `MEASURE`.
Candidate coordinates/order and all 520 accepted windows matched. All 147
fixed-plan return values also matched the previous full audit with Bootstrap.
Only rejected candidates were additionally replayed with the previous Newton
function on the same numerical input.

| Stopping outcome | Fixed plan | MEASURE |
|---|---:|---:|
| Normal convergence | 144 | 140 |
| Roundoff fallback | 1 | 1 |
| Iteration limit | 2 | 3 |
| Line-search failure | 0 | 3 |
| Accepted candidates passing max/MAD selection | 75 | 74 |

The three MEASURE line-search failures have different numerical characteristics:

- The known March 17 11:51:28 R candidate (max/MAD ≈ 24.27) has predicted gain
  `2.68e−16`, step norm `1.65e−8`, and four negative eigenvalues. The previous
  Newton function accepts it in four iterations.
- A March 16 05:55:44 T candidate (max/MAD ≈ −2.19) also stops near numerical
  resolution: predicted gain `2.23e−15`, step norm `1.54e−7`, and four negative
  eigenvalues. The best trial has relative change `−9.84e−16`. Its step and
  predicted gain exceed the fallback bounds; the previous Newton function
  accepts it in nine iterations. It is below the ordinary max/MAD cutoff.
- Another March 16 11:51:28 T candidate (max/MAD ≈ −2.71) has step norm about
  0.58, one positive eigenvalue, and a best trial that decreases power by about
  0.41%. The evidence does not support rescuing it as the same roundoff case.
  The three iteration-limit failures also retain appreciable steps/gains;
  failures cannot simply be relabeled as convergence.

Of the 75 selected events in the fixed-plan run, one uses the March 16 roundoff
fallback. Among the other 74 normally converged events, 11 have a final positive
gain at most `8 * eps_double`. This indicates possible sensitivity; it does not
establish that all 11 are misclassified or will disappear. Curvature, starting
points and iteration paths differ across peaks: some stop at the normal
tolerance, some obtain a positive gain at numerical resolution, some use the
fallback, and some are rejected. Reliable catalogs require distinguishing
adequately converged points from failures with unresolved updates. One run per
planning policy over three days does not estimate loss rates in other periods.

Build/run/summary scripts, full logs and candidate records are local under
`docs/benchmarks/pre-bootstrap-stop-scan-20260917/`; generated sources are under
`build-clang/pre-bootstrap-stop-scan/`. Both locations are Git-ignored.

## Follow-up and evidence

### Recommended rollout policy at the time of this investigation

Follow-up: the requested Newton-only rollback has now been validated and applied.
See [FP32 rollback verification](newton-fp32-rollback-20260917_en.md).
The proposal and unchanged-source statements in this investigation describe
the state before that follow-up.

To preserve the established analysis while fixing demonstrated defects, the
interim production candidate should restore only the previous Newton 4x4
arithmetic and stopping behavior. Retain the FP64 covariance inverse: its
isolated three-day comparison removed NaN while preserving acceptance and all
non-covariance output columns
([covariance validation](parameter-covariance-precision-20260917_en.md)).
Self-term subtraction, peak-search fixes and parallelization are separate
changes; do not revert the whole commit. Recheck the restored Newton variant
over three days and multiple FFT plans before adoption. Agreement on the two
frozen inputs does not establish repeatability of the entire previous version.

Revisit the FP64 Newton stopping decision on a separate branch. Where the
current gradient, curvature and full step establish adequate convergence,
consider allowing successful termination without requiring an additional
positive power increment. Separate gradient and parameter-step tolerances also
appear in, for example,
[SciPy BFGS gtol/xrtol](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-bfgs.html),
but its defaults should not be transferred to this problem. Select tolerances
from parameter normalization and required accuracy, retaining the curvature and
power checks described below.

Before adopting the FP64 variant, check the known accepted/rejected inputs,
saddle/flat-region rejection tests, March 15–17 under multiple FFT plans, and a
separate data period. Explain event-level acceptance changes and compare
location, slowness, energy ratios and common events in `all` / `selected`.
Making the missing event disappear under fixed FFT planning, or matching only
row counts, is insufficient. At the end of this investigation no rollback,
branch creation or production configuration change had been performed; the
linked follow-up subsequently applies only the Newton rollback.

### Candidate changes to the FP64 version

The stopping rule should handle improvements at numerical resolution consistently
while retaining rejection of saddles and unresolved large parameter steps. The
saved accepted/rejected inputs provide a regression case. Fixed FFT plans can
improve repeatability; this investigation has not changed the production policy
or selected a new convergence threshold.

### Proposed improvements (not implemented)

1. **Revise the stopping decision first.** Require improvement for an actual
   update, but permit stopping at the current point at numerical resolution.
   Alongside a negative-definite Hessian, finite quantities and negligible
   predicted gain, check the full Newton step before backtracking against
   tolerances for each parameter. `sqrt(epsilon)` alone does not express the
   required parameter accuracy; inspect normalization and physical units before
   selecting tolerances. Do not tune a bound to this single event or merely
   change `S1 > S0` to `S1 >= S0`. Preserve rejection of saddles, singular
   curvature, unresolved meaningful steps and meaningful power decreases.
2. **Make FFT plans reusable.** For a fixed build and runtime environment,
   candidates are `FFTW_ESTIMATE` for verification or saved measured wisdom.
   Evaluate wisdom reuse for performance, recording FFT size/type, array
   alignment, FFTW version and effective plan. The fixed-plan tests here used
   `ESTIMATE | UNALIGNED`; they do not establish the performance or exact
   repeatability of `ESTIMATE` alone. `UNALIGNED` restricts SIMD, so measure
   before adopting it as a default
   ([FFTW planner flags](https://fftw.org/doc/Planner-Flags.html)). Fixed plans
   improve repeatability; the stopping rule still needs verification against
   rounding differences in other environments.
3. **Verify and record acceptance changes.** Replay the saved accepted/rejected
   inputs and retain the existing rejection tests for saddles, flat regions and
   large steps, then rerun March 15–17 with multiple plans. Keep peak detection
   and max/MAD thresholds unchanged; compare changed events, location, slowness,
   energy ratios and common events between `all` and `selected`. Fix Bootstrap
   randomness separately when comparing uncertainty estimates. Distinguish
   normal convergence, roundoff, line-search failure, unsuitable curvature and
   iteration-limit stops to make numerical rejections traceable.

Build scripts, seven run logs, plan descriptions, input differences and frozen
replays are local under `docs/benchmarks/pre-bootstrap-repeatability-20260917/`.
Generated sources and snapshots are under `build-clang/pre-bootstrap-repeatability/`.
Both locations are Git-ignored. `build_probe.py` builds the diagnostic,
`run_probe.py MODE REPEAT [THREADS]` runs a fresh process, and `run_frozen.py`
reproduces the recorded accepted/rejected Newton inputs.
