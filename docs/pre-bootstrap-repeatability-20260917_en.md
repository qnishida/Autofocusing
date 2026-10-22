# Pre-Bootstrap repeatability and Newton stopping

This historical investigation concerns the FP64 Newton experiment. The adopted
implementation subsequently [restored FP32 Newton iteration](newton-fp32-rollback-20260917_en.md)
while retaining the FP64 covariance inverse.

## Reproduction

Separate Bootstrap randomness from earlier numerical differences. Freeze the
same spectra, station geometry, weights and initial parameters, and replay the
CPU grid/Newton path without Bootstrap. Preserve preceding waveform/QC history.
Within a process, repeated identical numerical inputs should be checked separately
from FFT-planning differences between fresh processes.

A targeted historical replay returned the same Newton result on repeated identical
in-memory inputs, while fresh `FFTW_MEASURE` processes could accept or reject a
near-converged fit. Fixed `FFTW_ESTIMATE | FFTW_UNALIGNED` controls reproduced the
accepted path in the tested cases. This is not a guarantee for arbitrary inputs.

## Where acceptance changes

Minute FFT differences can reach a strict `S1 > S0` line search when the expected
gain is already at floating-point resolution. If no trial is representably
better, rejection can occur before the normal convergence predicate is tested.
That is distinct from Bootstrap sampling variation or a large parameter shift.
The [FFTW FAQ](https://fftw.org/faq/section3.html) discusses planning-related
run-to-run numerical differences.

The FP64/roundoff experiments did not establish unchanged acceptance for all
candidates. Weak unselected events can reach different solutions even when
strong common events agree. The
[precision experiment](newton-precision-20260917_en.md) and
[roundoff study](newton-roundoff-20260917_en.md) describe these limits.

## Validation guidance

Record accepted and rejected candidates, iterations, line-search gains and
curvature, rather than comparing only emitted catalog rows. Use independent
input periods before adopting a stopping-policy change. Do not loosen production
tolerances merely to make a diagnostic replay pass.

Current synthetic Newton recovery is checked by `newton_wavefront_recovery`.
Frozen spectra, per-candidate traces and machine-specific replay scripts belong
in the analysis workspace. They are not public-source-only reproduction inputs.
