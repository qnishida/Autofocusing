# Optional early event selection (2026-09-17)

Status: the measurements below use the earlier FP64 Newton implementation.
The subsequent [FP32 Newton restoration](newton-fp32-rollback-20260917_en.md)
retains this feature and rechecks `all` / `selected` equivalence. Its CPU test
outputs 143 full-catalog rows and the same 75 selected rows. The timing results
below predate that restoration and have not been remeasured with it.

The new `selected` mode skips per-candidate parameter searches, Newton fitting,
Bootstrap and spectral-matrix evaluation when the initial-grid max/MAD fails
an experiment's cutoff. The default `all` mode retains full-catalog behavior.
Configuration and usage are described in
[Scripts/README.md](../Scripts/README.md#compute-only-events-above-maxmad-thresholds).

Peak discovery, search bounds, local contrast, candidate limits and horizontal
duplicate exclusions are unchanged. All U seeds remain available for horizontal
seed generation even when their own U fits are skipped. The cutoffs are strictly
greater than 7 for R/T and 35 for U by default, configurable independently.
The predicate uses full-precision initial-grid values; printed values can round
to the cutoff. This mode applies the stated max/MAD rule, not a new assessment
of physical source validity.

The executable reports its mode/thresholds through `--event-selection-info` and
the startup log. The launcher verifies support before running `selected` mode,
records the effective settings in the manifest, and rejects unsupported binaries
before creating a run. Per-window/component counters record candidates, skipped
fits, attempted fits and emitted rows. Full-catalog output remains available for
later re-selection, and the output columns and numerical estimators are unchanged.

## Correctness checks

All 17 registered CTests pass, including new strict-threshold/finite-input tests,
real executable configuration and manifest checks, rejection of old binaries,
and the existing peak, covariance, fitting, parallelism and workspace checks.
An existing test's `/usr/bin/time -l` needs OS access unavailable in the sandbox;
its targeted rerun outside the sandbox passed without changing the test.

Real-data comparisons use March 15–17, 2005, CPU, three components, 0.1–0.25 Hz,
±0.165 s/km, dp=0.005 s/km, eight threads. A generated test binary fixes Bootstrap
seeds and FFTW planning; production seed and planning policies are unchanged.

| Measure | all | selected |
|---|---:|---:|
| Initial candidates | 147 | 147 |
| Candidates skipped before fitting | 0 | 72 |
| Candidates fitted | 147 | 75 |
| Bootstrap calls / emitted rows | 145 | 75 |
| Selected R / T / U | 31 / 10 / 34 | 31 / 10 / 34 |

- Default output reproduces all 145 archived rows byte-for-byte.
- Selected output is byte-identical to the 75 selected rows of that full output,
  including all 38 columns, covariances and energy matrices.
- Traces of all 147 candidate coordinates, order and initial peak/MAD values
  are identical between modes, including U-derived horizontal candidates.
- Accepted waveform windows, waveform hashes and source-catalog hash agree.
- A separate horizontal-only March 15 comparison preserves all seven candidates
  and all seven output rows. All already pass the default thresholds, so this
  case verifies compatibility rather than a speed benefit.

Production Bootstrap remains time-seeded. Faster execution can select different
random samples; the exact comparisons above require the deterministic test
controls. Full-precision cutoff handling is covered by boundary tests, while
the real-data selected subset is sufficiently far from rounding ambiguities.

## Three-day timing of the installed executable

March 15–17, 2005 was measured with the installed production executable,
Apple M4 Max, Metal slant stack, CPU fitting/Bootstrap, and 12 OpenMP threads.
The three-component band/grid and cutoffs are as above; `AUTOFOCUSING_GPU_POWER`
is `off`. Production `FFTW_MEASURE` and time-derived Bootstrap seeds are retained.
One warmup per mode is excluded, followed by three alternating pairs in fresh
processes. Timing includes startup, FFT planning, input/output and computation,
but excludes the workspace launcher. OS cache is retained; these are not cold
disk timings. No builds or other test jobs were run concurrently.

| Production measurement | all | selected |
|---|---:|---:|
| Three-day wall-time median | 97.970 s | 52.265 s |
| Three-day wall-time range | 97.708–98.040 s | 52.206–52.443 s |
| Per-day average from the three-day median | 32.657 s | 17.422 s |
| Bootstrap-stage median, three days | 24.080 s | 12.665 s |
| Output rows in the three repetitions | 144 / 141 / 143 | 75 / 75 / 75 |
| Rows passing max/MAD in those outputs | 75 / 74 / 75 | 75 / 75 / 75 |

This is a **46.65% decrease in wall time (1.874× speedup)**. All runs load all
three files, accept the same 520 waveform windows and produce 147 initial
candidates; `selected` skips 72 fits. Its output contains R/T/U = 31/10/34 rows
(17, 29 and 29 rows on the individual days). Input and executable hashes remain
unchanged. The per-day figure is an average across these three days, not the
measured time for each individual day.

**Production repeatability remains an open issue.** The second measured `all`
run loses one above-cutoff R event that appears in the other `all` runs and all
`selected` runs: interval start `20050317T115128`, max/MAD ≈ 24.2705,
px/py ≈ −0.00964356/0.0880855 s/km. The initial candidate counts remain the same;
the difference occurs in fitting/output acceptance before Bootstrap. It is not
a cutoff-boundary case or an effect explained by Bootstrap randomness. Because
it also occurs between identical `all` configurations, it is not specific to
switching selection mode. A subsequent
[pre-Bootstrap investigation](pre-bootstrap-repeatability-20260917_en.md)
reproduced the target's rejection with zero Bootstrap calls: FFT-plan-dependent
roundoff changes the outcome of the Newton stopping rule. The production stopping
rule has not yet been corrected.

For events present in both members of each pair, printed px/py and all five
spectral-matrix entries match exactly. Latitude differs by at most 0.0001° in
the checked production comparisons; Bootstrap covariance and its mean/standard
deviation vary with the time-derived samples. Consequently, the fixed-control
checks above must not be read as a guarantee of bitwise-identical production
catalogs. No solver, threshold or FFT/RNG policy was changed for this benchmark.

A separate Metal replay of all three days with the deterministic test executable
returns 145 full rows and 75 selected rows. The selected subset matches in **all
38 columns byte-for-byte**; all 147 full-precision initial candidate traces and
520 accepted waveform windows also match. The R event identified above is present
in both modes. This confirms early-selection equivalence for this input under
fixed test controls. The subsequent investigation above isolates a reproducible
FFT/stopping-rule cause for the production repeatability issue.

## One-day timing with fixed test controls

Apple M4 Max, Metal slant stack with CPU fitting/Bootstrap, 12 threads,
March 15, 2005, three components and the same band/grid as above. Each mode
has one excluded warmup followed by three alternating measured pairs in fresh
processes. OS cache is retained; no builds or tests run concurrently.
This measurement uses the deterministic test executable described above, with
fixed Bootstrap seeds and `FFTW_ESTIMATE | FFTW_UNALIGNED`.

| Median for one day | all | selected |
|---|---:|---:|
| Whole run | 32.260 s | 13.252 s |
| Bootstrap stage | 7.853 s | 2.809 s |
| Initial parameter-grid stage | 3.050 s | 1.055 s |
| Emitted rows | 46 | 17 |

Whole-run time decreases by **58.92% (2.434× speedup)**. Measured ranges are
31.831–32.406 s and 13.142–13.400 s. Each `all` run exactly reproduces the archived
46-row Metal catalog, and each `selected` run matches its 17-row selected subset
byte-for-byte. Candidate counts and accepted waveform windows agree in every
pair. The eliminated work includes fitting and matrix calculation as well as
Bootstrap. This is a one-day, three-pair comparison on this host; CUDA was not
measured. The speed benefit depends on how many candidates fail the cutoffs.

## Local evidence

Detailed source/build hashes, generated probe sources, catalogs, seed traces and
timing logs are retained locally under
`docs/benchmarks/early-event-selection-20260917/` and
`build-clang/event-selection-probe/`, both excluded from Git.
The local `build_probe.py` builds from the configured compilation database;
`run_checks.py accuracy` repeats the data comparisons, and `run_checks.py timing`
runs the paired timing protocol. The waveform and source catalog remain external.

The installed-executable three-day timing, catalogs, comparison scripts and the
repeatability exception are retained under
`docs/benchmarks/event-selection-mar15-17-20260917/` (Git-ignored).
Run `run_benchmark.py` there for the paired protocol and `analyze_outputs.py` for
the common-event and missing-event comparisons. Metal execution requires access
to the GPU outside the macOS sandbox.
The fixed-control three-day Metal replay is saved separately under
`docs/benchmarks/event-selection-mar15-17-controlled-20260917/`; its
`run_check.py` verifies the complete selected output and candidate traces.
