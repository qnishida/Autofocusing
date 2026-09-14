# Ordered CPU fitting parallelism

## Scope and branch workflow

The baseline is CUDA commit `977736f`; `6a7a4bd` adds profiling and an independent
frozen-kernel regression without changing the numerical kernels. Development
uses `gpu` as the combined CPU/CUDA/Metal integration branch. The three changes
are evaluated independently in `test/cpu-hessian-parallel`,
`test/cpu-objective-parallel`, and `test/cpu-rotation-parallel`.

Only passing changes are merged into `gpu`. CPU-only and CUDA checks run on
Linux. Merging the combined GPU history into `main` additionally requires the
planned M4 Max CPU/Metal build, kernel/event regression and performance checks.
The Linux results do not imply that Metal has been tested. No push or binary
installation is part of this work.

## Implementation and invariants

The adopted Hessian implementation distributes independent time-window stacks
over OpenMP workers. At most 16 windows of intermediate results are retained.
Station order and the frequency phase recurrence are unchanged. Final sums
are accumulated serially in the original window/frequency/component order;
there is no OpenMP floating-point reduction. The objective experiment uses the
same ordering approach, but is not included in `gpu` because it regressed the
real initial-grid path.

Single-thread, fewer-than-two-window and already-parallel Hessian callers use the
original serial kernel. The outer initial-grid evaluation therefore does not
create nested worker teams. `cal_S` remains the baseline implementation in `gpu`. Serial and parallel implementations are both
checked against frozen kernels from `977736f`; their matching arithmetic must
be maintained when changing scientific formulas in the future.

Rotation partitions window/station pairs. Each worker owns separate R/T
spectra and weights, retaining the scalar expressions and leaving E/N/U data
unchanged. Existing geographic-coordinate calculation is retained.

The code uses double precision and keeps search parameters, convergence rules,
production random seeds, FFTW plans and GPU kernels unchanged. No fast-math,
native-architecture compiler change, next-day prefetch or GPU objective policy
change is included. The recommended CUDA configuration remains
`AUTOFOCUSING_GPU_POWER=bootstrap`.

## Verification protocol

`test_cpu_parallel` compares objective values (with and without bias), all
gradient/Hessian entries, and every rotated spectrum/weight with an independent
frozen reference. It covers windows 0/1/3/16/17/48, frequencies 1/31/155,
stations 1/2/31/650, zero/nonuniform weights, coherent/random inputs, horizontal
and 3c rotation, and concurrent outer OpenMP callers. Thread counts are
1/4/8/16. The componentwise bound is `1e-12 * abs(reference) + 1e-30`, adopting
the existing CPU regression tolerance for the new derivative checks. Exact
equality is additionally reported, not assumed across compilers or machines.

`AUTOFOCUSING_PROFILE=1` now reports `#CPU_PROFILE` call counts and inclusive
times for objective, Hessian, rotation and fitting. `serial_caller` measures
elapsed call time, including any internal parallel work. `outer_parallel`
is the sum of elapsed worker-call times, **not** elapsed pipeline time.
Fitting includes objective/Hessian calls: do not add these categories together.
The existing segment/load/power profiles remain the source of stage wall time.

The real-event driver uses fixed FFTW flags and Bootstrap seeds only in the
test binary. Paired comparisons use the same backend/power mode before and
after each CPU change. Event fields, including Bootstrap, must match at output
precision; only the existing convergence-residual absolute bound `1e-12`
applies. Accepted windows and initial grid choices must match too. This does
not reuse the relaxed CPU-versus-FP32-GPU bounds for CPU implementation changes.

Performance uses the supplied January 1–3, 2004 horizontal data (20 events),
five alternating fresh-process trials per binary, 16 threads, and one excluded
warmup per binary. Caches are not purged. No agent-started builds, other benchmarks or
regressions run concurrently with adoption measurements. Five-day qualification
(29 events) may omit warmups and its times are not adoption measurements.

Adoption requires at least 10% reduction in the target serial-caller stage and
no more than 5% whole-run regression. One-thread and already-parallel grid
benchmarks must not show a reproducible regression exceeding 5%. The combined
configuration must improve whole-run time. Real 3c input and Metal GPU checks
remain separate qualifications.

## Reproduction

Set `OMP_NUM_THREADS=16` to enable the intended worker count; the new CPU work
uses OpenMP in CPU, CUDA and Metal builds. No new backend setting is required.
Hessian uses at most 16 workers for its window batches. A one-thread invocation
uses the original serial Hessian kernel. Rotation uses the configured team size.


Run from the repository root so the velocity model is found. Use separate build
directories or preserve each deterministic executable before changing branches.
The baseline and changed executables must include the same profiling code.

```bash
cmake --build build-cuda -j 4
OMP_NUM_THREADS=4 ctest --test-dir build-cuda --output-on-failure
build-cuda/test_cpu_parallel --benchmark hessian 16 9 48
build-cuda/test_cpu_parallel --benchmark objective 16 9 48
build-cuda/test_cpu_parallel --benchmark rotation 16 9 48
build-cuda/test_cpu_parallel --benchmark grid 16 9 48
python3 tests/build_io_probe.py build-cuda
python3 tests/run_cpu_parallel_events.py BASELINE_PROBE CURRENT_PROBE \
  "$HINET_ROOT" "$CMT_CATALOG" NEW_OUTPUT_DIRECTORY \
  --backend cuda --power bootstrap --threads 16 --days 3 --repeats 5
```

For correctness qualification use `--days 5 --repeats 1 --warmups 0` and repeat
with `--backend cpu --power off`. CUDA power modes `off`, `grid`, and `all` are
checked as regression paths as well. Use `--target hessian|objective|rotation`
for individual adoption runs; their Boolean adoption results are in the JSON
report and are evaluated separately from correctness failures.

## Results

Linux measurements use a Threadripper 3990X and RTX PRO 2000 Blackwell, GCC
13.3 Release, with native-architecture optimization disabled. These are warmed
NAS measurements, not cold-I/O comparisons with a Mac.

| Independent change | Whole-run median, before → after | Target stage, before → after | Initial-grid change | Decision |
| --- | --- | --- | --- | --- |
| Hessian (`989a4e4`) | 39.343 → 27.817 s | 13.298 → 1.552 s | +3.77% | Adopt |
| Objective (`ff619dd`) | 39.389 → 37.886 s | 2.469 → 0.432 s | +15.83% | Hold |
| Rotation (`d7df22e`) | 39.342 → 38.795 s | 0.830 → 0.375 s | −3.11% | Adopt |

All five pairs in each row produced byte-identical event outputs. The objective
experiment accelerates its target stage but violates the existing-grid 5%
regression gate. A direct-dispatch follow-up (`81426e8`, one diagnostic pair)
still showed +14.55% grid time. Its whole-run timing is not an adoption result.
Both experiments are preserved on `test/cpu-objective-parallel`, without merging
them into `gpu`. The reason for the grid regression remains unresolved.

Nine-pair synthetic Hessian timings gave 1.00×/3.00×/5.38×/5.98× speedups at
1/4/8/16 threads. Rotation gave 1.00×/1.72×/2.05×/3.22×. Earlier Hessian variants
that shared a lambda or put dispatch inside the serial loop caused measurable
single-thread regressions; the adopted version retains a separate original
serial kernel. Peak process RSS in the three-day Hessian comparison was about
2.10 GiB before and 2.11 GiB after; this includes input buffers and GPU runtime.

The expanded frozen-kernel regression compares 205,618,624 scalar components:
all are exactly equal (`max_abs=0`, `max_relative=0`) in both Linux CPU-only and
CUDA builds. Both configurations pass all seven CTest cases. CUDA power passes
216 cases (`max_float_scaled=0.000181899`, scaled by uncorrected double power);
CUDA slant stack passes the 32-case matrix, with zero peak mismatches and maximum
scaled error `4.88329e-6`. Those GPU errors describe the existing FP32 backend,
not errors introduced by this CPU change.

Hessian and rotation independently produce byte-identical outputs for all
29 events in the five-day CUDA/bootstrap qualification. The five-day checks may
run alongside correctness builds/tests and are not performance measurements.
The combined implementation also passes five-day CPU/off and CUDA
`off`/`bootstrap`/`grid`/`all` comparisons, each with 29 events and byte-identical
event files relative to its same-backend baseline, including Bootstrap fields.
The integrated synthetic checks retain one-thread performance (Hessian 1.000×,
rotation 1.002×). Existing outer-grid timings are within 2.2% of the frozen
reference at 1/4/8/16 threads.

### Combined result

The adopted Hessian + rotation combination (`c3be328`; final test coverage in
`951a7f8`) reduces the three-day whole-run median from **39.479 s to 27.317 s**:
**1.445× throughput, 30.81% less elapsed time**. Five measured runs per binary,
excluding warmups, range from 39.381–39.592 s before and 27.166–27.367 s after.
All outputs are byte-identical. The whole-run and initial-grid gates pass.

| Stage | Before | After |
| --- | --- | --- |
| Hessian, inclusive | 13.335 s | 1.520 s |
| Rotation, inclusive | 0.820 s | 0.395 s |
| Fitting, inclusive | 14.172 s | 4.016 s |
| Initial grid, wall | 1.478 s | 1.457 s |
| Bootstrap, wall | 3.486 s | 1.823 s |

Hessian is 8.77× faster and rotation 2.08× faster in the integrated real-data
measurement. These inclusive stages overlap and must not be added. This gain
is relative to the existing CUDA version before this CPU-parallelization work;
it is not a comparison with the original unoptimized CPU-only program or Mac.

[Machine-readable measurements](benchmarks/cpu-fitting-parallel-20260914.json)
include every paired run, binary hashes, event hashes, RSS, per-stage totals,
independent experiments, the held objective change, and synthetic timings.
The ignored `build-cuda/cpu-parallel/` directory retains local fixed drivers,
logs and event files for inspection.

### Cumulative comparison with the original CPU implementation

A separate same-host comparison uses `393a93e`, before the first CPU slant-stack
optimization, against the integrated CUDA/bootstrap implementation (`c3be328`).
Both use the 3990X with 16 OpenMP threads; the current version also uses the
RTX PRO 2000 Blackwell. GCC 13.3 Release flags and inputs are identical. Only
the test drivers fix Bootstrap seeds and FFTW plans, as described above.

| Version | Whole-run median | Measured range |
| --- | --- | --- |
| Original CPU (`393a93e`) | 155.805 s | 155.647–155.992 s |
| Current CUDA + CPU optimizations | 27.416 s | 27.267–27.417 s |

The cumulative speedup is **5.683×**, or **82.40% less elapsed time**. This
includes CPU slant-stack and I/O optimization, GPU acceleration, and the latest
CPU parallelization. It is a direct paired comparison, not multiplied speedups
from different machines or workloads.

The protocol uses January 1–3, 2004 horizontal data, one excluded warmup per
binary and three measured pairs with alternating order. Profiling is disabled
in both binaries; wall time includes startup, input loading, GPU transfers and
output. No agent-started builds or tests run alongside the timings. OS caches
are not purged, so these results do not represent cold NAS throughput.

All runs retain 20 events, 433 accepted windows and identical initial candidate
choices. Each backend produces repeatable event bytes. CPU-versus-GPU outputs
pass the existing bounds: beam maximum/MAD `rtol=1e-4`, Bootstrap columns 21–33
`rtol=1e-3`, residual `atol=1e-12`, and other fields exact. Observed maximum
relative differences are `9.00657e-6` for beam maximum/MAD and `1.19376e-5` for
Bootstrap. The legacy build also passes its two CTest cases.

[The full measurement record](benchmarks/cpu-legacy-total-20260914.json) retains
all eight runs, including excluded warmups, event and binary hashes, RSS,
accepted windows and summary statistics. Its `event_file` entries refer to
local artifacts under the ignored `build-cuda/cpu-parallel/legacy-total-events/`
directory. The detached legacy checkout and generated deterministic drivers
remain under `build-cuda/`; production sources were not modified for this run.

Mac qualification was initially deferred on 2026-09-14, then partially completed
later that day. The [M4 Max Metal report](metal-cpu-parallel-performance-20260914.md)
records passing builds, regressions and three-day real-event comparisons:
Hessian + rotation parallelism gives 1.600× whole-run speedup, and the direct
pre-optimization CPU comparison gives 8.887× cumulative speedup. `main` remains
unchanged. The subsequent [five-day and one-thread Metal qualification](metal-main-qualification-20260914.md)
passes both requested follow-up checks: all 29 five-day events are byte-identical,
and five alternating one-thread pairs show no >5% regression in whole-run,
Hessian or rotation time. Real three-component event data remains unqualified;
synthetic rotation and GPU tests cover both horizontal and three-component layouts.
