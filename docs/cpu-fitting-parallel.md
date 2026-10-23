# Ordered CPU fitting parallelism

The CPU fitting kernels are shared by CPU, CUDA and Metal builds. See
[current behavior](../README.md#current-behavior-and-limits) for the adopted configuration.

## Implementation and invariants

The Hessian distributes independent time-window stacks over OpenMP workers.
It retains at most `min(32, window_count, omp_get_max_threads())` windows of
intermediate results, then accumulates them
serially in the original window/frequency/component order. Station order,
phase recurrence and double-precision accumulation are unchanged. One-thread,
fewer-than-two-window and already-parallel callers use the original serial kernel.
There is no OpenMP floating-point reduction or nested worker team.

RT rotation partitions window/station pairs. Each worker owns separate output
spectra and weights; E/N/U inputs and scalar expressions are preserved.
The scalar objective `cal_S` remains serial. An objective-parallelization
experiment was not adopted because it slowed the initial-grid path.

Three-component CPU Bootstrap evaluates independent resamples in batches of
at most 16 workers. Resampling and final statistics retain their original order.
Production seeds remain time-dependent; fixed seeds and FFTW plans belong only
to validation drivers. The adopted Newton iteration uses FP32 4x4 operations,
while the covariance inverse uses FP64; see the
[precision record](newton-fp32-rollback-20260917_en.md).

## Hessian scratch storage

`DELTAP_ENABLE_HESSIAN_PADDING` defaults to `ON` on Linux and `OFF` on other
platforms. It adds 31 unused frequency elements to each parallel Hessian scratch
row. Active bins, arithmetic and accumulation order are unchanged. Use
`-DDELTAP_ENABLE_HESSIAN_PADDING=OFF` to retain the unpadded layout.

Hessian's limit is an upper bound, not a fixed team size. A 24-window Hessian
can use 24 workers when `OMP_NUM_THREADS` permits it, while smaller teams reserve
only their needed scratch slots. CPU Bootstrap remains capped at 16 workers.
The build uses the system allocator on all platforms.

## Verification protocol

`test_cpu_parallel` compares objectives, gradients, Hessians, rotated spectra
and weights against independently frozen serial kernels. Coverage includes
zero/nonuniform weights, coherent/random data, horizontal/3c modes, partial
batches, 24/32-thread Hessians and already-parallel callers. The componentwise bound is
`1e-12 * abs(reference) + 1e-30`; exact equality is additionally reported.

`AUTOFOCUSING_PROFILE=1` reports inclusive objective, Hessian, rotation and
fitting times. `serial_caller` describes the calling context, including any
internal parallel work. `outer_parallel` sums worker-call durations and is
not pipeline wall time. Fitting contains objective/Hessian work, so these
categories must not be added together.

Real-event comparisons use the same backend, inputs, QC history, seed and FFTW
controls. Compare every output field at printed precision, with only the existing
absolute `1e-12` exception for the convergence residual. CPU/GPU error bounds
are not a substitute for same-backend implementation-regression checks.

## Reproduction

Configure a Release build as described in the manual, then run:

```bash
cmake --build build --target test_cpu_parallel -j 2
OMP_NUM_THREADS=4 ctest --test-dir build -R '^cpu_parallel_equivalence$' --output-on-failure
build/test_cpu_parallel --benchmark hessian 16 5 24
build/test_cpu_parallel --benchmark objective 16 5 24
build/test_cpu_parallel --benchmark rotation 16 5 24
```

The benchmark uses synthetic data. Real-waveform checks require external
waveforms and a catalog. The public helpers `tests/build_io_probe.py` and
`tests/run_cpu_parallel_events.py` provide deterministic same-backend comparisons;
use their help and a new output directory outside the source checkout.

## Results

Historical Linux measurements used a Threadripper 3990X, RTX PRO 2000 Blackwell,
GCC 13.3 Release, native-architecture optimization disabled and 16 threads.
The input was January 1–3, 2004, horizontal data; caches were warm. Each
implementation had one excluded warmup and five alternating measured runs.
These are comparisons of the stated implementations, not cold-storage tests.

| Independent change | Whole-run median, before → after | Target stage, before → after | Decision |
| --- | ---: | ---: | --- |
| Hessian (`989a4e4`) | 39.343 → 27.817 s | 13.298 → 1.552 s | Adopted |
| Objective (`ff619dd`) | 39.389 → 37.886 s | 2.469 → 0.432 s | Not adopted: grid time increased 15.83% |
| Rotation (`d7df22e`) | 39.342 → 38.795 s | 0.830 → 0.375 s | Adopted |

All five pairs in each row produced byte-identical event outputs. The combined
Hessian/rotation implementation (`c3be328`) reduced the whole-run median from
39.479 to 27.317 s: 1.445x throughput and 30.81% less elapsed time.
Five-day checks covered 29 events with byte-identical same-backend outputs.
These results apply to the tested revisions and inputs.

## Peak-search integration validation, 2026-09-17

The three-component Bootstrap integration compared corrected serial and parallel
implementations on March 15–17, 2005 with eight CPU threads, fixed seeds and
FFTW controls. All 145 events and 147 initial candidates matched; selection
retained 75 events. Fourteen selected CTests passed.

On M4 Max, a separate March 15 CPU timing comparison used one excluded warmup
and five alternating trials per version. All 12 outputs had the same 46 events.
Median whole time fell from 110.92 to 56.14 s; Bootstrap fell from 64.73 to
10.05 s. These results qualify that workload, not arbitrary events or GPU paths.

See the [Metal comparison](metal-cpu-parallel-performance-20260914_en.md) and
[five-day qualification](metal-main-qualification-20260914_en.md) for backend-specific
validation summaries. Detailed run records and machine-specific procedures are
maintained separately from the public source distribution.
