# CPU slant-stack optimization before Metal

Implemented on `perf/cpu-slant-stack`, based on `393a93e`. These measurements were made before creating a `metal` branch. Slant-stack arithmetic stays double/complex<double>;
FFT/QC, nonlinear fitting and bootstrap still run on the existing CPU paths.
No GPU implementation, precision reduction or slowness-grid change is included.

## Measurement and outcome

The real-data comparison uses only 2004-01-01 (653 loaded stations), the same
CMT catalog and grid, and sequential runs on this Mac. Baseline is `393a93e`
with only the same stage timers added. Builds use Homebrew LLVM/Clang, Release
`-O3`, C++14, and no fast-math. OS caches were not purged. These are local
measurements, not cold-disk or multi-day throughput guarantees.

| Stage | Baseline, 4 threads | Optimized, 4 threads | Optimized, 16 threads |
| --- | ---: | ---: | ---: |
| FFT and QC | 0.388 s | 0.402 s | 0.193 s |
| Packing | 0.061 s | 0.068 s | 0.068 s |
| Slant stack and power reduction | 75.285 s | 26.632 s | 8.549 s |
| Candidate search and subsequent fitting | 3.972 s | 4.003 s | 3.813 s |
| Total segment processing | 79.770 s | 31.146 s | 12.667 s |

At the same four threads, slant stack is **2.83× faster** and segment processing
is **2.56× faster**. The 16-thread column combines code improvements with more
threads and must not be interpreted as an optimization-only speedup. Totals
exclude input loading and startup. Initialization plus reading took about
2.15 s in the baseline; the initial bottleneck was slant stacking.

Kernel-only comparisons use 650 active stations, 155 frequency bins and a
67×67 grid. Each row averages three calls; old and new paths see the same
synthetic spectra and coordinates. The reference preserves the original
Boost-array loops, including the original always-three-component work.

| Threads | Horizontal old/new | Horizontal speedup | 3c old/new | 3c speedup |
| --- | ---: | ---: | ---: | ---: |
| 4 | 0.628 / 0.221 s | 2.85× | 0.627 / 0.262 s | 2.40× |
| 8 | 0.324 / 0.113 s | 2.88× | 0.322 / 0.133 s | 2.41× |
| 16 | 0.198 / 0.069 s | 2.85× | 0.195 / 0.083 s | 2.36× |

The benchmark reference allocates its legacy grid scratch for each call,
whereas the old application allocated it once per segment and cleared it
per window. Use the real-day comparison above as the end-to-end evidence.
Detailed measurements, including rejected experiments, are stored in
[the JSON report](benchmarks/cpu-slant-stack-20260908.json).

## Changes retained

- Generate a station/grid phase recurrence once and apply it to all active
  components. Horizontal mode does not compute the zero U stack.
- Give each grid cell one OpenMP owner, which performs both complex summation
  and R/T/U power reduction. There is no cross-thread floating-point reduction.
- Replace the full-grid complex intermediate with small per-thread scratch.
  The old default-grid array occupied 33,398,160 bytes (31.85 MiB). New
  horizontal sum scratch is 4,960 bytes per worker, plus shared per-frequency
  power arrays. This excludes the retained input/output and downstream buffers.
- Traverse adjacent py cells together and use `schedule(dynamic, 8)` to reduce
  scheduling overhead while balancing workers. Static scheduling was slower
  at higher thread counts in these measurements.
- Expose a compact internal CPU kernel interface in `src/slant_stack.h` with
  contiguous ENU input, active station count, grid parameters and RTU output.
  It preserves the existing caller's accumulation semantics. This is also a
  clear boundary for a later Metal implementation; it is not a GPU backend.
- Add `AUTOFOCUSING_PROFILE=1` stage logs and a reusable comparison script.

The mathematical expression and the order of station/frequency accumulation
within each cell are retained. Slowness increment, range, source selection,
CMT masking, rotation and fitting rules were not altered.

## Validation

- CTest: **3/3 passed**, including existing horizontal input and launcher tests.
- Slant-stack regression: both modes; one/four threads; offset and centered
  grids; unused station capacity; repeated accumulation; fewer-than-two-station
  guard. The production-size and enlarged-grid benchmarks also compare every
  cell to the independent legacy reference.
- Maximum absolute cell difference in measured cases: **0**. This is an
  observation on this build; tests allow scale-relative error up to 1e-12 plus
  1e-30 for portability. No general bitwise guarantee is claimed across
  compilers or different FFTW plans.
- Real day: all versions accept **47, 19, 26, 42** windows (134 total) and
  output the same two R events. Event times, locations, slownesses, deterministic
  powers and flags match at output precision. U-only fields remain `nan`.
- The near-zero optimizer residual is compared with absolute tolerance 1e-12.
  Existing bootstrap results use a time-based random seed, so covariance and
  bootstrap-power columns 21–32 are excluded from cross-run identity checks.
- The earlier valid 3c spectrum/QC comparison against `149995f` still passes
  for 1599 values at rtol=1e-10, atol=1e-18. No real 3c catalog regression was run.

## Reproducing checks and profiling

Build with the Clang setup in [the manual](../manual.md), then:

```bash
cmake --build build-clang -j 4
ctest --test-dir build-clang --output-on-failure
build-clang/test_slant_stack --benchmark 650 33 4 3
```

The benchmark arguments are active stations, grid half-width, threads and
repetitions. A half-width of 33 gives 67×67 cells; 67 gives 135×135. The latter
was checked in a static-scheduling experiment, with exact cell agreement and
roughly four times the runtime as expected from the cell count. This does not
change the production grid.

For a bounded real run, create an input tree containing only the day being
measured. The local verification tree `build-gcc/one-day/` contains one
read-only symlink to the original HDF5 file; its name does not indicate the
compiler used for execution. Use a fresh result ID for each run:

```bash
AUTOFOCUSING_PROFILE=1 OMP_NUM_THREADS=4 build-clang/src/cal_ccf_clang \
  2004 tilt_horizontal cpu-check build-gcc/one-day ../moment_loc_76_24 \
  build-clang/performance/results horizontal \
  > build-clang/performance/cpu-check.log 2>&1
python3 Scripts/compare_profiles.py baseline.log cpu-check.log
python3 tests/compare_event_results.py baseline.dat cpu-check.dat
```

The parent launcher also inherits `AUTOFOCUSING_PROFILE` and `OMP_NUM_THREADS`.
Its normal archive path scans the full configured period, so use the bounded
input tree for benchmarking. The optimized executable has been installed to
`bin/cal_ccf_clang` during verification. The parent configuration was unchanged;
the CPU changes were subsequently approved for merging into local `main`.

## Deferred work and decisions

A separate phase buffer was slower (four-thread horizontal 0.341 s versus
0.221 s for the retained path), so it was rejected. No fast-math, explicit
low-precision arithmetic or full station/grid/frequency phase cache was added.

The existing macOS prefetch function does not issue the POSIX hint because
`POSIX_FADV_WILLNEED` is unavailable in this build; it is not an asynchronous
HDF5 loader. Adding real overlap would require a consecutive-day I/O benchmark,
a defined buffer/memory budget and review of HDF5 thread safety. It was deferred
because the measured slant stack dominated the baseline. Loading and the
roughly four-second fitting stage will matter more after further acceleration.

The statistical concerns documented in
[the earlier verification notes](horizontal-verification.md) are unchanged.
They should be handled separately from a performance-equivalence change.

## Ablation: phase sharing versus intermediate buffer

A follow-up comparison restores one feature at a time in the optimized
kernel. These are synthetic-data kernel timings: 650 stations, 67×67 cells,
155 bins, four threads, means of five calls. Scheduling (`dynamic, 8`), the
contiguous cell order, and skipping U in horizontal mode remain fixed.

| Variant | Horizontal | 3c |
| --- | ---: | ---: |
| Final optimized kernel | 0.218 s | 0.261 s |
| Repeat phase generation separately for each active component | 0.413 s | 0.614 s |
| Restore full-grid complex buffer and separate serial power reduction | 0.221 s | 0.271 s |
| Restore both | 0.438 s | 0.633 s |

All cases had maximum absolute cell error 0 versus the legacy reference.
At this size, sharing the phase recurrence and combining the component work
is the dominant improvement. Restoring the intermediate array and the serial
reduction had a much smaller timing effect. This supports attribution to
phase sharing/component-loop restructuring, not specifically to the initial
sin/cos calls alone. The recurrence across all frequency bins is also shared.

The buffer experiment restores allocation, storage, and separate reduction
together, so it does not isolate pure allocation cost. Small differences of
a few milliseconds should not be overinterpreted. Avoid treating the gains
as additive percentages. Intermediate-buffer removal still reduces memory
footprint as the grid grows; larger-grid timing attribution was not measured
in this ablation. Horizontal U skipping is common to all ablation variants
and contributes separately to the comparison with the old always-3c loop.

Reproduce without editing production code:

```bash
python3 tests/ablate_slant_stack.py build-clang 4 5
```

This generates temporary variants under the build directory, compiles using
the existing compilation database, and compares each against the legacy
reference. Raw summaries are included in the JSON report above. The production
kernel was not changed during this experiment.
