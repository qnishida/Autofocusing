# Metal: CPU parallelization and cumulative speedup from before CPU optimization

Detailed benchmark logs and run results are retained locally and are not
distributed in this repository. References to detailed records below describe
those local artifacts.

[Japanese](metal-cpu-parallel-performance-20260914_jp.md)

Measured on 2026-09-14 on Apple M4 Max running macOS 26.6.2.
**Parallelizing the Hessian and RT rotation gave a 1.600× overall speedup; the cumulative
speedup from before CPU optimization was 8.887×.**
The first compares code changes with identical Metal settings; the second compares the
old CPU version with the current Metal version.

## Conditions and versions

- Input: Hi-net tilt horizontal components, 2004-01-01 through 03; 20 events and 433 accepted windows per trial.
- Compilers: Homebrew LLVM/Clang 23.1.1 (C/C++), Apple Clang (Objective-C++).
  All versions used Release `-O3 -DNDEBUG`, native architecture disabled, 16 OpenMP threads, and dynamic teams disabled.
- GPU settings: `AUTOFOCUSING_BACKEND=metal`, `AUTOFOCUSING_GPU_POWER=bootstrap`.
  The initial grid runs on the CPU. The power setting matches the existing CUDA comparison.
- Separate processes ran in alternating order on the same machine. One initial warmup per version was excluded.
  OS caches were not cleared, and no other agent-started builds or tests ran during measurement.
- Only the validation driver fixes FFTW to `ESTIMATE | UNALIGNED` and the bootstrap seed to
  `1837 + 104729 * replicate`. Production source was unchanged.
- Wall time includes startup, input, GPU transfers and synchronization, and output.

| Version | Commit | Purpose |
| --- | --- | --- |
| Before CPU optimization | `393a93e030d74ace231dbd46ea85175a0795f59a` | CPU-only baseline for the cumulative comparison |
| Before Hessian and rotation parallelization | `6a7a4bd4d40a20f1383f683fe233776aed674efa` | Existing Metal implementation with shared profiling added; numerical kernels unchanged |
| Current gpu | `dabb03727454421f4fbb0b42a61604539fe99a6d` | After parallelization; target of both comparisons |

`393a93e` predates the first CPU slant-stack optimization and is the same baseline used
for the cumulative CUDA comparison. `6a7a4bd` is the same baseline used for the CUDA
CPU-parallelization comparison.

## (i) Effect of Hessian and RT-rotation parallelization

Both versions use Metal/Bootstrap with profiling enabled. Values are medians of five runs per version.

| Operation | Before parallelization | After parallelization | Speedup |
| --- | ---: | ---: | ---: |
| Total wall time | 21.039 s | 13.150 s | **1.600×** |
| Hessian | 8.945 s | 1.417 s | **6.311×** |
| RT rotation | 0.342 s | 0.103 s | **3.312×** |
| Fitting | 10.083 s | 3.592 s | 2.807× |
| Bootstrap stage | 1.702 s | 0.660 s | 2.577× |
| Initial-grid stage | 1.144 s | 1.060 s | 1.080× |

Total elapsed time fell by **37.50%**. Measured ranges were 20.980–21.105 s before
parallelization and 13.099–13.430 s after. Both total time and initial-grid time passed
the regression checks. Hessian and rotation each met the criterion of at least a 10%
reduction in the target stage.
Hessian, rotation, and fitting use `serial_caller` inclusive time and overlap, so they
must not be summed. This breakdown measures the two adopted changes together, rather
than experiments enabling each change separately.

Event output was byte-identical in all 12 trials, including warmups.
Accepted windows and initial candidates also matched. The comparison treats this as
a change within the same backend and does not relax the bootstrap tolerance for GPU error.

All trials, hashes, and stage timings
are saved alongside the logs and event files for each trial.

## (ii) Cumulative effect from before CPU optimization

The old version uses CPU/off and the current version uses Metal/Bootstrap.
Profiling is disabled for both. Values are medians of three runs per version.

| Version | Median total time | Measured range |
| --- | ---: | ---: |
| Before CPU optimization, `393a93e` | 121.214 s | 120.622–121.527 s |
| Current Metal + CPU optimizations | 13.639 s | 13.443–13.660 s |

**8.887× speedup, with an 88.75% reduction in elapsed time.** This direct comparison on
the same machine includes CPU slant-stack optimization, shared I/O, Metal acceleration,
and the Hessian and rotation parallelization measured here.
It is not a product of speedup ratios from separate measurements.
The current-version median is reported separately because the trials and profiling
conditions differ from (i).

All eight trials, including warmups, matched in their 20 events, 433 accepted windows,
and initial candidates. Event output was byte-identical across repetitions of each version.
All 38 columns were checked across CPU and GPU using the existing tolerances:

- Beam maximum and MAD (columns 12–13): `rtol=1e-4`, maximum relative difference `1.5593e-5`.
- Bootstrap (columns 21–33): `rtol=1e-3`, maximum relative difference `2.8718e-5`.
- Convergence residual (column 15): `atol=1e-12`.
- Other columns: exact output-string equality. Existing NaNs for missing U in horizontal mode are preserved;
  non-finite values are rejected in columns compared numerically.

All trials, hashes, accuracy, and accepted windows
are saved alongside the logs and event files for each trial.

## Build and reproduction

Prepare each commit in a separate complete checkout. Use the same build configuration
as the [earlier Metal verification](gpu-metal-verification-20260914.md).
Select Apple Clang for Objective-C++. The old CPU version does not need Metal build options.
The current `build_io_probe.py` can also be copied into the old version's `tests/` directory.
Do not replace other reference-version sources or numerical kernels with current versions.

```bash
# Run after building each checkout.
python3 tests/build_io_probe.py build

# From the gpu checkout root, specify an output directory that does not exist.
python3 tests/run_cpu_parallel_events.py \
  "$BEFORE_PROBE" "$CURRENT_PROBE" "$HINET_ROOT" "$CMT_CATALOG" \
  parallel-events --backend metal --power bootstrap --threads 16 --days 3 --repeats 5

# Reusable comparison script added during this validation.
python3 tests/run_legacy_gpu_events.py \
  "$LEGACY_PROBE" "$CURRENT_PROBE" "$HINET_ROOT" "$CMT_CATALOG" \
  cumulative-events --backend metal --power bootstrap --threads 16 --repeats 3
```

CTest passed for all three versions: 2/2 for the old version, 7/7 before parallelization,
and 7/7 for the current version. The current version also passed 32 Metal slant-stack cases,
216 power cases, and batch-boundary tests.
The first launcher test for the old version failed because the archive copy lacked Git
metadata; it passed when rerun with Git metadata for the target commit. No source fix was needed.
Metal and `/usr/bin/time -l` ran outside the agent sandbox on this host.
Build and CTest logs
are also saved. The comparison scripts were checked with synthetic cases for identical
output, bootstrap tolerance boundaries, event-identity mismatches, and non-finite-value rejection.

The existing CUDA record is at `gpu@dabb037:docs/cpu-fitting-parallel.md`.
On 3990X + RTX PRO 2000, CUDA achieved 1.445× from parallelization and 5.683× from before
CPU optimization. The 1.600× and 8.887× values here compare versions on M4 Max;
they are not a direct CUDA-versus-Metal performance comparison.
Real three-component data, a five-day before/after comparison, the full archive, and
additional single-thread measurements were outside this measurement scope.
These results do not mean that every main-merge criterion has been completed.

The follow-up [five-day and single-thread Metal validation](metal-main-qualification-20260914_en.md)
confirmed byte-identical output over five days before and after parallelization and passed
the 5% regression checks for single-thread total time, Hessian, and rotation.
Integration into main has not yet been performed.

The working branch `cuda` was retained. No production code edits, installation, commits,
or pushes were performed.
