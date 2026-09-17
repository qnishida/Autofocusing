# CUDA implementation and validation

Detailed benchmark logs and run results are retained locally and are not
distributed in this repository. References to detailed records below describe
those local artifacts.

CUDA was added to `cuda` from `b9620e7`, retaining the CPU and Metal backends.
The common dispatch accepts `AUTOFOCUSING_BACKEND=cpu|metal|cuda`; CPU is the
default and explicit GPU requests fail when unavailable.

## Implementation

Slant stacking uses the Metal equations and component/station/frequency layout:
a frequency-wise self/cross-power correction, followed by one 256-thread block
per slowness cell. Station sums remain ordered; frequency contributions use a
shared-memory tree reduction. Horizontal mode leaves U untouched. GPU powers
are checked for finiteness, promoted to double and added to the existing output.

The optional power objectives use separate phase, weighted-power
and reduction kernels. At most 32 candidates are processed at once; spectra
and coordinates transfer once per batch API call. Shared geometry/weights use
the existing interpretation, including Bootstrap's common phase. Potential
grid winners retain the existing `2e-3 * abs(max) + 1e-30` band for parallel CPU
double re-evaluation and original strict-max selection order.

CUDA uses reusable device buffers and pinned host staging buffers, explicit
copies, and one stream per backend state. A mutex serializes access to each
state; RAII releases resources and drains work on exceptions. Allocation,
copy, launch and synchronization failures are reported. Dimensions and FP32
conversion are checked before use. CUDA compilation disables fast math, FMA
fusion and flush-to-zero and requests precise division/square root. This is
FP32 computation for horizontal objectives, with the existing `5e-4` validation
bounds. Three-component Bootstrap uses the same kernels instantiated for FP64;
its synthetic checks use `1e-12` scaled by uncorrected CPU double power plus
the existing `1e-30` absolute floor.

FFTW, HDF5 I/O, production seeds, sampling, final fitting/derivatives and CPU
objectives are unchanged. GPU Bootstrap supports horizontal and 3c input;
GPU initial-grid objectives remain restricted to horizontal input.
Metal's existing source and math settings are retained; its legacy power
environment variable remains available as described in the manual.

## Local checks (2026-09-13)

Ubuntu 24.04, GCC 13.3, CMake 3.28.3, CUDA Toolkit 13.2.86, driver 595.84,
RTX PRO 2000 Blackwell (`sm_120`), Ryzen Threadripper 3990X. Synthetic benchmark
comparisons use 16 OpenMP threads. All sources are built in Release mode.

- Unmodified CPU baseline: original CTest 4/4 passed.
- CPU-only and CUDA-enabled builds: CTest 6/6 passed, including backend settings,
  legacy precedence, unavailable builds and the real-event comparison driver.
- Slant-stack synthetic regression: all 32 horizontal/3c, random/coherent,
  centered/offset, 1/31/155/257-bin conditions pass. Maximum scaled error
  `4.89e-6`, maximum relative L2 error `3.77e-6`, no peak mismatches.
  Accumulation, empty-frequency/one-station guards and verification dispatch
  are exercised, as are invalid CUDA dimensions.
- Power regression: 216 cases pass; maximum error scaled by uncorrected CPU
  double power is `1.819e-4` (limit `5e-4`). The 65-candidate tests cover all four
  geometry/weight-sharing combinations across chunk boundaries, plus zero/tie
  behavior and invalid/overflow dimensions.
- Compute Sanitizer memcheck: zero errors for both GPU test executables.
- Metal GPU regression is not runnable on this Linux host.

The working tree initially lacked several tracked CPU sources and its velocity
model. Initial validation used a complete temporary checkout of the base commit
plus the CUDA changes. After explicit approval, the 11 missing tracked files
were restored from HEAD; the Git-root CUDA build, CTest 6/6 and both GPU
regressions also passed. The pre-existing `repo/` directory was left untouched.

## Synthetic performance

One warmup, five alternating CPU/GPU calls. These are per-call mean times;
CUDA includes packing, allocation/reuse, transfers, synchronization and double
accumulation. Device initialization was separately about 0.16–0.17 seconds.

| Grid | Horizontal CPU / CUDA | Speedup | 3c CPU / CUDA | Speedup |
| --- | ---: | ---: | ---: | ---: |
| 67×67 | 99.36 / 4.72 ms | 21.06× | 127.95 / 6.51 ms | 19.64× |
| 133×133 | 386.13 / 17.07 ms | 22.63× | 495.25 / 23.74 ms | 20.86× |
| 201×201 | 884.09 / 38.90 ms | 22.73× | 1135.99 / 54.07 ms | 21.01× |

All six benchmark peak sets agree with CPU; maximum scaled error is `2.50e-5`
or less. The input has 650 active stations, 653 allocated slots and 155 bins.
Grid spacings are .005, .0025 and .0025 s/km respectively.

For initial fitting grids (650 stations, 3 windows, five trials), median
CPU/CUDA speedups are 1.04–1.08× for 75 candidates, 3.10–3.33× for 750, and
3.99–4.53× for 7500. Each range covers finer and wider search variants. All
selected winners agree. The standard synthetic grid does **not** meet the
10% stage-reduction adoption condition; this must be assessed separately on
real events. Tiny slant-stack inputs can also be slower on GPU.

Raw synthetic timings, precision metrics and source hashes are in
the machine-readable record.

## Real-event qualification and reproduction

The supplied horizontal tilt data and CMT catalog pass the five-day
2004-01-01–05 comparison: CPU and CUDA detect the same 29 events; CUDA
`off`, `bootstrap`, `grid` and `all` preserve initial candidates and event
identity. Maximum Bootstrap relative difference is `1.19e-5`. The unchanged
CPU path is byte-identical to the deterministic `b9620e7` CPU driver across
all 29 events. A separate instrumented run checks all 740 windows: maximum
scaled error `9.993e-6`, relative L2 error `6.332e-6`, zero peak mismatches,
and byte-identical event output to CUDA/all. Peak reusable device capacities
are 1,590,904 bytes for slant stacking and 70,634,472 bytes for power; these
exclude the CUDA context/driver and host staging memory. Days 4–5 were not used to adjust any tolerances or refinement
settings. Qualification timings overlapped other validation/build activity
and are not adoption measurements. See the
qualification record.


### Repeated real-data performance

Three days (2004-01-01–03), 20 events per run, 16 CPU threads, five fresh
processes per mode in alternating order. No other agent builds, checks or
benchmarks ran alongside these measurements; OS caches were not purged.
Every run passed output checks, including CPU/GPU comparisons. The wall
clock includes initialization, input loading, packing, transfers and sync.
Median seconds:

| Backend / power mode | Whole run | Slant stack | Bootstrap | Initial grid |
| --- | ---: | ---: | ---: | ---: |
| CPU | 108.266 | 40.242 | 34.020 | 1.464 |
| CUDA / off | 69.747 | 1.967 | 34.051 | 1.454 |
| CUDA / bootstrap | 39.141 | 1.976 | 3.466 | 1.451 |
| CUDA / grid | 70.344 | 1.981 | 33.023 | 3.000 |
| CUDA / all | 40.742 | 1.982 | 3.383 | 2.987 |

**Recommended: `AUTOFOCUSING_BACKEND=cuda` with
`AUTOFOCUSING_GPU_POWER=bootstrap`.** This is 2.77× faster overall than CPU.
Slant stacking alone is 20.46× faster at its stage; Bootstrap is 9.82× faster
than CUDA/off's CPU Bootstrap stage. Both pass the 10% stage-reduction and
5% maximum whole-run regression gates.

CUDA initial-grid evaluation is accurate but fails the stage-speed gate
(3.000 s versus 1.454 s). It remains available as an explicit experimental
option; it is not recommended for the current standard grid. `all` passes
the aggregate grid+Bootstrap gate but is slower than `bootstrap` alone.
The default backend remains CPU and the default GPU power mode remains off.

Peak process RSS across the five runs (includes host memory, not standalone VRAM):

- cpu: 1.391 GiB.
- off: 1.926 GiB.
- bootstrap: 2.099 GiB.
- grid: 2.059 GiB.
- all: 2.086 GiB.

Full trial ranges, stage timings, event hashes and adoption flags are in the
performance record.
These results cover the supplied five horizontal days. Real 3c events and the
full archive were not qualified; 3c coverage here is synthetic slant stacking.

### Commands

Use an existing complete checkout, run from its root (for the velocity model),
and provide a new output directory for each run:

```bash
python3 tests/build_io_probe.py build-cuda
python3 tests/run_metal_power.py \
  build-cuda/io-fixed-current/cal_ccf_io_probe "$HINET_ROOT" "$CMT_CATALOG" \
  build-cuda/cuda-qualification --backend cuda --threads 16 \
  --days 5 --repeats 1 --compare-cpu
python3 tests/run_metal_power.py \
  build-cuda/io-fixed-current/cal_ccf_io_probe "$HINET_ROOT" "$CMT_CATALOG" \
  build-cuda/cuda-timings --backend cuda --threads 16 \
  --days 3 --repeats 5 --compare-cpu
```

The historical script filename is retained for compatibility. It supports
Linux RSS measurement, selected GPU/off baselines, CPU/off comparisons,
alternating mode order, explicit thread counts and optional expected event
counts. It never assumes that CUDA reproduces Metal's 20/29 event counts.
CPU/GPU event comparison uses the existing beam maximum/MAD `rtol=1e-4` and
residual `atol=1e-12` rules. Power-mode comparisons require the same initial
candidates and event identities; only Bootstrap/all allow columns 21–33 to
vary by `rtol=1e-3, atol=1e-30`. New nonfinite Bootstrap values fail validation.
The fixed seed and FFTW plan flags apply only to the test driver.

GPU kernel-only timing is available separately with
`AUTOFOCUSING_CUDA_PROFILE=1`. Profiled or CPU-verification runs must not be
used for adoption timing. Real-event acceptance requires unchanged event
outputs under the documented comparison and a target-stage reduction of at
least 10%, with total runtime no more than 5% worse than its baseline.

## Three-component Bootstrap, 2026-09-17

CUDA `bootstrap` and `all` now evaluate three-component Bootstrap on the GPU.
The event already supplies one selected R, T or U spectrum, so the batched
objective layout is unchanged. Initial-grid fitting remains on CPU in 3c mode;
`off` and `grid` retain the CPU OpenMP Bootstrap implementation. Peak search,
sampling order, seeds, sample count, Hessians and final statistical reductions
are unchanged. This is opt-in; the default power mode remains `off`.

Simply enabling the existing FP32 path failed the established real-event bound:
one T event on 2005-03-15 had Bootstrap mean power `8.05874e-21` instead of
`8.03698e-21`, exceeding `rtol=1e-3`. No tolerance was relaxed. Three-component
CUDA power now uses FP64 spectra, coordinates, weights, phase evaluation,
partial sums, bias subtraction and output. Horizontal CUDA objectives retain
FP32. Metal three-component Bootstrap retains CPU OpenMP evaluation.

On Threadripper 3990X / RTX PRO 2000 Blackwell, both CPU-only and CUDA builds
pass all 12 CTest cases. The CUDA objective test passes 432 FP32/FP64 cases,
with maximum errors scaled by uncorrected CPU power of `1.819e-4` and
`3.088e-14`, respectively. An additional 1,212 FP64 values cover three selected
component spectra, 1/3/17/48 windows and 101-entry Bootstrap batches spanning
32-entry chunks; maximum scaled error is `1.263e-14`. The FP64 bound is
`1e-12 * uncorrected_power + 1e-30`: the existing CPU relative coefficient is
applied to the same uncorrected-power scale as the GPU tests, avoiding division
by a nearly cancelled corrected power.
Compute Sanitizer memcheck reports zero errors.

The fixed-seed Denoise 2005-03-15–17 comparison uses CUDA slant stacking in both
runs and compares CPU OpenMP Bootstrap (`off`) with GPU FP64 Bootstrap
(`bootstrap`). All 59 event records (28 R, 16 T, 15 U), including covariance,
Bootstrap and spectral matrices, are byte-identical. Initial candidates and
all 520 accepted windows match. These are correctness runs, not adoption
timings; other correctness checks ran during qualification.

### Repeated one-day timing

The Denoise 2004-04-01 workload (6 events, 137 windows, 0.1–0.25 Hz,
67×67 slant-stack grid) uses 16 OpenMP threads, dynamic teams disabled, fixed
seeds/FFTW plans, one excluded warmup per mode and five alternating fresh-process
measurements. CUDA slant stacking is enabled in both modes. Caches are not
purged; no other agent-started builds or tests ran alongside these measurements.

| Stage | CPU OpenMP Bootstrap (`off`) | GPU FP64 Bootstrap (`bootstrap`) |
| --- | ---: | ---: |
| Whole-run median | 14.240 s | 13.891 s |
| Bootstrap median | 0.855 s | 0.551 s |
| Initial-grid median | 0.356 s | 0.362 s |
| Slant-stack median | 0.811 s | 0.811 s |

Bootstrap is 1.551× faster and whole-run time decreases by 2.45%. Whole-run
ranges are 14.090–14.292 s and 13.842–13.941 s. The target-stage and whole-run
gates pass. All 12 event files, including warmups, are byte-identical, with
matching candidates and windows. This measures the incremental Bootstrap gain
over OpenMP with GPU slant stacking already enabled, not the combined CPU-to-GPU
speedup. Larger Bootstrap workloads may give a different whole-run benefit.

Reproduce with the deterministic driver and a fresh output directory:

```bash
python3 tests/build_io_probe.py build-cuda
python3 tests/run_metal_power.py \
  build-cuda/io-fixed-current/cal_ccf_io_probe "$HINET_ROOT" "$CMT_CATALOG" \
  build-cuda/bootstrap-timing-new --backend cuda --components 3c \
  --start-date 2004-04-01 --days 1 --threads 16 \
  --modes off bootstrap --repeats 5 --warmups 1 --expected-events 6
```

For the independent 59-event qualification use `--start-date 2005-03-15
--days 3 --repeats 1 --warmups 0 --expected-events 59`. Detailed logs, event
files and JSON reports remain in the ignored build directory.
