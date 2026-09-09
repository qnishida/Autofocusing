# Metal slant-stack backend

Implemented on `metal`, based on CPU commit `51f60a7`. The default runtime
backend remains double/OpenMP CPU. `AUTOFOCUSING_BACKEND=metal` selects float
GPU slant stacking for either horizontal or three-component spectra. Loading,
rotation, FFT/QC and subsequent fitting retain their existing double CPU paths.
The installed `bin/cal_ccf_clang` includes both backends. Parent workspace
configuration was not changed; enable Metal explicitly as described in
[the manual](../manual.md#11-metal-gpu-slant-stacking).

## Implementation

The backend packs only active stations/components into reusable shared Metal
buffers, converting complex<double> to float2. A GPU kernel computes station
self-power/cross-power corrections once per frequency. Another assigns one
256-thread group to each slowness cell, with frequencies distributed among
threads. Each thread sums stations in order, shares the phase rotation across
E/N/U, and computes corrected R/T/U power. A threadgroup reduction sums power
across frequencies. Frequencies beyond 256 use another iteration; partial
frequency groups are zero-filled for reduction. U input and output are omitted
in horizontal mode.

GPU phase factors use sine/cosine at each frequency instead of the CPU's serial
frequency recurrence. Safe Metal math is requested, and neither the CPU
compiler's fast-math options nor downstream arithmetic are changed. Host code
waits for successful command completion, checks finite output, promotes float
powers to double, then adds them to the caller's existing array. The existing
daily/segment accumulation and normalization semantics are preserved.

Device, queue, pipelines and buffers persist across windows; calls are
serialized around shared buffers. Shader source is embedded in the executable
and compiled at startup using Apple's
[source-library API](https://developer.apple.com/documentation/metal/mtldevice/makelibrary(source:options:)).
This avoids an external shader file or the offline `metal` compiler. Explicit
Metal requests fail on unavailable devices, shader/command errors, or builds
without Metal. There is no implicit CPU fallback.

## Measurements

Measured on Apple M4 Max, Homebrew LLVM/Clang Release `-O3`, 16 OpenMP threads,
2026-09-08. CPU comparator is the optimized double kernel from `51f60a7`.
GPU times include float packing, buffer allocation/reuse, submission, waiting
and double accumulation. Kernel benchmarks warm both paths once and average
five calls. One-time device/shader initialization was about 0.04 seconds with
warm driver caches; the first probe was about 0.15 seconds. Cold compilation,
other GPUs and sustained thermal behavior were not benchmarked.

Synthetic spectra: fixed seed 1837, 650 active stations, 653 allocated stations,
155 bins (102–256), station coordinates within ±1500 km. Each mode uses the
same inputs for CPU and GPU.

| Grid | Step / maximum (s/km) | Horizontal CPU / GPU | Speedup | 3c CPU / GPU | Speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| 67×67 | 0.005 / 0.165 | 69.12 / 7.34 ms | 9.42× | 82.22 / 8.07 ms | 10.18× |
| 133×133 | 0.0025 / 0.165 | 269.43 / 19.46 ms | 13.84× | 321.40 / 27.49 ms | 11.69× |
| 201×201 | 0.0025 / 0.25 | 617.77 / 40.96 ms | 15.08× | 769.20 / 56.97 ms | 13.50× |

Real-data comparison uses Hi-net tilt-derived velocity on 2004-01-01, 653
loaded stations, the existing CMT catalog and horizontal mode. A bounded
input tree contains only that day's HDF5 symlink. Runs are sequential, with
unchanged FFT/QC and accepted-window counts 47/19/26/42 (134 total). OS caches
were not purged. Segment totals exclude startup and input loading.

| Stage, default 67×67 grid | CPU | Metal | Speedup |
| --- | ---: | ---: | ---: |
| FFT/QC | 0.199 s | 0.205 s | — |
| Packing | 0.070 s | 0.067 s | — |
| Slant stack | 8.408 s | 0.522 s | 16.11× |
| Candidate search/fitting | 3.817 s | 3.798 s | — |
| Total segment processing | 12.536 s | 4.628 s | 2.71× |

Synthetic and real throughput differ because real windows have different
station counts and spectra, and sequential GPU submissions have a different
execution cadence from alternating CPU/GPU benchmarks. The real comparison
uses the same input and 16 threads on both paths; speedup is not attributed to
increasing CPU threads. Fitting and input loading now limit overall throughput.

For the wider 201×201 grid, a separate CPU/GPU pair processed the same day:

| Stage, 201×201 grid | CPU | Metal | Speedup |
| --- | ---: | ---: | ---: |
| Slant stack | 75.081 s | 3.678 s | 20.42× |
| Candidate search/fitting | 1.817 s | 1.790 s | — |
| Total segment processing | 77.229 s | 5.805 s | 13.30× |

This wider-grid GPU run used the installed executable. It has about nine
times as many cells as the default grid, yet its slant-stack time is below
the default-grid CPU time. Fitting times differ between grids because the
candidate sets differ; compare backends within each grid, not across grids.

## Numerical validation and scope

Maximum scaled error means `max(abs(GPU-CPU))/max(abs(CPU))` over the RTU
array, not a per-cell relative error near zero. Relative L2 error is the norm
of the difference divided by the CPU norm. Verification limits are 5e-4 for
both metrics, with a 1e-30 absolute floor on maximum error.

- Synthetic regression: 32 combinations of 2c/3c, random/coherent plane-wave
  inputs, centered/offset grids and 1/31/155/257 frequency bins. It checks
  finite values, inactive station capacity, pre-existing result accumulation,
  repeated buffer reuse, unchanged U output in horizontal mode and the
  fewer-than-two-stations guard. All passed; coherent peaks matched.
- Production-size benchmarks above: all R/T/U peak cells matched. Maximum
  scaled error was at most 2.27e-5, and relative L2 error at most 1.66e-5.
- Real default grid: all 134 windows checked against double CPU, maximum
  scaled error 9.56e-6, maximum relative L2 5.89e-6, zero R/T peak mismatches.
- Real 201×201 grid: all 134 windows checked, maximum scaled error 1.43e-5,
  maximum relative L2 8.45e-6, zero R/T peak mismatches.
- Default-grid event catalogs: same two R events; times, source positions,
  slownesses and flags match at output precision. Beam MAD changes in the last
  printed digit; beam maximum/MAD are checked at rtol=1e-4. The convergence
  residual uses absolute tolerance 1e-12. Bootstrap fields (columns 21–32)
  retain time-based seeding and are excluded, as in CPU comparisons.
- Wide-grid event catalogs: CPU and GPU detect the same one R event, passing
  the same event comparison. The count differs from the default grid on both
  backends. The default CPU catalog also passes the original strict comparison
  against the pre-Metal `51f60a7` run.
- Existing CTest suite: 3/3 passed with Metal compiled and with
  `DELTAP_ENABLE_METAL=OFF`. The CPU-only executable clearly rejects Metal.
  This is a macOS CPU-only build check, not an Ubuntu build measurement.
- Installed binary: real wide-grid GPU run and parent-launcher fixture passed.
  The launcher fixture runs from another working directory with an empty
  archive; the real-day verification separately exercises actual GPU stacks.

Float results are not bitwise identical. Near-tied peaks/thresholds can change
with other data; one day is not a full-archive scientific validation. Real 3c
catalogs were not tested; 3c kernel correctness and performance were tested.
Tiny problems can be slower on GPU because submission overhead dominates.
Finer/wider grids intentionally change the candidate search and can change
event counts even with identical CPU/GPU results at a given grid.

## Reproducing

Build/install as in the manual. Run GPU tests outside a sandbox that denies
Metal device access:

```bash
ctest --test-dir build-clang --output-on-failure
build-clang/test_metal_slant_stack
# Arguments: stations, half-width, CPU threads, repetitions, optional step.
build-clang/test_metal_slant_stack 650 33 16 5
build-clang/test_metal_slant_stack 650 66 16 5 .0025
build-clang/test_metal_slant_stack 650 100 16 5 .0025
```

For a real day, first create a bounded `HINET_ROOT/YYYY/MMDD/` input tree.
The local `build-gcc/one-day` directory is only such a symlink tree; its name
has no relation to the compiler used. Use unique result IDs on every run:

```bash
AUTOFOCUSING_BACKEND=metal AUTOFOCUSING_PROFILE=1 OMP_NUM_THREADS=16 \
  build-clang/src/cal_ccf_clang 2004 tilt_horizontal metal-check \
  build-gcc/one-day ../moment_loc_76_24 build-clang/metal-verification/results \
  horizontal > build-clang/metal-verification/metal-check.log 2>&1
```

Repeat with `AUTOFOCUSING_BACKEND=cpu` and a different result ID/log. Add
`AUTOFOCUSING_SLOWNESS_STEP=.0025 AUTOFOCUSING_SLOWNESS_MAX=.25` to both for
the wider grid. For every-window accuracy checks, separately run Metal with
`AUTOFOCUSING_VERIFY_METAL=1`; those timings include CPU verification and must
not be used as GPU speed measurements.

```bash
python3 Scripts/compare_profiles.py cpu.log metal.log
python3 tests/compare_event_results.py cpu.dat metal.dat --metal
```

The `--metal` comparison only relaxes beam maximum/MAD; CPU comparison defaults
remain unchanged. Machine-readable timings and validation summaries are in
[the measurement artifact](benchmarks/metal-slant-stack-20260908.json).
