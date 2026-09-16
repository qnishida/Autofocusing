# Metal power-evaluation optimization

[Japanese](metal-power-optimization_jp.md)

## Scope and usage

Based on `metal` (`f745119`), which already includes the shared I/O improvements,
`perf/metal-power` implements experimental objective-function evaluation for bootstrap
and the initial fitting grid. This is separate from the slant-stack `(px, py)` grid.

```bash
export AUTOFOCUSING_BACKEND=metal
export AUTOFOCUSING_METAL_POWER=bootstrap  # off / bootstrap / grid / all
export OMP_NUM_THREADS=16
```

The default is `off`. `bootstrap` evaluates 100 samples and the final corrected power
on Metal; `grid` evaluates 35 distance candidates and 40 curvature candidates; `all`
enables both. The CPU backend always uses the existing double-precision path.
Because the added paths have not yet been validated with real three-component data,
they are restricted to horizontal mode. Three-component mode uses the existing CPU
objective function, while retaining the existing three-component Metal slant stack.

The final fitting objective, gradient, Hessian, random-number generation, sample count,
surrounding statistical calculations, FFT, and I/O remain unchanged.
Normal bootstrap runs retain time-based seeds, so exact equality of bootstrap values
across separate runs is not required. Comparisons use fixed seeds in validation drivers.

## Implementation

`metal_power_batch` is an internal API accepting parameter sets, spectra, station
coordinates, and weights. GPU arithmetic uses FP32; results are returned to the CPU
as double. Returning double does not remove GPU numerical error.

- Spectra are prepared once per batch call, and GPU buffers are reused.
- At most 32 candidates are evaluated at a time, limiting temporary GPU memory for expanded grids.
- Bootstrap phases are shared across samples, generated for each group of 32 candidates.
- Stations are summed for each frequency and time window. Power and self-term corrections are aggregated separately using hierarchical reductions.
  Metal fast math is disabled.
- For the grid, candidates within `2e-3 * abs(max) + 1e-30` of the GPU maximum are reevaluated in parallel using CPU double precision.
  Final selection uses the original order and strict comparisons.
  This margin is an engineering criterion with headroom over the largest synthetic-test error, not a general error bound.
- `AUTOFOCUSING_PROFILE=1` emits `#POWER_PROFILE stage=grid|bootstrap total_s=...`.
  Timings include preparation, transfers, synchronization, and reevaluation. Bootstrap time also includes the existing Hessian calculation.

## Accuracy validation

- All four existing CTests passed in both Metal-enabled and Metal-disabled builds.
- The existing Metal slant-stack regression tests passed on a real GPU.
- All 216 synthetic cases for the new kernel passed: 2/31/650 stations, 1/3 windows, 5/90/175 degrees,
  random/coherent waveforms, amplitudes `1e-12/1/1e6`, masks and nonuniform weights, with/without correction.
  The maximum error relative to uncorrected double-precision power was `1.73e-4`, against a threshold of `5e-4`.
- Indexing across the 32-candidate batch boundary was checked with 65 distinct parameter and weight sets.
  Zero signals, tied candidates, and disabling the path for CPU/three-component mode were also checked.
- The off, bootstrap, grid, and all modes were compared against 29 events from the existing Metal version for 2004-01-01 through 05.
  The final two days were independent validation days, unused in tuning the candidate-reevaluation margin.
  The first three days contain 20 events and the final two contain nine.
  These must not be confused with the 30 events in the earlier CPU investigation.
- All initial-candidate logs, event identities and counts, and non-bootstrap columns matched.
  The absolute tolerance for convergence residuals is `1e-12`. Bootstrap-related columns 21–33
  passed with `rtol=1e-3, atol=1e-30`; the maximum relative difference was approximately `2.87e-5` (0.0029%).
  Existing NaNs for missing U in horizontal mode are preserved; new non-finite values are rejected.
- The existing `compare_event_results.py` excludes bootstrap columns and is therefore not used for this check.
  `run_metal_power.py` also compares those columns.

## Performance and adoption results (2026-09-12)

Apple M4 Max, 16 CPU threads. The 20 events from 2004-01-01 through 03 were measured
five times per path in alternating order. Each trial ran in a separate process and
included initialization. OS and driver caches were not cleared. Medians are below.

| Path | Total time | Bootstrap stage | Grid stage |
| --- | ---: | ---: | ---: |
| off | 41.895 s | 22.912 s | 1.076 s |
| bootstrap | 20.641 s | 1.670 s | 1.068 s |
| grid | 41.484 s | 22.743 s | 0.859 s |
| all | 20.417 s | 1.668 s | 0.845 s |

Enabling both gives a **2.05× overall speedup**. Individual stages improve by 13.72×
for bootstrap and 1.25× for the grid.
Both candidates meet the adoption criteria: at least a 10% reduction in the target stage
and no more than a 5% increase in total time.
Grid selections and non-bootstrap output for real events match the original Metal version.
The maximum relative difference in bootstrap-related values remained within approximately
0.0029% across all repetitions.

- Maximum process RSS for off: 2.970 GB.
- Maximum process RSS for bootstrap: 3.158 GB.
- Maximum process RSS for grid: 3.194 GB.
- Maximum process RSS for all: 3.193 GB.

These are whole-process peaks, not GPU memory alone.
See the [repeated-measurement JSON](benchmarks/metal-power-events-20260912.json) for speed,
ranges, CPU time, and memory details.

For grid-only evaluation with synthetic data, the standard 75 candidates yielded
0.89–0.97×, making the GPU slightly slower. With 750 candidates the speedup was
1.61–1.91×, and with 7500 candidates it was 2.38–2.45×.
Every case selected the same candidate as the CPU. The GPU is not necessarily faster
for every small input. This expanded evaluation uses synthetic data; changes to real-data
search settings require separate validation. See the [expanded-grid JSON](benchmarks/metal-power-grid-20260912.json).

Both candidates are selected for adoption and integration into `metal`. The default remains `off`.
Leave `main`, the prefetch TODO, and installed binaries unchanged; do not push.

## Validation records

- [Kernel tests](benchmarks/metal-power-kernel-20260912.json)
- [Initial five-day comparison](benchmarks/metal-power-qualification-20260912.json)
- [Comparison after parallelizing double-precision reevaluation](benchmarks/metal-power-refined-qualification-20260912.json)
- [29 events from the original Metal version](benchmarks/metal-power-reference-events-20260912.dat)

Initial comparison timings may include interference from concurrent builds, so they
are not used for the performance adoption decision.
After parallelizing double-precision reevaluation without numerical changes, the
independent validation days were compared again as well.

## Reproduction

```bash
cmake -S . -B build-power -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++" \
  -DCMAKE_OBJCXX_COMPILER="$(xcrun -f clang++)" \
  -DDELTAP_ENABLE_METAL=ON
cmake --build build-power -j4
OMP_NUM_THREADS=4 ctest --test-dir build-power --output-on-failure
OMP_NUM_THREADS=4 build-power/test_metal_power
OMP_NUM_THREADS=4 build-power/test_metal_slant_stack
build-power/test_metal_power --benchmark-grid
python3 tests/build_io_probe.py build-power --reference f745119
python3 tests/build_io_probe.py build-power
python3 tests/run_metal_power.py \
  build-power/io-fixed-current/cal_ccf_io_probe "$HINET_ROOT" "$CMT_CATALOG" \
  build-power/qualification-new --days 5 --repeats 1 \
  --legacy build-power/io-fixed-reference/cal_ccf_io_probe
python3 tests/run_metal_power.py \
  build-power/io-fixed-current/cal_ccf_io_probe "$HINET_ROOT" "$CMT_CATALOG" \
  build-power/benchmark-new --days 3 --repeats 5
```

Real-GPU tests require GPU access. Specify a new output directory.
Fixed FFTW plans and seeds apply only to test drivers.
Do not run measurements concurrently with other builds or performance tests.
`--benchmark-grid` measures 1×, 10×, and 100× candidate counts five times each with
650 stations and three windows. It separately tests refinement within the original range
and expansion of the range (distance 1–179 degrees, curvature −2e-5 to 4e-5, excluding
the upper endpoint), comparing CPU/GPU candidate selection in the two-stage search.
Normal operational search ranges and increments remain unchanged.
