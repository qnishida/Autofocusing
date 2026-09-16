# Shared CPU I/O optimization

Detailed benchmark logs and run results are retained locally and are not
distributed in this repository. References to detailed records below describe
those local artifacts.

[Japanese](cpu-io-optimization_jp.md)

## Implementation (2026-09-12)

`perf/cpu-io` is based on the CPU-optimized `main` (`51f60a7`).
Filter parallelization is in `ab175dc`; shared HDF5 metadata is in `46dd8e0`.

- HDF5 reads remain single-threaded. After reading, OpenMP processes independent `hp_filt` operations for each station and component.
- `load_h5` combines daily initialization and loading, sharing the file handle, station enumeration, and component metadata.
- Station-name lookup is indexed, and the station array is reserved in advance and constructed in place.
- The existing `init_station` / `read_h5` APIs are preserved, as are unit conversion, filter coefficients and arithmetic precision, station order, selection radius, and component QC.
- Missing or non-finite coordinate attributes now raise explicit exceptions instead of allowing undefined coordinate use. File handles are released on exceptions as well.

This work excludes next-day prefetching, concurrent HDF5 calls, GPU processing, and changes to cache capacity.
`OMP_NUM_THREADS` also controls filtering. `AUTOFOCUSING_PROFILE=1` emits
`#LOAD_PROFILE init_s=... read_decode_copy_s=... filter_s=... total_s=...`.
`read_decode_copy_s` includes HDF5 reading, decompression, conversion, copying, and remaining metadata processing; it is not pure SSD wait time.
The previous `#Read data: init_station/read_h5` lines are consolidated into `#Read data: load_h5`.

## Validation results

All four CTests passed in both the CPU and Metal-enabled builds.
`tests/io_reference.h` freezes the loader from `51f60a7` for comparison.
For synthetic horizontal two-component and three-component data, runs with 1, 4, and 16 threads
matched in station order, coordinates, accepted station counts, waveform metadata, and every
stored waveform sample (bitwise waveform equality).
Checks also covered the legacy APIs, repeated loading, missing components, invalid units and
sampling, non-orthogonal axes, missing files, no stations within the radius, missing coordinates,
and handle cleanup on exceptions.

The Metal validation branch `test/metal-cpu-io` reapplies the old Metal commit
`27fccd3` on top of the shared I/O changes. All 32 real-GPU cases passed on Apple M4 Max,
with zero peak mismatches and a maximum scaled error of approximately `5.05e-6`.
This checks for regressions in existing Metal calculations; it does not measure I/O speedup.

After connecting the SSD, full-waveform comparisons also passed for 2004-01-01 through 03,
2014-01-01, and 2024-01-01. Accepted station counts were 653, 653, 652, 725, and 725,
respectively. Each date was checked with 1, 4, and 16 threads, including repeated loading
and the legacy APIs. See the real-waveform validation record.

Medians for cached input from 2004-01-01, with five runs per implementation in alternating order:

| Threads | Before | After | Speedup |
| --- | ---: | ---: | ---: |
| 1 | 2.119 s | 2.080 s | 1.02× |
| 4 | 2.124 s | 1.298 s | 1.64× |
| 8 | 2.117 s | 1.144 s | 1.85× |
| 16 | 2.121 s | 1.080 s | 1.96× |

Input time fell by approximately 49% with 16 threads. Actual disk-read counters were
0 bytes in every trial, so these results do not measure USB bandwidth or uncached input.
The largest peak RSS was approximately 1.035 GB before and 1.011 GB after the change.
See the measurement JSON for medians, ranges, CPU time, and raw RSS data.

### End-to-end comparison

With FFTW plans and bootstrap seeds fixed in validation builds, 2004-01-01 through 03
was measured twice per implementation, reversing the execution order.
For each backend, all output bytes for 20 events (including bootstrap results) and all
433 accepted windows matched.

| Backend | Median total time before | After | Time reduction |
| --- | ---: | ---: | ---: |
| CPU, 16 threads | 71.225 s | 68.325 s | 4.1% |
| Metal, 16 CPU threads | 45.210 s | 42.325 s | 6.4% |

The integration criteria were met: at least a 10% improvement in multithreaded input time
and no more than a 5% increase in total time.
This comparison uses the version before the I/O optimization, not the version before
CPU slant-stack optimization. Normal builds retain FFTW_MEASURE and time-based bootstrap seeds.

- CPU measurements and comparisons, CPU event output
- Metal measurements and comparisons, Metal event output

## Reproduction

Use Homebrew LLVM for C/C++ and Apple Clang for Metal Objective-C++.
Explicitly selecting the latter avoids the Objective-C method-call linker errors observed
with LLVM 23 in this environment.

```bash
cmake -S . -B build-io -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++"
cmake --build build-io -j4
OMP_NUM_THREADS=4 ctest --test-dir build-io --output-on-failure
```

For the Metal validation branch, add the following configure options and use a separate build directory:

```bash
-DDELTAP_ENABLE_METAL=ON -DCMAKE_OBJCXX_COMPILER="$(xcrun -f clang++)"
# After building, run in an environment with GPU access.
OMP_NUM_THREADS=4 ./build-metal-io/test_metal_slant_stack
```

Example measurement after connecting the SSD (`test_io` exists in both branches):

```bash
python3 tests/benchmark_io.py build-io/test_io \
  /Volumes/Seismic_Data/hdf5/Hi-net_tilt/2004/0101/20040010000.h5 \
  docs/benchmarks/cpu-io-2004001.json
```

The script first compares all waveforms, then measures each implementation five times
with 1, 4, 8, and 16 threads in alternating order with a warm cache.
It saves wall-time medians, ranges, and ratios, CPU time, actual disk-read bytes, and peak
process RSS to JSON. macOS RSS is reported in bytes.
Each trial runs in a separate process with identical FFTW initialization conditions.
This is not a cold-cache or USB-bandwidth measurement.

## Reproducing event output

Configure with `CMAKE_EXPORT_COMPILE_COMMANDS=ON` and build normally, then fix the FFTW
plans and bootstrap seeds only in the test driver.
Use `51f60a7` as the reference commit for CPU and `27fccd3` for Metal.
Use a build directory for the same backend.

```bash
python3 tests/build_io_probe.py build-io --reference 51f60a7
python3 tests/build_io_probe.py build-io
python3 tests/run_io_events.py \
  build-io/io-fixed-reference/cal_ccf_io_probe \
  build-io/io-fixed-current/cal_ccf_io_probe \
  /Volumes/Seismic_Data/hdf5/Hi-net_tilt \
  ../moment_loc_76_24 build-io/io-events-cpu
```

For Metal, replace the build paths and add `--backend metal` to the run command.
Specify an output directory that does not yet exist. The script creates symlinked input
for only three days and runs each version twice, in before/after and after/before order.
It checks for 20 events, 12 segments, matching accepted-window counts and all output bytes
including bootstrap results, and no more than a 5% increase in total time.
Run logs, event files, and `report.json` remain in the specified output directory.
Save important aggregates in `docs/benchmarks/`.

## TODO: next-day prefetching (deferred, 2026-09-12)

At the user's direction, further prefetch optimization is deferred because current I/O
performance is sufficient. Cached input and preprocessing currently take approximately
1.08 seconds per day with 16 threads. Reassess with new measurements if NAS use, larger
datasets, or faster computation makes input waiting a significant share of total time.

The current `prefetchFile()` requests OS prefetching on Ubuntu/Linux through
`posix_fadvise(..., POSIX_FADV_WILLNEED)`.
That call is compiled out in the current macOS build, which only opens and closes the
next day's file. It was not removed by the CPU/I/O optimization.
macOS also provides `fcntl(..., F_RDADVISE, ...)`, but the current code does not implement it.
Neither OS hint manages completion of next-day waveform decompression and filtering.

- [ ] When work resumes, consider a mechanism shared by CPU and Metal that overlaps the current day's analysis with the next day's loading and preprocessing.
- [ ] With sufficient RAM in the current environment, first consider keeping one next-day dataset in RAM.
  Allocate capacity according to actual data rather than fixing it at 1 GB, and manage both retained dataset count and memory limits.
  Current compressed files are approximately 0.5 GB/day; expanded horizontal waveforms are approximately 0.9 GB/day, plus working memory.
- [ ] Separate daily shared state such as station-center coordinates, and avoid concurrent calls to non-thread-safe HDF5 functions.
  Limit and measure prefetch CPU usage to assess contention with the current day's analysis.
- [ ] Across multiple consecutive days, compare waveform and event-output equality, startup and steady-state timing, and peak RAM.
  Distinguish cache conditions when measuring NAS performance.

Temporary storage on a local SSD is an alternative when RAM is constrained or when
avoiding another transfer for later reanalysis. It has not been judged faster than RAM
prefetching under the current conditions with sufficient RAM.
This TODO records considerations only; no prefetch code has been changed.

## Branch integration

Fast-forward `main` to the shared I/O commits that passed real-data validation, then
reapply the existing Metal development commits on that head.
Align `perf/cpu-io` with `main` and `test/metal-cpu-io` with `metal`.
Preserve the old `metal` (`27fccd3`) as `backup/metal-before-cpu-io-20260912`.
Do not push or install binaries.

```text
51f60a7  CPU slant-stack optimization
  └─ Shared CPU I/O optimization and validation  (main, perf/cpu-io)
       └─ Metal backend                         (metal, test/metal-cpu-io)
```

Investigation materials are saved at `dd8bc8c` on the original
`investigate/metal-float-candidates` branch.
Code used only for GPU/FP32 investigation was not incorporated into the CPU branch.
