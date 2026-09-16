# gpu branch: Metal verification (2026-09-14)

Detailed benchmark logs and run results are retained locally and are not
distributed in this repository. References to detailed records below describe
those local artifacts.

Verified unmodified `origin/gpu` commit `dabb03727454421f4fbb0b42a61604539fe99a6d`
on Apple M4 Max, macOS 26.6.2, with 16 OpenMP threads. The working checkout
remained on `cuda`; no source changes, installation, commits or pushes were made.

## Build

Homebrew LLVM/Clang 23.1.1 for all languages failed to link:
`_objc_msgSendClass$new$_OBJC_CLASS_$_MTLCompileOptions` and
`_objc_msgSendClass$stringWithUTF8String:$_OBJC_CLASS_$_NSString` were unresolved.
Selecting Apple Clang for Objective-C++ resolved the failure without source edits.
This matches the existing Metal power reproduction configuration.

```bash
cmake -S . -B build-apple-metal \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++" \
  -DCMAKE_OBJCXX_COMPILER="$(xcrun --find clang++)" \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DDELTAP_ENABLE_METAL=ON -DDELTAP_ENABLE_CUDA=OFF
cmake --build build-apple-metal -j 4
ctest --test-dir build-apple-metal --output-on-failure
build-apple-metal/test_metal_slant_stack
build-apple-metal/test_metal_power
```

Metal requires execution outside the agent sandbox on this host. Inside it,
Metal reported no accessible GPU and the qualification helper failed under
`/usr/bin/time -l`. The launcher also requires Git metadata; the initial archive
copy lacked it. Restoring metadata and rerunning outside the sandbox passed all
7 CTest tests. Neither initial failure required a source change.

## Results

- Slant stack: 32 synthetic cases covering horizontal and 3c passed;
  maximum scaled error `5.04548961e-6`, maximum relative L2 error
  `3.46418047e-6`, and zero peak mismatches.
- Power: 216 synthetic cases passed, maximum scaled error
  `1.7298449269270586e-4` (threshold `5e-4`); batched/chunk-boundary checks passed.
- Real horizontal data, 2004-01-01 through 2004-01-03: 20 events for CPU/off
  and Metal off/bootstrap/grid/all. CPU versus Metal/off passed the existing
  GPU event comparison (Bootstrap excluded). All Metal modes preserved initial
  choices and accepted windows; non-Bootstrap fields matched Metal/off, allowing
  convergence-residual absolute error `1e-12`. Maximum Bootstrap relative
  difference was approximately `2.87e-5` (threshold `1e-3`).

The real-data driver fixes FFTW planning and Bootstrap seeds for comparison;
production sources and binaries were not modified.

```bash
python3 tests/build_io_probe.py build-apple-metal
python3 tests/run_metal_power.py \
  build-apple-metal/io-fixed-current/cal_ccf_io_probe \
  "$HINET_ROOT" "$CMT_CATALOG" event-check \
  --days 3 --repeats 1 --modes off bootstrap grid all \
  --compare-cpu --expected-events 20
```

Single-run wall times were CPU/off 64.872 s, Metal/off 34.749 s,
Metal/bootstrap 13.249 s, Metal/grid 34.434 s, and Metal/all 13.108 s.
These are functionality-check timings, not repeated performance qualification.
Real 3c data, five-day validation, and paired pre/post CPU-change performance
qualification were not performed; this does not complete all `main` merge gates.

Raw report,
CTest,
slant stack,
power, and
initial linker failure
are retained alongside real-event logs and outputs.
