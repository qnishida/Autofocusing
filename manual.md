# Autofocusing Manual

This manual describes the practical steps needed to build and run Autofocusing
with external seismic array waveform data in HDF5 format and a Global
CMT-derived earthquake catalog.

## 1. Build and Install

Run CMake from the repository root, not from `src/`.

```bash
mkdir -p build
cmake -S . -B build
cmake --build build -j
cmake --install build
```

The install step places the main executable under:

```text
bin/cal_ccf_gcc
```

The executable suffix depends on the compiler selected by CMake. For example,
GCC produces `cal_ccf_gcc`.

## 2. Required Input Data

Autofocusing requires two external inputs:

- `HINET_ROOT`: root directory of the Hi-net HDF5 waveform archive.
- `CMT_CATALOG`: whitespace-separated earthquake catalog generated from Global
  CMT data.

The repository also includes `vel_Nishida2008.dat`, a station-correction
velocity model used at runtime. Keep this file in the repository root when
running the executables from this directory.

## 3. Hi-net HDF5 Waveform Archive

The analysis uses a Hi-net waveform archive converted from WIN/WIN32 seismic
data downloaded from NIED. The conversion used:

> Nishida, K. (2026). win32conv: A modern C implementation of WIN/WIN32
> seismic data conversion tools (Version 1.0.0) [Software].
> https://github.com/qnishida/win32conv.
> https://doi.org/10.5281/zenodo.19879657

See `NOTICE` for data provenance and citation notes.

The workflow is not limited to Hi-net in principle. Other seismic array data,
such as USArray data, can be used if they are converted to the HDF5 layout
expected by `cal_ccf_gcc`. This repository does not currently include a
maintained conversion script for those other data sets.

`cal_ccf_gcc` expects the waveform archive to be arranged by year and day:

```text
HINET_ROOT/
└── 2004/
    ├── 0101/
    │   └── one_hinet_file.h5
    ├── 0102/
    │   └── one_hinet_file.h5
    └── ...
```

For each day, the program looks under:

```text
HINET_ROOT/YYYY/MMDD/
```

Missing day directories are skipped. If a day directory exists, it is processed
only when it contains exactly one regular HDF5 file. Empty day directories,
directories containing multiple files, and missing files are skipped.

Malformed or unreadable HDF5 files may produce HDF5 errors instead of being
treated as ordinary missing data.

## 4. Global CMT Catalog Input

`CMT_CATALOG` must use the format expected by the C++ code:

```text
yy jday hour minute second longitude latitude depth moment
```

Raw Global CMT files and converted local catalogs are not included in this
repository. Use the helper and notes under `Scripts/GlobalCMT/` to regenerate
the local catalog input from Global CMT data, and cite Global CMT following
their guidance.

## 5. Running Autofocusing

For a new workspace, copy `Scripts/run.sh` to the parent of `repo/` as an
editable regular file, and copy `Scripts/local_config.example.sh` to the
parent's `local_config.sh`. Use `cp -n` for both to preserve existing files.
The local configuration is outside this repository. Launch with
`../run.sh` from the repository or an absolute path from elsewhere.

The script resolves paths relative to its workspace and changes into `repo/`
for runtime model files. It discovers installed `cal_ccf_gcc` or
`cal_ccf_clang`; `AUTOFOCUSING_BIN` can override the executable.

```bash
export HINET_ROOT="/Volumes/Seismic_Data/hdf5/Hi-net_tilt"
export CMT_CATALOG="moment_loc_76_24"
export COMPONENT_MODE="horizontal"
export PARAM_ID="tilt_horizontal"
export START_YEAR="2004"
```

These settings belong in the parent's `local_config.sh`, which is sourced by
the launcher. Relative input/output paths refer to the parent directory.
`RESULTS_ROOT` defaults to its `results/` directory. `OMP_NUM_THREADS` optionally
limits parallelism. The launcher does not build the program automatically.

### Horizontal input and processing

`horizontal` selects a complete `LE/LN` pair, falling back to `E/N`. U is not
read or required. Waveforms must be 2 Hz integer samples with `unit="nm/s"`;
`sample * sensitivity * 1e-9` converts to m/s. This mode accepts velocity
already converted from tilt, and does not perform that conversion again.

The two sensor axes must be orthogonal, with E azimuth equal to N azimuth
plus 90 degrees (within 0.01 degrees). The existing N azimuth rotation converts
sensor coordinates to geographic E/N. Already rotated waveforms must carry
geographic azimuths E=90 and N=0 to avoid a second rotation.

Each component must provide the existing sensitivity, azimuth, sampling,
length, and time metadata. Invalid lengths, rates, units, or incomplete pairs
are rejected. The existing requirement for more than half a day's samples
is retained. The `gap` attribute does not describe individual missing intervals;
this change retains the existing waveform/QC handling of gaps.

Horizontal amplitude checks use E and N only. Stability checks use their mean
band power instead of U. Existing thresholds, frequency band, window length,
and array selection are retained. Candidate searches use R and T only; no
vertical candidates or vertical-derived horizontal seeds are used.

### Direct invocation and compatibility

```text
cal_ccf_<compiler> YYYY <param-id> <git-version> <hinet-root> <cmt-catalog> [output-root] [3c|horizontal]
```

Omitting the optional arguments retains `output/` and `3c`. The parent launcher
uses `horizontal` and a separate result directory by default. `cal_ccf_eq` has
not been extended with horizontal mode.

The start year defaults to 2004 in the launcher. The existing scan stops at
2024-12-31 or its original 366 × 20.75-day iteration bound, whichever comes
first. This change does not add support for 2025 or a date-range interface.

## 6. Output

The launcher writes `results/<param-id>/<git-version>/` under the workspace.
The filename retains `YYYY_<segment length>_<min frequency>-<max frequency>.dat`.
A direct invocation defaults to `output/` if no output root is supplied.

The 38-column format remains unchanged. Horizontal events use component
numbers R=0 and T=1. Columns 34–35 contain R/T matrix diagonal values;
columns 36–38 (U power and real/imaginary R–U cross-power) are `nan` in
horizontal mode. These values indicate missing observations, not zero power.

Each invocation truncates its output file. Use distinct parameter IDs for
separate experiments. Git identifiers include a `-dirty` suffix for tracked
uncommitted changes; they do not uniquely identify each uncommitted edit.

## 7. Citation

If you use Autofocusing, please cite the Zenodo release:

```text
https://doi.org/10.5281/zenodo.20742021
```

See `CITATION.cff` for machine-readable citation metadata.


## 8. Troubleshooting

If `hdf5.h` is missing during compilation, install the HDF5 development package
or point CMake to the HDF5 installation prefix:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/hdf5/prefix
cmake --build build -j
```

If CMake warns that no top-level `project()` or `cmake_minimum_required()` is
present, it was probably run from `src/`. Reconfigure from the repository root:

```bash
rm -rf build src/CMakeCache.txt src/CMakeFiles
cmake -S . -B build
cmake --build build -j
```

If the parent launcher cannot find an installed executable, run the install step:

```bash
cmake --install build
```

If no output is produced for expected days, check that each day directory exists
and contains exactly one regular HDF5 file under `HINET_ROOT/YYYY/MMDD/`.

## 9. macOS with Homebrew libraries

Use Clang with Homebrew's C++ libraries. GCC can compile this project on macOS
but mixing its libstdc++ with Homebrew Boost's libc++ can break filesystem
operations at runtime. A consistent Homebrew LLVM build is:

```bash
cmake -S . -B build-clang \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++" \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-clang -j 4
ctest --test-dir build-clang --output-on-failure
cmake --install build-clang
```

This requires the dependencies listed in README plus LLVM/OpenMP. Current
Boost/Eigen require C++14. The C language is enabled for HDF5 discovery;
Boost.System is header-only in modern Boost and is not requested separately.

The synthetic regression check can additionally compare the current 3c loader
against a specified Git revision:

```bash
python3 tests/compare_legacy.py build-clang 149995f
```

## 10. CPU performance profiling

Slant stacking uses double precision with OpenMP. To enable per-segment timings,
set `AUTOFOCUSING_PROFILE=1`; `OMP_NUM_THREADS` controls worker count. These
variables can be exported before invoking the launcher or in the parent local
configuration. Timings cover FFT/QC, packing, slant stack, fitting and total
segment processing; existing read logs report initialization/loading separately.

Compare equivalent runs using `python3 Scripts/compare_profiles.py old.log new.log`.
Use a bounded input tree and a distinct parameter/version ID when benchmarking
so an experiment does not launch the full archive or replace previous results.
The CPU optimization, numerical checks, measurements and Metal follow-up
boundary are documented in [the performance report](docs/cpu-slant-stack-performance.md).

### 共通 I/O の処理と検証

日単位の HDF5 読み込み後に波形フィルタを CPU 並列処理します。
スレッド数は `OMP_NUM_THREADS`、内訳の計測ログは `AUTOFOCUSING_PROFILE=1` で設定できます。
実装範囲と検証状況は [共通 CPU I/O 最適化](docs/cpu-io-optimization.md) を参照してください。

## 11. Metal GPU slant stacking

On macOS, CMake builds the optional Metal backend by default. Use the Clang
configuration above; Apple Metal/Foundation frameworks and an accessible Metal
GPU are required. Shaders compile from embedded source at startup, so the
standalone `metal` compiler/full Xcode installation is unnecessary. On other
platforms the Metal backend is disabled by default. To build only the portable
CPU path, configure with `-DDELTAP_ENABLE_METAL=OFF`.

The runtime default is still `AUTOFOCUSING_BACKEND=cpu`. To use the GPU, add
these exports to the parent's `local_config.sh`, or export them before launching
if that file does not override them:

```bash
export AUTOFOCUSING_BACKEND=metal
export OMP_NUM_THREADS=16
export PARAM_ID=tilt_horizontal_metal
```

Both `horizontal` and `3c` are supported. Only slant stacking and its power
reduction use float on the GPU. Input loading, rotation, FFT, QC and subsequent
estimation retain the existing double CPU path. GPU powers are promoted and
accumulated into the existing double result array. Startup logs identify the
backend/device; an explicit Metal request fails if no GPU is accessible or the
backend was not built. It does not silently switch to CPU. Sandboxed execution
may require GPU access permission even when the same executable works normally
in a terminal. This backend applies to `cal_ccf`, not `cal_ccf_eq`.

Grid settings apply equally to CPU and Metal, without recompilation:

| Variable | Default | Meaning |
| --- | --- | --- |
| `AUTOFOCUSING_SLOWNESS_STEP` | `0.005` | Grid spacing, s/km |
| `AUTOFOCUSING_SLOWNESS_MAX` | `0.165` | Maximum absolute px/py, s/km |

For example, spacing `0.0025` and maximum `0.25` produce a 201×201 grid
instead of the default 67×67. The upper limit is rounded down to a grid multiple;
the actual grid is logged. Positive finite settings and a half-width from 1 to
1024 are required. Finer/wider grids change the search itself and may change
detected events. Use distinct `PARAM_ID` values for backend/grid experiments;
the output format and file-overwrite behavior remain as described in section 6.

`AUTOFOCUSING_PROFILE=1` measures normal processing, including float packing,
submission, GPU completion and result accumulation in `stack_s`. Startup shader
compilation is outside that stage. For numerical diagnosis only, export
`AUTOFOCUSING_VERIFY_METAL=1`: every window also runs the double CPU kernel,
logs maximum scaled/L2 error and peak mismatches, and fails on error exceeding
5e-4 (with a 1e-30 absolute floor for maximum error). Unset this variable for
speed measurements; setting it to `0` still enables the check.

Validation commands, measured speedups and numerical limits are recorded in
[the Metal report](docs/metal-slant-stack-performance.md).

### Optional Metal fitting objectives

For validated horizontal inputs, `AUTOFOCUSING_METAL_POWER` selects additional
GPU power evaluations: `off` (default), `bootstrap`, `grid`, or `all`.
Set `AUTOFOCUSING_BACKEND=metal` as well. CPU and three-component runs retain
the existing CPU objectives. Final gradient fitting and Hessians are unchanged.
See [Metal power evaluation](docs/metal-power-optimization.md) for numerical
limits, measurements and reproduction commands. Use a distinct result ID when
enabling these options; they do not change the existing overwrite policy.


## 12. NVIDIA CUDA GPU backend

CUDA supports the same `cal_ccf` stages as Metal: horizontal/3c slant stacking,
and optional horizontal Bootstrap/initial-grid objectives. CPU remains the
runtime default. FFT, I/O, final fitting and 3c fitting objectives remain on CPU.
GPU arithmetic is FP32; returned doubles do not restore the lost precision.

The CUDA build requires CMake 3.18 or newer, a CUDA Toolkit compatible with the
host compiler/GPU, and the NVIDIA driver. For RTX PRO 2000 Blackwell, use the
installed CUDA 13.2 Toolkit with compute capability 12.0:

```bash
cmake -S . -B build-cuda -DDELTAP_ENABLE_METAL=OFF \
  -DDELTAP_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120 \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-cuda -j4
ctest --test-dir build-cuda --output-on-failure
build-cuda/test_cuda_slant_stack
build-cuda/test_cuda_power
```

Use the architecture of your device for other NVIDIA GPUs. CPU-only builds
need no CUDA Toolkit: set `DELTAP_ENABLE_CUDA=OFF`. Metal retains its existing
Apple-only build setting. CUDA kernels are compiled into the executable;
there is no shader file or runtime compiler requirement.

Set these in the environment or the launcher's parent `local_config.sh`:

```bash
export AUTOFOCUSING_BIN="/absolute/path/to/build-cuda/src/cal_ccf_gcc"
export AUTOFOCUSING_BACKEND=cuda
export AUTOFOCUSING_GPU_POWER=bootstrap
export OMP_NUM_THREADS=16
export PARAM_ID=tilt_horizontal_cuda
```

`AUTOFOCUSING_BACKEND` accepts `cpu`, `metal`, or `cuda`. An unavailable or
unbuilt explicit GPU backend fails with a diagnostic; it never silently falls
back to CPU. CUDA uses the current runtime device (normally device 0);
`CUDA_VISIBLE_DEVICES` can restrict device visibility. GPU tests need access to
the host GPU device nodes, including when running inside a container/sandbox.

`AUTOFOCUSING_GPU_POWER` accepts `off` (default), `bootstrap`, `grid`, or `all`.
For the measured RTX PRO 2000 workload, `bootstrap` is recommended: CUDA grid
evaluation passes accuracy checks but is slower than CPU grid evaluation.
The setting applies to the selected GPU only in horizontal mode. The legacy
`AUTOFOCUSING_METAL_POWER` is still accepted for Metal when the new setting is
absent; the new setting takes precedence, including `off`. The legacy setting
has no effect on CUDA. CPU and 3c fitting always retain CPU objectives.

For accuracy checks, `AUTOFOCUSING_VERIFY_CUDA=1` evaluates each slant-stack
window on both CPU and CUDA, reporting scaled maximum error, relative L2 error
and peak mismatches. Existing limits of `5e-4` apply; peak mismatches are
reported for inspection. `AUTOFOCUSING_PROFILE=1` includes packing, transfers
and synchronization in stage timings. Separately, `AUTOFOCUSING_CUDA_PROFILE=1`
reports kernel-only CUDA-event time and reusable device-buffer capacities.
Do not enable either verification or CUDA-event profiling during throughput
measurements. The latter adds event instrumentation and log overhead.

See [CUDA implementation and validation](docs/cuda-performance.md) for tests,
measurements and the remaining qualification limits.

## 13. Workspace and Git management

The workspace is `Autofocusing/`; its only active Git checkout is `repo/`.
The parent contains `run.sh`, `local_config.sh`, the catalog (or a link to
its original), and `results/`. It has no `.git` or duplicate source tree.
Open `repo/` as the development workspace and run Git/build commands there.
Build products stay inside ignored directories such as `repo/build-cuda/`.
When relocating a checkout, configure a fresh build rather than reusing a
CMake cache containing the old absolute paths.

The tracked scripts are templates. The parent launcher can be edited for an
analysis; it is not replaced by `git pull`. Keep routine environment choices
in `local_config.sh`. For example, the following executable path is relative
to the parent workspace, so it survives moving the whole workspace:

```bash
export AUTOFOCUSING_BIN="repo/build-cuda/src/cal_ccf_gcc"
export RESULTS_ROOT="results"
```

After updating source templates, compare them with the parent copies using
`diff -u Scripts/run.sh ../run.sh` and
`diff -u Scripts/local_config.example.sh ../local_config.sh`. A difference
is expected for local settings; apply only the changes needed for the analysis.
Do not copy templates over existing customized files automatically.

GitHub synchronizes commits in `repo/`, not the parent workspace or uncommitted
source edits. Directory relocation neither creates commits nor pushes them.
Check `git rev-parse --show-toplevel` and `git status` before Git operations.
Preserve unfinished changes before switching branches or updating. Fetch to
inspect upstream changes, and use fast-forward-only pulls when ready. A
divergence needs deliberate integration; do not reset or force-push to make
the checkouts appear synchronized. Branches `main`, `metal` and `cuda` remain
separate until explicitly integrated.

For an SSH remote, `git ls-remote --heads origin` checks read access without
changing the checkout. `Permission denied (publickey)` indicates an SSH
authentication problem, independent of the directory layout. Cached
`origin/*` refs are not proof of the current GitHub state. Fix access using
the intended account/key before fetching or pushing; do not change the remote
or generate replacement keys automatically.
