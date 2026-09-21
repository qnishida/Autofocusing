# Autofocusing

Autofocusing is a C++ implementation of the auto-focusing workflow used to
estimate source locations from seismic array data in HDF5 format.

This repository contains the analysis code and small helper scripts. Large
waveform inputs, Global CMT raw catalog files and CMT catalogs converted 
for use with Autofocusing are not distributed under the Autofocusing 
software license.

For detailed usage, input layout, and troubleshooting notes, see [the manual](manual.md).
The [documentation guide](#documentation-guide) links current usage and the validation history.

## Method Background

For a compact Python demonstration of the underlying centroid single-force and
auto-focusing ideas, see `CSF_supplement`:

```text
https://github.com/qnishida/CSF_supplement
```

The auto-focusing method was introduced for a global centroid single-force
catalog of P-wave microseisms in:

> Nishida, K., & Takagi, R. (2022). A global centroid single force catalog of
> P-wave microseisms. Journal of Geophysical Research: Solid Earth, 127,
> e2021JB023484. https://doi.org/10.1029/2021JB023484

A preprint describing an improved application to S-wave microseisms is
available as:

> Kato, S., Nishida, K., & Takagi, R. (2026). Origin of S-wave microseisms.
> https://doi.org/10.22541/essoar.15002647/v1

## Requirements

The build uses CMake and C++14. Dependency discovery uses `find_package` where
available and `find_path`/`find_library` for FFTW3 and GeographicLib.

Required libraries:

- Boost
- Eigen
- FFTW3
- GeographicLib
- HDF5
- OpenMP

If a dependency is installed in a non-standard location, pass its prefix or
specific CMake cache variables at configure time, for example:

```bash
cmake ../CMakeLists.txt -DCMAKE_PREFIX_PATH=/path/to/prefix
```

For machine-specific optimization, configure with:

```bash
cmake ../CMakeLists.txt -DDELTAP_ENABLE_NATIVE_ARCH=ON
```

## Build

Use an out-of-source build:

```bash
mkdir -p build
cd build
cmake ../CMakeLists.txt
cmake --build . -j
cmake --install .
```

The install step places the main executable in `bin/`:

```text
bin/cal_ccf_gcc
```

depending on the compiler suffix selected by CMake.

## Input Data

`cal_ccf_gcc` requires two external inputs:

- `HINET_ROOT`: root directory of the Hi-net HDF5 waveform archive.
- `CMT_CATALOG`: a whitespace-separated earthquake catalog generated from
  Global CMT data.

The Hi-net HDF5 waveform archive used for the analysis was converted from
WIN/WIN32 seismic data downloaded from NIED. The conversion used:

> Nishida, K. (2026). win32conv: A modern C implementation of WIN/WIN32
> seismic data conversion tools (Version 1.0.0) [Software].
> https://github.com/qnishida/win32conv.
> https://doi.org/10.5281/zenodo.19879657

See `NOTICE` for data provenance and citation notes.

Although the published analysis used Hi-net data, the workflow can be adapted
to other seismic array data if they are converted to the HDF5 layout expected by
the code.

The repository includes `vel_Nishida2008.dat`, a station-correction velocity
model derived from Nishida et al. (2008, JGR, doi:10.1029/2007JB005395).
Keep this file in the repository root when running the executables from this
directory.

The catalog format expected by the C++ code is:

```text
yy jday hour minute second longitude latitude depth moment
```

Global CMT raw catalog files and converted catalogs are not included in this
repository. See `Scripts/GlobalCMT/README.md` for the conversion command and
citation notes.

## Run

The workspace separates public source, analysis definitions and generated results:

```text
Autofocusing/
├── repo/                         # This public Git repository
├── run.sh                        # Delegates to repo/Scripts/run.sh
├── local_config.sh               # Machine-local settings
├── analysis/                     # Separate, optionally private Git repository
│   └── primary-microseisms/       # config.sh and analysis notes
└── results/                      # Generated outputs, outside Git
```

From the source checkout, create missing directories and templates:

```bash
bash Scripts/setup_workspace.sh --experiment primary-microseisms --init-git
# Edit ../local_config.sh and ../analysis/primary-microseisms/config.sh first.
bash Scripts/run.sh --experiment primary-microseisms
```

Setup preserves existing files and initializes only `analysis/` when `--init-git`
is requested. Do not initialize Git in the parent workspace. Existing customized
parent launchers are preserved and need a deliberate migration to the new runner.
Use `--workspace PATH` to select a workspace other than the checkout's parent.
Bash, Git and Python 3.7+ are required for the scripts.

The template still uses the nominal **0.1–0.25 Hz** band and ±0.165 s/km px/py
range. The experiment name does not select a frequency band. To analyze
primary microseisms, set `AUTOFOCUSING_FREQ_MIN=0.05` and
`AUTOFOCUSING_FREQ_MAX=0.1` in the experiment's `config.sh`. Both endpoints round
down to FFT bins; requested and effective bands are recorded. See
[frequency settings](Scripts/README.md#analysis-frequency-band).
The launcher defaults to horizontal-only analysis. Input paths and inclusive
`START_DATE`/`END_DATE` values must be configured in the workspace.
`--experiment NAME` is required: omitting it prints usage and an example and
exits without starting an analysis. See
[Scripts/README.md](Scripts/README.md) for setting precedence,
run records, setup options and migration, and [manual.md](manual.md) for the
scientific input requirements and build details.

Direct execution remains available from the repository root:

```bash
./bin/cal_ccf_clang 2004 tilt_horizontal example-run \
  /path/to/hdf5 /path/to/moment_loc_76_24 ../results horizontal 2004-01-01 2004-01-07
```

Here `example-run` is an output-directory label; use a new label for each direct run
and record the source release or commit separately. This example also works with
a source ZIP without Git history.

The compiler determines the executable suffix (`gcc` or `clang`). The optional
trailing arguments are `[output-root] [3c|horizontal] [START_DATE END_DATE]`.
The date pair requires both preceding options, and `YYYY` must match the start
date year. Omitting all optional arguments preserves three-component mode and
the `output/` destination.

The launcher scans the configured date interval, including both endpoints.
Missing/empty/multiple-file daily directories are skipped until the next day
within that interval. Input directories use `HINET_ROOT/YYYY/MMDD/` and must
contain exactly one regular file, interpreted as HDF5. If no input is loaded,
the log says so and the result is empty. Explicit dates remove the historical
2024 cutoff; year-only direct calls retain their old bounds. Rebuild and install
the executable to use the new date arguments. See
[date settings](Scripts/README.md#analysis-dates-and-missing-days).

## Output

The parent launcher writes to:

```text
results/<experiment>/<UTC-timestamp>/YYYY_2048_0.099609-0.250000.dat
```

Here `YYYY` is the start year, not the year of each individual event.
The experiment name is explicitly selected with `--experiment`.
Horizontal mode preserves the 38-column format, emits R=0 and T=1 events, and
writes `nan` for U-related matrix values (columns 36–38). Each launcher execution
allocates a new run directory and records settings,
source/analysis commits and logs alongside the results. Direct executable calls
still overwrite the `.dat` file when their destination is reused.

Generated build directories, binaries, and outputs are not tracked.
[Verification notes](docs/horizontal-verification.md) describe the input
checks, tests, and remaining limitations.

## Citation

If you use this software, please cite the released version you used. For example,
for version 1.0.0:

> Nishida, K. (2026). Autofocusing (Version 1.0.0) [Software].
> https://doi.org/10.5281/zenodo.20742021

The DOI above refers to all versions. Citation metadata in `CITATION.cff` are
being prepared for version 1.1.0, which has not yet been released. Set
`date-released` to the actual publication date when preparing the release.

## License and Attribution

Autofocusing is distributed under the GNU General Public License, version 2 or, at
your option, any later version. See `LICENSE`.

Third-party code, external library licenses, and data provenance notes are
summarized in `NOTICE`. Global CMT raw catalog files and CMT catalogs 
converted for use with Autofocusing are not distributed under the Autofocusing 
software license. Use `Scripts/GlobalCMT/` to regenerate the local catalog input 
from Global CMT data and cite Global CMT following their guidance.


## Optional GPU backends

`cal_ccf` supports CPU (default), Apple Metal and NVIDIA CUDA. CUDA is opt-in at
build time (`DELTAP_ENABLE_CUDA=ON`) and runtime (`AUTOFOCUSING_BACKEND=cuda`).
See [CUDA setup](manual.md#12-nvidia-cuda-gpu-backend) and
[validation results](docs/cuda-performance.md). Metal settings remain supported.

CPU Hessian evaluation and RT rotation are also parallelized with OpenMP.
See [CPU parallelization and validation](docs/cpu-fitting-parallel.md) for
thread settings and numerical checks. The [M4 Max Metal measurements](docs/metal-cpu-parallel-performance-20260914_en.md)
show 1.60× speedup from CPU parallelization and 8.89× cumulative speedup against
the pre-optimization CPU implementation, with the remaining validation scope documented.


## Current behavior and limits

This summary describes source commit `5153085`; it is not a new release designation.

- Spectral matrices use the corrected complex same-station subtraction. Peak search
  includes the boundary and candidate-skipping fixes.
- Newton 4×4 arithmetic and stopping behavior use the restored FP32 path. Objective
  and Hessian accumulation remain double precision, and the parameter-covariance
  inverse uses FP64. The experimental FP64 Newton roundoff fallback is not active.
- Early event selection is optional: `all` retains full-catalog behavior; `selected`
  applies the configured initial-grid max/MAD thresholds before fitting and bootstrap.
- Production bootstrap seeds depend on time, and FFT-plan differences can affect
  numerical results. Exact comparisons in the records use stated deterministic
  controls. Weak candidates can still have unusable covariance estimates; finite
  values alone do not establish a valid uncertainty interpretation.
- CPU, Metal and CUDA validation covers the revisions, inputs and modes identified
  in each record. Historical speedups and agreement checks are not measurements of
  every later revision or GPU configuration.

The calculation records linked below provide the evidence and limitations for
these statements. Runtime options and input/output contracts are maintained in
the manual and script guide.

## Documentation guide

### Usage

- [Manual](manual.md): build, input layout, outputs and CPU/GPU configuration.
- [Workspace and script guide](Scripts/README.md): experiment setup, dates, frequency
  bands, early selection and run records.
- [Global CMT conversion](Scripts/GlobalCMT/README.md): external catalog preparation.

### Current calculation behavior and limitations

| Topic | Validation record |
| --- | --- |
| Complex spectral-matrix correction | [English](docs/spectral-matrix-3c-2005-verification_en.md) / [日本語](docs/spectral-matrix-3c-2005-verification_jp.md) |
| Peak-search corrections | [English](docs/peak-search-validation-20260917_en.md) / [日本語](docs/peak-search-validation-20260917_jp.md) |
| FP64 parameter covariance | [English](docs/parameter-covariance-precision-20260917_en.md) / [日本語](docs/parameter-covariance-precision-20260917_jp.md) |
| Adopted FP32 Newton restoration | [English](docs/newton-fp32-rollback-20260917_en.md) / [日本語](docs/newton-fp32-rollback-20260917_jp.md) |
| Optional early event selection | [English](docs/early-event-selection-20260917_en.md) / [日本語](docs/early-event-selection-20260917_jp.md) |

### Performance and validation history

| Topic | Record |
| --- | --- |
| Horizontal components | [Verification](docs/horizontal-verification.md) |
| Analysis frequency band | [Verification](docs/frequency-band-verification.md) |
| Workspace and dates | [Recorded verification](Scripts/README.md#workspace-verification-recorded-on-2026-09-16) |
| Slowness range and initial seeds | [Current behavior](Scripts/README.md#slowness-range-and-initial-seeds); [earlier inspection](docs/peak-search-validation-20260917_en.md#earlier-inspection-2026-09-16) |
| CPU slant stacking | [Performance](docs/cpu-slant-stack-performance.md) |
| CPU fitting and bootstrap parallelism | [Performance and validation](docs/cpu-fitting-parallel.md) |
| Shared CPU I/O | [English](docs/cpu-io-optimization_en.md) / [日本語](docs/cpu-io-optimization_jp.md) |
| Metal slant stacking | [Performance](docs/metal-slant-stack-performance.md) |
| Metal objective evaluation | [English](docs/metal-power-optimization_en.md) / [日本語](docs/metal-power-optimization_jp.md) |
| CUDA | [Implementation and validation](docs/cuda-performance.md) |
| Metal integration-branch check | [Recorded check](docs/metal-main-qualification-20260914_en.md#earlier-integration-branch-check-2026-09-14) |
| CPU parallelism with Metal | [English](docs/metal-cpu-parallel-performance-20260914_en.md) / [日本語](docs/metal-cpu-parallel-performance-20260914_jp.md) |
| Metal five-day and single-thread checks | [English](docs/metal-main-qualification-20260914_en.md) / [日本語](docs/metal-main-qualification-20260914_jp.md) |
| Three-component CPU/Metal comparison | [English](docs/cpu-gpu-3c-2005-verification_en.md) / [日本語](docs/cpu-gpu-3c-2005-verification_jp.md) |
| Matrix ratios after threshold selection | [English](docs/matrix-ratio-threshold-2005_en.md) / [日本語](docs/matrix-ratio-threshold-2005_jp.md) |
| Matrix-estimate distributions | [English](docs/matrix-estimate-distributions_en.md) / [日本語](docs/matrix-estimate-distributions_jp.md) |
| Historical catalog agreement | [English](docs/historical-catalog-agreement_en.md) / [日本語](docs/historical-catalog-agreement_jp.md) |
| Catalog consistency investigation | [English](docs/catalog-consistency-investigation_en.md) / [日本語](docs/catalog-consistency-investigation_jp.md) |
| FP64 Newton experiment, subsequently superseded | [English](docs/newton-precision-20260917_en.md) / [日本語](docs/newton-precision-20260917_jp.md) |
| Roundoff fallback, subsequently removed | [English](docs/newton-roundoff-20260917_en.md) / [日本語](docs/newton-roundoff-20260917_jp.md) |
| FFT-plan repeatability investigation | [English](docs/pre-bootstrap-repeatability-20260917_en.md) / [日本語](docs/pre-bootstrap-repeatability-20260917_jp.md) |
| Bootstrap/peak-search integration before the FP32 restoration | [English / 日本語](docs/cpu-fitting-parallel.md#peak-search-integration-validation-2026-09-17) |

### Reading the validation records

Records preserve the observations, failures and branch/install status at the time
of each experiment. Use their stated revisions, dates and conditions; a branch
name alone does not identify an exact tested commit. Later adoption notes do not
extend the scope of the original measurements.

- **Public-source checks:** the CMake/CTest suite provides synthetic and script
  checks after building with the required dependencies. GPU executable checks
  additionally require the corresponding hardware. See the manual for commands.
- **External-data checks:** real-waveform comparisons require separately obtained
  waveform and catalog inputs matching the recorded component mode and period.
- **Historical/local reproductions:** commands using old commits require a Git
  clone containing those commits, rather than a source ZIP alone. Saved probe
  sources, deterministic audit setups, catalogs and scripts under
  `docs/benchmarks/` are not part of the current source distribution. A command
  depending on those artifacts cannot be reproduced from public source alone;
  its record identifies the necessary local prerequisites. Historical commands
  may also require the original source revision rather than the latest checkout.

For the contents of generated run directories and what to review before sharing,
see [run records](Scripts/README.md#what-each-run-records).
