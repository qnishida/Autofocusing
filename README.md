# Autofocusing

Autofocusing is a C++ implementation of the auto-focusing workflow used to
estimate source locations from seismic array data in HDF5 format.

This repository contains the analysis code and small helper scripts. Large
waveform inputs, Global CMT raw catalog files and CMT catalogs converted 
for use with Autofocusing are not distributed under the Autofocusing 
software license.

For detailed usage, input layout, and troubleshooting notes, see `manual.md`.

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

The local analysis workspace uses this layout:

```text
Autofocusing/
├── run.sh
├── local_config.sh
├── moment_loc_76_24
├── repo/
└── results/
```

Only `repo/` is a Git repository. Its root on GitHub contains `src/`,
`Scripts/`, and the other source files directly; `repo/` is a local checkout
directory name, not another directory to commit. Do not initialize Git in
the parent `Autofocusing/` workspace.

From `repo/`, copy the templates once, preserving existing analysis files:

```bash
cp -n Scripts/run.sh ../run.sh
cp -n Scripts/local_config.example.sh ../local_config.sh
../run.sh
```

The parent `run.sh` is an editable regular file, not a symlink. Keep machine
paths and backend/thread settings in `local_config.sh`, and customize the
parent launcher when an analysis needs different execution steps. Git pulls
update the templates inside `repo/` only. Compare template changes with
`diff -u Scripts/run.sh ../run.sh` and incorporate the relevant changes manually;
do not overwrite a customized parent launcher or configuration.

Run Git commands from `repo/`. Moving a checkout on disk does not publish any
changes to GitHub: intentional source changes need a commit and a separate
push. Before updating a checkout, inspect `git status` and preserve unfinished
work; then use `git fetch origin` to inspect upstream changes and
`git pull --ff-only` when the working tree is ready. Keep machine-specific
settings, catalogs and results in the parent workspace. See
[workspace management](manual.md#13-workspace-and-git-management).

The launcher defaults to horizontal-only analysis of `Hi-net_tilt`, starting
in 2004. Set `HINET_ROOT`, `CMT_CATALOG`, `START_YEAR`, `COMPONENT_MODE`, and
`PARAM_ID` in the parent `local_config.sh`. Relative paths are resolved against
that parent. For three-component data, set `COMPONENT_MODE=3c` and use a distinct
parameter ID. See [manual.md](manual.md) for mode details and macOS build notes.

Direct execution remains available from the repository root:

```bash
./bin/cal_ccf_clang 2004 tilt_horizontal "$(git describe --tags --always --dirty)" \
  /path/to/hdf5 /path/to/moment_loc_76_24 ../results horizontal
```

The compiler determines the executable suffix (`gcc` or `clang`). The two
optional trailing arguments are `[output-root] [3c|horizontal]`; omitting both
preserves three-component mode and the `output/` destination.

The existing daily scan ends at its built-in limit (2024-12-31, also bounded
by 366 × 20.75 days from the start year). It has not been extended to 2025.
Daily input directories use `HINET_ROOT/YYYY/MMDD/` and must contain exactly
one regular HDF5 file. Missing/empty/multiple-file days are skipped.

## Output

The parent launcher writes to:

```text
results/<param-id>/<git-version>/YYYY_2048_0.099609-0.250000.dat
```

Here `YYYY` is the start year, not the year of each individual event.
The default parameter ID is `tilt_horizontal`. Horizontal mode preserves the
38-column format, emits R=0 and T=1 events, and writes `nan` for U-related
matrix values (columns 36–38). Use a separate parameter ID for each experiment;
reusing the same destination overwrites the existing `.dat` file.

Generated build directories, binaries, and outputs are not tracked.
[Verification notes](docs/horizontal-verification.md) describe the input
checks, tests, and remaining limitations.

## Citation

If you use this software, please cite the Zenodo release:

> Nishida, K. (2026). Autofocusing (Version 1.0.0) [Software].
> https://doi.org/10.5281/zenodo.20742021

Citation metadata are also provided in `CITATION.cff`.

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
