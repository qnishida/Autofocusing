# DeltaP Manual

This manual describes the practical steps needed to build and run DeltaP with
external Hi-net waveform data and a Global CMT-derived earthquake catalog.

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

DeltaP requires two external inputs:

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

## 5. Running DeltaP

The standard wrapper is `run.sh`:

```bash
HINET_ROOT=/path/to/hdf5 \
CMT_CATALOG=/path/to/moment_loc_76_24 \
./run.sh
```

For repeated local runs, create a git-ignored `local_config.sh` in the
repository root:

```bash
export HINET_ROOT="/path/to/hdf5"
export CMT_CATALOG="/path/to/moment_loc_76_24"
```

Then run:

```bash
./run.sh
```

The direct command is:

```bash
./bin/cal_ccf_gcc 2004 param_set_A "$(git describe --tags --always)" \
  /path/to/hdf5 /path/to/moment_loc_76_24
```

The first argument, `YYYY`, is the analysis start year. For example, `2004`
starts the daily scan at `2004-01-01`. The current `cal_ccf_gcc` workflow scans
forward day by day until the built-in end date.

## 6. Output

`run.sh` writes results under:

```text
output/<param-id>/<git-version>/
```

The main output file name follows:

```text
YYYY_<segment length>_<min frequency>-<max frequency>.dat
```

Here `YYYY` is the analysis start year passed to `cal_ccf_gcc`. With the
default `run.sh`, `<param-id>` is `param_set_A`, and `<git-version>` is resolved
with `git describe --tags --always`.

Generated outputs are not distributed as part of the source repository. If
final data products are published separately, include their provenance,
generation command, column definitions, units, and relation to the manuscript.

## 7. Troubleshooting

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

If `./run.sh` cannot find `./bin/cal_ccf_gcc`, run the install step:

```bash
cmake --install build
```

If no output is produced for expected days, check that each day directory exists
and contains exactly one regular HDF5 file under `HINET_ROOT/YYYY/MMDD/`.
