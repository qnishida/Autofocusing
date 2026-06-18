# DeltaP

DeltaP is a C++ implementation of the auto-focusing workflow used to estimate
source locations from seismic array data in HDF5 format.

This repository contains the analysis code and small helper scripts. Large
waveform inputs, Global CMT catalog files, converted earthquake catalogs, and
generated outputs are not distributed as part of the source tree.

For detailed usage, input layout, and troubleshooting notes, see `manual.md`.

## Requirements

The build uses CMake and C++11. Dependency discovery uses `find_package` where
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

Set the input paths and run:

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

When `local_config.sh` exists, `run.sh` loads it automatically.

The direct command is:

```bash
./bin/cal_ccf_gcc 2004 param_set_A "$(git describe --tags --always)" \
  /path/to/hdf5 /path/to/moment_loc_76_24
```

The first argument, `YYYY`, is the analysis start year. For example, `2004`
starts the daily scan at `2004-01-01`. The current `cal_ccf_gcc` workflow scans
forward day by day until the built-in end date.

For each day, the program looks for a Hi-net directory under:

```text
HINET_ROOT/YYYY/MMDD/
```

Missing day directories are skipped. If a day directory exists, it is processed
only when the directory contains exactly one regular HDF5 file. Empty day
directories, day directories with multiple files, and missing files are skipped.
Malformed or unreadable HDF5 files may still cause HDF5 errors instead of being
treated as ordinary missing data.

## Output

`run.sh` writes results under:

```text
output/<param-id>/<git-version>/
```

The main output file name follows:

```text
YYYY_<segment length>_<min frequency>-<max frequency>.dat
```

Here `YYYY` is the analysis start year passed to `cal_ccf_gcc`. For the default
`run.sh` settings, `<param-id>` is `param_set_A` and `<git-version>` is resolved
with `git describe --tags --always`.

## Repository Notes

- `Scripts/GlobalCMT/` contains only the conversion helper and instructions
  needed to regenerate the local CMT catalog input.
- `local_config.sh`, `bin/`, `build/`, waveform data, generated catalogs, and
  output products are intentionally not tracked.

## License and Attribution

DeltaP is distributed under the GNU General Public License, version 2 or, at
your option, any later version. See `LICENSE`.

Third-party code, external library licenses, and data provenance notes are
summarized in `NOTICE`. Global CMT raw catalog files and converted DeltaP
catalog files are not distributed under the DeltaP software license. Use
`Scripts/GlobalCMT/` to regenerate the local catalog input from Global CMT data
and cite Global CMT following their guidance.
