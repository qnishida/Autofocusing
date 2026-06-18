# Global CMT catalog conversion

This directory contains the small conversion utility used to prepare the
earthquake catalog input for DeltaP from the Global Centroid Moment Tensor
(Global CMT) catalog.

The Global CMT catalog itself, downloaded NDK files, and converted catalog
products are not distributed in this repository. The repository software
license applies to the scripts and source code here, not to Global CMT data or
derived data products.

Converted catalog files are intentionally not distributed. Regenerate the
catalog from Global CMT NDK files and pass the generated file via
`CMT_CATALOG`.

## Input data

Download the required NDK catalog files from Global CMT:

- https://www.globalcmt.org/CMTfiles.html

For the published analysis, record the exact source files, date downloaded, and
time span used in your run notes. Cite Global CMT following their citation
guidance:

- https://www.globalcmt.org/CMTcite.html

## Convert NDK to DeltaP catalog format

The C++ code expects a whitespace-separated catalog with this format:

```text
yy jday hour minute second longitude latitude depth moment
```

Generate it from one or more concatenated NDK files:

```bash
cd Scripts/GlobalCMT
cat /path/to/*.ndk | perl ndk_to_moment_loc.pl > /path/to/moment_loc_76_24
```

Then pass the generated file to DeltaP:

```bash
HINET_ROOT=/path/to/hdf5 \
CMT_CATALOG=/path/to/moment_loc_76_24 \
./run.sh
```

or call the executable directly:

```bash
./bin/cal_ccf_gcc 2004 param_set_A "$(git describe --tags --always)" \
  /path/to/hdf5 /path/to/moment_loc_76_24
```
