# Global CMT catalog conversion

This directory contains the small conversion utility used to prepare the
earthquake catalog input for Autofocusing from the Global Centroid Moment Tensor
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

## Convert NDK to Autofocusing catalog format

The C++ code expects a whitespace-separated catalog with this format:

```text
yy jday hour minute second longitude latitude depth moment
```

From the repository root, generate it from one or more concatenated NDK files:

```bash
cat /path/to/*.ndk | perl Scripts/GlobalCMT/ndk_to_moment_loc.pl > /path/to/moment_loc_76_24
```

Then prepare an experiment using the [workspace setup guide](../README.md#setup).
Set `HINET_ROOT` and `CMT_CATALOG` in the workspace's `local_config.sh`, and set
both `START_DATE` and `END_DATE` in the experiment's `config.sh`. Configuration
files override inherited environment values, so configure the files actually
used by the runner. From the repository root:

```bash
bash Scripts/setup_workspace.sh --experiment primary-microseisms
# Edit ../local_config.sh and ../analysis/primary-microseisms/config.sh.
bash Scripts/run.sh --experiment primary-microseisms
```

The experiment name does not select a frequency band; choose the analysis
settings explicitly as described in the workspace guide.

Alternatively, from the repository root, call the executable directly with an
explicit date interval and component mode:

```bash
./bin/cal_ccf_gcc 2004 param_set_A example-run \
  /path/to/hdf5 /path/to/moment_loc_76_24 ../results 3c 2004-01-01 2004-01-07
```

Use `cal_ccf_clang` for a Clang build. `example-run` is an output-directory label,
so the command also works with a source ZIP without Git history. Use a new label
for each direct run and record the source release or commit separately. The
waveform and catalog inputs must cover the selected interval.
