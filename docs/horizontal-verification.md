# Horizontal-only implementation and verification

Verified 2026-09-08 on branch `feature/horizontal-two-component`, based on
`149995f`. Verification was performed before committing the changes.
No full-archive analysis was started.

## Input contract and local workspace

- Waveforms: `/Volumes/Seismic_Data/hdf5/Hi-net_tilt/YYYY/MMDD/*.h5`.
- CMT: `../moment_loc_76_24`, supplied by the user, outside the source repository.
- The user confirmed that the data are velocity converted from tilt, still in
  sensor coordinates. Apply the HDF5 `cmpaz` rotation once, after the FFT.
- On January 1 of 2004, 2014 and 2024, both N.AGMH and N.AGWH have LE/LN,
  2 Hz, 172800 integer samples, sensitivity 0.1, unit `nm/s`, and gap=0.
  LE/LN azimuths are 97/7 degrees and 93/3 degrees respectively.
- Conversion to m/s remains `integer * sensitivity * 1e-9`; no tilt-to-velocity
  operation is performed by Autofocusing.
- Parent launcher: `../run.sh`; configuration: `../local_config.sh`; production
  destination: `../results/tilt_horizontal/<git-version>/`. Tracked deployment
  originals are in `Scripts/`.

## Implementation decisions

Direct CLI gains optional `[output-root] [3c|horizontal]`, defaulting to
`output/` and `3c`. The parent launcher selects horizontal mode. LE/LN is
preferred to E/N; complete component pairs are required. Metadata, dataset
extent, rate, velocity unit, and orthogonal azimuths are checked. Invalid
components cannot produce an accepted station. Three-component unit behavior
is retained for compatibility with legacy inputs without a unit attribute.

Horizontal spectra drive amplitude and stability checks; the latter uses
mean E/N band power. Thresholds, array selection, frequency band, window
length, CMT exclusion, and core search/optimization are retained. U candidate
search and U-derived horizontal seeds are skipped. R/T event numbers remain
0/1. The output retains 38 columns, with U-related columns 36–38 explicitly
`nan`. Component buffers are initialized so missing U cannot introduce
uninitialized values or nonzero weights.

The existing `prm_init` in gradient fitting was uninitialized. It is now
initialized from the incoming parameters: before this fix, real-data output
had unintended `nan` in the fit-change diagnostics (columns 9–10). This fixes
these diagnostics in both modes and does not alter the optimization path.

## Build and automated checks

Use the Homebrew LLVM/Clang procedure in `manual.md`. C++14 is needed by the
installed dependencies. C is enabled for HDF5 detection, AppleClang is named
correctly, and the obsolete separate Boost.System dependency is removed.
Explicit filesystem and standard-library includes fix modern Boost/compiler
compatibility; macOS skips the unavailable POSIX prefetch hint.

A GCC build linked against this Mac's Homebrew Boost compiled but failed
filesystem behavior at runtime (no days found, output in the wrong location).
The Clang build worked. Tests check path joining and actual output placement
to catch this toolchain mismatch. The installed executable is
`bin/cal_ccf_clang`; the launcher prefers it when present.

- `ctest --test-dir build-clang --output-on-failure`: **2/2 passed**.
- `horizontal_components`: LE/LN reading, E/N fallback, velocity scaling/filter,
  independent azimuth rotation reference, missing/malformed components,
  unsupported units/rate/nonorthogonal axes, preserved U in 3c, horizontal
  stability without U, geographic R/T rotation and U-missing beam matrix.
- `parent_launcher`: deployed layout invoked from another cwd, relative
  configuration paths, actual executable output under results, missing CMT,
  invalid mode, and legacy direct CLI destination.
- `python3 tests/compare_legacy.py build-clang 149995f`: **passed**, 1599
  spectral and band-power values for a valid 3c fixture match the pre-change
  loader with rtol=1e-10, atol=1e-18. FFTW_MEASURE may choose different plans;
  comparison therefore allows floating-point rounding, rather than requiring
  byte identity. This is not a full 3c event-catalog regression.
- Shell syntax and `git diff --check`: passed.

## Real-data smoke tests

All reads were read-only. The 2014 and 2024 tests used
`build-clang/test_horizontal --real <file>`.

| Day | Loaded stations | First-window spectra | First-window amplitude QC |
| --- | ---: | ---: | ---: |
| 2014-01-01 | 725 | 725 | 683 |
| 2024-01-01 | 725 | 725 | 647 |

A one-day input tree under `build-gcc/one-day/2004/0101/` contains only a
symlink to the original 2004-01-01 HDF5 file. The `build-gcc` path is just the
scratch input location; final verification used the Clang executable:

```bash
OMP_NUM_THREADS=4 bin/cal_ccf_clang 2004 tilt_horizontal final-verification \
  build-gcc/one-day ../moment_loc_76_24 build-clang/one-day-results horizontal
```

The process exited 0, loaded 653 stations, and accepted 47, 19, 26 and 42
windows across its four segments (134 total, following the existing segment
iteration). It emitted two R events and no T events on this day. Each output
row has 38 columns; all numeric values through column 35 are finite, and
columns 36–38 are `nan`. Fit-change diagnostics are finite after the
initialization fix. This demonstrates execution and output integrity, not
independent confirmation of the physical source locations.

The verified launcher and configuration were deployed to the parent directory
without overwriting existing files. The one-day `.dat` result and `run.log`
are preserved in `../results/tilt_horizontal/verification-20260908/` (relative
to the repository root). The launcher and configuration match their tracked
originals. Production analysis has not been started.

## Limits and follow-up context

- Start year remains 2004 by default; the existing 2024-12-31 cutoff and
  366 × 20.75-day iteration bound are unchanged. No 2025 support was added.
- Gap metadata does not specify sample-level missing intervals. Existing gap
  and threshold behavior is retained and warrants dataset-specific review
  before production interpretation.
- Existing search behavior was not redesigned. In particular, `ssRTU` is
  cleared once per day and reused across segments, and the legacy beam-matrix
  bias loop subtracts cross-component terms from diagonal entries. These
  deserve a separate scientific review; this change does not claim to validate
  their statistical interpretation or the complete existing algorithm.
- The same parameter/version/start-year combination overwrites its result
  file. `-dirty` records tracked uncommitted changes but is not a unique
  identifier for successive uncommitted versions.
