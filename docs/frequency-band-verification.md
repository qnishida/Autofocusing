# Configurable frequency band — verification, 2026-09-16

`cal_ccf` accepts `AUTOFOCUSING_FREQ_MIN` / `AUTOFOCUSING_FREQ_MAX` in Hz,
with defaults 0.1 / 0.25. Experiment configuration is recorded by the runner,
which queries `--frequency-info` before allocating output to reject invalid
settings and unsupported older binaries. The returned FFT bins and effective
band are stored in `manifest.json`; startup emits the same JSON in `run.log`.

Both endpoints retain legacy floor-to-bin conversion at 1/1024 Hz. For the
primary-microseism band 0.05–0.1 Hz, bins 51–102 are included, corresponding to
0.0498046875–0.099609375 Hz. Default bins remain 102–256. Values must be finite,
positive, ordered, below 1 Hz Nyquist, and resolve above DC. A single-bin band
is allowed when distinct requested endpoints round to the same bin.

The analysis buffers, fitting and peak-frequency diagnostic use the selected
band. Spectrum storage retains at least 266 bins for fixed QC and grows for a
higher selected band. QC's upper limit is independently capped at its previous
value, so extra analysis bins do not change QC power. The peak-frequency
calculation now excludes stored bins above the analysis band, including on
otherwise default runs; its zero-power result stays inside the selected band.

## Verification

Built `cal_ccf_clang` and the affected C++ tests in `build-clang/`, retaining
Metal support and the documented Homebrew C/C++ / Apple Objective-C++ compiler
configuration. Ran:

```bash
ctest --test-dir build-clang -R 'frequency_|cpu_parallel_equivalence|io_equivalence|horizontal_components|slant_stack_equivalence|workspace_scripts|parent_launcher|analysis_date_range' --output-on-failure
```

All 9 CTest entries passed:

- Frequency processing: default/primary FFT bins, fixed QC power for synthetic
  multi-frequency waveforms at low and extended analysis bands, sufficient
  FFT storage, and exclusion of strong out-of-band power from fitting and peak
  frequency (including endpoint and zero-power cases).
- Frequency configuration (3 Python cases): metadata/log/filename agreement,
  invalid/nonfinite/reversed/DC/Nyquist settings rejected before run allocation,
  near-Nyquist and single-bin conversion.
- Slant-stack reference equivalence: both default and primary bands, horizontal
  and 3c modes, offsets, multiple thread counts and repeated accumulation.
- Existing CPU fitting, horizontal component and I/O equivalence tests.
- Workspace integration (10 cases), parent launcher and date-range regressions.
  The workspace tests also cover config precedence, frequency metadata and an
  unsupported binary being rejected without starting an analysis.

Installed the tested executable with `cmake --install build-clang`. The installed
`bin/cal_ccf_clang --frequency-info`, with 0.05/0.1 supplied in the environment,
returned the expected bins and effective endpoints. `git diff --check` passed.

## Limits

Tests used CPU calculations, synthetic waveforms and temporary empty archives;
no real-data primary-microseism analysis or GPU numerical qualification was run.
QC/stability bands and thresholds, 0.03 Hz high-pass preprocessing and earthquake
exclusion remain unchanged; changing frequency does not retune them. Existing
experiment configs were preserved and need explicit frequency settings.
`cal_ccf_eq` is unchanged. The separate slowness peak-selection issues recorded
in [slowness-seed-selection.md](slowness-seed-selection.md) were not modified.
