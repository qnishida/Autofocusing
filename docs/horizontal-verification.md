# Horizontal-only implementation and verification

See the [manual](../manual.md) for current CLI and date handling. The original
horizontal implementation was compared against source revision `149995f`.

## Input contract

Horizontal inputs must already be velocity waveforms in sensor coordinates.
The program does not perform a tilt-to-velocity conversion. Apply the HDF5
`cmpaz` rotation once, after the FFT. LE/LN is preferred to E/N, and a complete
orthogonal pair is required. Dataset extents, rate, sensitivity, velocity units
and component azimuths are checked before accepting a station.

The velocity conversion remains `integer * sensitivity * 1e-9` for `nm/s` input.
Malformed or missing components cannot produce an accepted station. Three-component
unit behavior retains compatibility with legacy inputs lacking a unit attribute.
Waveforms and the earthquake catalog are external inputs.

## Implementation decisions

The direct CLI accepts a component mode, with `3c` as its legacy default.
The workspace launcher defaults to horizontal mode unless configured otherwise.
Horizontal spectra supply amplitude and stability checks; the latter uses mean
E/N band power. U candidate search and U-derived horizontal seeds are skipped.
R/T event numbers remain 0/1. Output keeps 38 columns, with U-related columns
36–38 explicitly `nan`. Missing U buffers and weights are initialized.

The `prm_init` fix initializes fit-change diagnostics from the incoming parameters.
It corrects columns 9–10 without changing the optimization path. Matrix arithmetic
has subsequent corrections documented in the
[spectral-matrix validation](spectral-matrix-3c-2005-verification_en.md).

## Build and automated checks

Use compatible compiler and library builds as described in the manual. The
`horizontal_components` CTest covers LE/LN reading, E/N fallback, scaling/filtering,
independent azimuth rotation, invalid metadata, preserved U in 3c mode, horizontal
QC, geographic R/T rotation and the U-missing matrix. Launcher tests cover path
joining, external result placement and direct-CLI compatibility.

```bash
ctest --test-dir build -R 'horizontal_components|parent_launcher' --output-on-failure
```

The historical synthetic 3c loader comparison against `149995f` checked 1599
spectral/band-power values with `rtol=1e-10`, `atol=1e-18`. FFTW plan-dependent
rounding is allowed; this is not a full event-catalog regression.

## Validation summary and limits

Historical real-data smoke checks verified loading, accepted windows and the
38-column output contract. They demonstrated execution and finite diagnostics,
not independent confirmation of physical source locations or full-archive accuracy.

Gap metadata does not identify every missing sample. QC thresholds, earthquake
exclusion and array geometry need dataset-specific review. Use explicit dates,
a new run identifier and versioned inputs; consult the manual for current date
bounds and provenance rather than relying on historical deployment settings.
Detailed per-day smoke results and deployment records are maintained separately.
