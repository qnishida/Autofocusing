# CPU / Metal agreement for real three-component data

Detailed run logs, output catalogs and per-event comparison records are retained
locally and are not distributed in this repository. References to these records
below describe local artifacts.

Verified 2026-09-17 on Apple M4 Max. For 2005-03-15 through 2005-03-17,
CPU and Metal produced the same 60 events, fitted parameters, bootstrap fields,
and spectral-matrix fields at catalog output precision. Both the original and
corrected matrix formulas passed the comparison independently. No difference
in the energy ratios computed from these printed matrix values was observed.

## Controlled comparison

- Requested frequency 0.1–0.25 Hz; actual bins 102–256 at 1/1024 Hz spacing
  (0.099609375–0.25 Hz, 155 bins).
- Slowness coordinates ±0.165 s/km, step 0.005 s/km, 67 × 67 grid;
  existing radial mask retained; approximately 684–685 stations, three components.
- Identical matrix-audit executable, waveform and source-catalog hashes, date
  boundaries, eight CPU threads, fixed seeds and FFTW ESTIMATE | UNALIGNED plans.
  These deterministic controls are test-only changes; production binaries were
  not replaced. Identical start dates avoid the restart-state confound described
  in the [catalog investigation](catalog-consistency-investigation_en.md).
- Explicit Metal backend, GPU power off, per-window CPU/GPU verification enabled.
  The log confirms GPU execution; all 520 accepted windows were checked.
- Predeclared catalog tolerances: max and MAD relative tolerance 1e-4; convergence
  residual absolute tolerance 1e-12; all other fields exact printed strings.

## Results

| Quantity | CPU versus Metal |
|---|---|
| Events | 60 each, same component/order/time |
| Per-window slant-stack peak mismatches | 0 / 520 |
| Fitted location and slowness | Exact at output precision |
| Bootstrap fields and scalar MS_max | Exact at output precision |
| SRR, STT, SUU, Re(SRU), Im(SRU) | Exact at output precision |
| Maximum relative difference in max | 0.0006852% |
| Maximum relative difference in MAD | 0.0017671% |
| Maximum relative difference in max/MAD | 0.0017671% |
| Selected R / T / U events | 28 / 10 / 15 in both backends |

Selection uses strict max/MAD > 7 for R/T and > 35 for U (vertical/V).
Only max (2 rows) and MAD (50 rows) differ among the 38 printed columns.
The table applies to both matrix formulas. Event relative differences use
100 × abs(GPU − CPU) / abs(CPU), computed from printed values.

Across the 520 windows, the largest slant-stack max-scaled error was
1.96957e-6 (0.0001970%); the largest relative L2 error was 4.44972e-6.
Max-scaled error divides the largest absolute GPU−CPU grid difference by the
largest absolute CPU grid value across components. It is not a pointwise
relative error, an SUU error bar, or event-to-event variability.

A separate 700-station, 155-bin synthetic kernel check also passed: the 3c
max-scaled error was 1.14745232e-5 (0.0011475%), relative L2 9.64068053e-6,
with no peak mismatches. The horizontal case passed as well; see the raw record.

## Scope and interpretation

Metal evaluates the initial slant stack in float, versus the CPU double kernel.
In three-component mode, subsequent fitting, bootstrap and matrix calculations
use the shared CPU path. That shared path itself includes float matrices, so
this is not a comparison of an entirely single-precision program with an
entirely double-precision program. Equal output values do not prove bitwise
identity of unprinted intermediate quantities.

The matrix subtraction bug affects both backends; CPU/GPU agreement does not
validate the original formula. The observed energy shifts after its correction
are separate from GPU rounding. This test found no CPU/GPU energy-ratio shift
in these 60 events, but does not cover CUDA, all dates, or all operating modes.
Verification computes both kernels and is not a performance benchmark.

## Records and reproduction

- GPU comparison report,
  log,
  original-formula catalog,
  corrected catalog.
- CPU baseline.
- 700-station kernel check.
- Artifact hashes and event-difference statistics.

With a local passing CPU audit directory, its matching binary and input files,
use a fresh output path. Replace CPU_RUN below with that directory:

```bash
python3 -B tests/run_cpu_gpu_events.py \
  CPU_RUN \
  build-clang/cpu-gpu-3c-repeat
./build-clang/test_metal_slant_stack 700 33 4 1
```

The script verifies executable and input hashes against the CPU baseline and
requires access to the Metal device. Logs retain original runtime paths;
detailed artifacts remain in the local run directories.
