# Documentation guide

See the [project overview](../README.md) and
[current behavior and limits](../README.md#current-behavior-and-limits).

## Usage

- [Manual](../manual.md): build, input layout, outputs and CPU/GPU configuration.
- [Workspace and script guide](../Scripts/README.md): experiment setup, dates, frequency
  bands, early selection and run records.
- [Global CMT conversion](../Scripts/GlobalCMT/README.md): external catalog preparation.

## Current calculation behavior and limitations

| Topic | Validation record |
| --- | --- |
| Complex spectral-matrix correction | [English](spectral-matrix-3c-2005-verification_en.md) / [日本語](spectral-matrix-3c-2005-verification_jp.md) |
| Peak-search corrections | [English](peak-search-validation-20260917_en.md) / [日本語](peak-search-validation-20260917_jp.md) |
| FP64 parameter covariance | [English](parameter-covariance-precision-20260917_en.md) / [日本語](parameter-covariance-precision-20260917_jp.md) |
| Adopted FP32 Newton restoration | [English](newton-fp32-rollback-20260917_en.md) / [日本語](newton-fp32-rollback-20260917_jp.md) |
| Optional early event selection | [English](early-event-selection-20260917_en.md) / [日本語](early-event-selection-20260917_jp.md) |

## Performance and validation history

| Topic | Record |
| --- | --- |
| Horizontal components | [Verification](horizontal-verification.md) |
| Analysis frequency band | [Verification](frequency-band-verification.md) |
| Workspace and dates | [Recorded verification](../Scripts/README.md#workspace-verification-recorded-on-2026-09-16) |
| Slowness range and initial seeds | [Current behavior](../Scripts/README.md#slowness-range-and-initial-seeds); [earlier inspection](peak-search-validation-20260917_en.md#earlier-inspection-2026-09-16) |
| CPU slant stacking | [Performance](cpu-slant-stack-performance.md) |
| CPU fitting and bootstrap parallelism | [Performance and validation](cpu-fitting-parallel.md) |
| Shared CPU I/O | [English](cpu-io-optimization_en.md) / [日本語](cpu-io-optimization_jp.md) |
| Metal slant stacking | [Performance](metal-slant-stack-performance.md) |
| Metal objective evaluation | [English](metal-power-optimization_en.md) / [日本語](metal-power-optimization_jp.md) |
| CUDA | [Implementation and validation](cuda-performance.md) |
| hiyodori GPU and OpenMP measurements, 2026-10-22 | [日本語](hiyodori-performance-20261022_jp.md) / [Summary JSON](hiyodori-performance-20261022.json) |
| Metal integration-branch check | [Recorded check](metal-main-qualification-20260914_en.md#earlier-integration-branch-check-2026-09-14) |
| CPU parallelism with Metal | [English](metal-cpu-parallel-performance-20260914_en.md) / [日本語](metal-cpu-parallel-performance-20260914_jp.md) |
| Metal five-day and single-thread checks | [English](metal-main-qualification-20260914_en.md) / [日本語](metal-main-qualification-20260914_jp.md) |
| Three-component CPU/Metal comparison | [English](cpu-gpu-3c-2005-verification_en.md) / [日本語](cpu-gpu-3c-2005-verification_jp.md) |
| Matrix ratios after threshold selection | [English](matrix-ratio-threshold-2005_en.md) / [日本語](matrix-ratio-threshold-2005_jp.md) |
| Matrix-estimate distributions | [English](matrix-estimate-distributions_en.md) / [日本語](matrix-estimate-distributions_jp.md) |
| Historical catalog agreement | [English](historical-catalog-agreement_en.md) / [日本語](historical-catalog-agreement_jp.md) |
| Catalog consistency investigation | [English](catalog-consistency-investigation_en.md) / [日本語](catalog-consistency-investigation_jp.md) |
| FP64 Newton experiment, subsequently superseded | [English](newton-precision-20260917_en.md) / [日本語](newton-precision-20260917_jp.md) |
| Roundoff fallback, subsequently removed | [English](newton-roundoff-20260917_en.md) / [日本語](newton-roundoff-20260917_jp.md) |
| FFT-plan repeatability investigation | [English](pre-bootstrap-repeatability-20260917_en.md) / [日本語](pre-bootstrap-repeatability-20260917_jp.md) |
| Bootstrap/peak-search integration before the FP32 restoration | [English / 日本語](cpu-fitting-parallel.md#peak-search-integration-validation-2026-09-17) |

## Reading the validation records

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
see [run records](../Scripts/README.md#what-each-run-records).
