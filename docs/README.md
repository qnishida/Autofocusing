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
| Complex spectral-matrix correction | [English](spectral-matrix-3c-2005-verification_en.md) |
| Peak-search corrections | [English](peak-search-validation-20260917_en.md) |
| FP64 parameter covariance | [English](parameter-covariance-precision-20260917_en.md) |
| Adopted FP32 Newton restoration | [English](newton-fp32-rollback-20260917_en.md) |
| Optional early event selection | [English](early-event-selection-20260917_en.md) |

## Performance and validation history

Public documentation is maintained in English.
Machine-specific measurement reports, benchmark plans, and real-data diagnostic
scripts belong in the separate analysis workspace. Detailed logs and generated
outputs belong in the workspace's results directory. Keep public records free
of internal paths and machine-local account information.

| Topic | Record |
| --- | --- |
| Horizontal components | [Verification](horizontal-verification.md) |
| Analysis frequency band | [Verification](frequency-band-verification.md) |
| Workspace and dates | [Regression checks](../Scripts/README.md#workspace-and-date-regression-checks) |
| Slowness range and initial seeds | [Current behavior](../Scripts/README.md#slowness-range-and-initial-seeds); [earlier inspection](peak-search-validation-20260917_en.md#earlier-inspection-2026-09-16) |
| CPU slant stacking | [Performance](cpu-slant-stack-performance.md) |
| CPU fitting and bootstrap parallelism | [Performance and validation](cpu-fitting-parallel.md) |
| Shared CPU I/O | [English](cpu-io-optimization_en.md) |
| Metal slant stacking | [Performance](metal-slant-stack-performance.md) |
| Metal objective evaluation | [English](metal-power-optimization_en.md) |
| CUDA | [Implementation and validation](cuda-performance.md) |
| Metal integration-branch check | [Recorded check](metal-main-qualification-20260914_en.md#earlier-integration-branch-check-2026-09-14) |
| CPU parallelism with Metal | [English](metal-cpu-parallel-performance-20260914_en.md) |
| Metal five-day and single-thread checks | [English](metal-main-qualification-20260914_en.md) |
| Three-component CPU/Metal comparison | [English](cpu-gpu-3c-2005-verification_en.md) |
| Matrix ratios after threshold selection | [English](matrix-ratio-threshold-2005_en.md) |
| Matrix-estimate distributions | [English](matrix-estimate-distributions_en.md) |
| Historical catalog agreement | [English](historical-catalog-agreement_en.md) |
| Catalog consistency investigation | [English](catalog-consistency-investigation_en.md) |
| FP64 Newton experiment, subsequently superseded | [English](newton-precision-20260917_en.md) |
| Roundoff fallback, subsequently removed | [English](newton-roundoff-20260917_en.md) |
| FFT-plan repeatability investigation | [English](pre-bootstrap-repeatability-20260917_en.md) |
| Bootstrap/peak-search integration before the FP32 restoration | [Validation](cpu-fitting-parallel.md#peak-search-integration-validation-2026-09-17) |

## Reading the validation records

Public summaries describe implementation behavior, verification methods, selected
aggregate results and their limits. Use the stated source revisions and conditions;
a branch name alone does not identify an exact tested commit. Detailed development
diaries, machine-specific commands and translations are kept outside this repository.

- **Public-source checks:** the CMake/CTest suite provides synthetic and script
  checks after building with the required dependencies. GPU executable checks
  additionally require the corresponding hardware. See the manual for commands.
- **External-data checks:** real-waveform comparisons require separately obtained
  waveform and catalog inputs matching the recorded component mode and period.
- **Historical/local reproductions:** commands using old commits require a Git
  clone containing those commits, rather than a source ZIP alone. Saved
  deterministic audit inputs and experimental drivers
  are not part of the source distribution. Historical replay may require those
  external materials and the exact original revision. Do not interpret an old
  result as a test rerun of the latest checkout.

For the contents of generated run directories and what to review before sharing,
see [run records](../Scripts/README.md#what-each-run-records).
