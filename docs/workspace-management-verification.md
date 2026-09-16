# Workspace management verification — 2026-09-16

The script changes separate experiment definitions from generated runs, retain
source/analysis provenance, and add explicit inclusive date ranges while keeping
year-only direct executable calls compatible. The
workflow and reproduction commands are in [Scripts/README.md](../Scripts/README.md).

## Checks performed

- `bash -n Scripts/run.sh Scripts/setup_workspace.sh`: passed.
- The `tests/workspace.py` suite (also run through CTest): all 9 integration tests passed. Coverage includes
  non-overwriting setup, independent Git initialization, refusal to initialize
  inside another checkout, settings precedence, paths with spaces, source and
  analysis revisions/diffs, uninitialized and unborn repositories, timestamp
  collisions, concurrent runs, missing inputs, child
  failure, Metal power alias precedence, and SIGINT/SIGTERM handling.
  Omitting `--experiment` now returns status 2 with usage and an example before
  sourcing configuration or creating results, even when `PARAM_ID` is set.
  Both generated and direct launchers are covered, along with `--workspace`
  alone, missing/empty experiment values, `--help` (status 0), and the Python
  helper's required argument check.
  Date settings are required, validated before run allocation, passed to the
  executable and saved in metadata; the start year is derived from START_DATE.
- Configured a fresh CPU build in `build-workspace-validation/` using the manual's
  Homebrew LLVM compiler discovery, with Metal/CUDA disabled. Built the
  `cal_ccf_clang` target successfully.
- `ctest --test-dir build-workspace-validation -R 'workspace_scripts|parent_launcher|analysis_date_range' --output-on-failure`:
  all three tests passed. The parent-launcher test used the freshly built executable
  with an empty waveform archive, checked run metadata and relative paths, and
  retained checks of direct CLI compatibility and invalid backend/grid settings.
- The `tests/date_range.py` suite passed all 6 checks against the real executable:
  leap days, year boundaries, a single-day interval, removal of historical caps
  for explicit dates, invalid/reversed dates, skipping missing/empty/multiple-file
  days, inclusive endpoints and exclusion of files outside the selected interval.
  Malformed singleton files intentionally verify that later existing inputs are
  reached and that corrupt data still fails instead of being silently skipped.
- Called the actual parent launcher without arguments from
  `analysis/primary-microseisms/`: usage and an example were displayed and the
  process exited with status 2; no analysis was started.
- `git diff --check`: passed.
- Rebuilt the Metal-enabled Release target in `build-clang/` and installed
  `bin/cal_ccf_clang`. The installed binary passed both
  `python3 -B tests/launcher.py . bin/cal_ccf_clang build-clang` and
  `python3 -B tests/date_range.py bin/cal_ccf_clang build-clang` (6 date tests).
  Tests explicitly selected CPU; Metal was compiled but no GPU analysis was run.

An initial attempt to reuse `build-cpu-only/` failed because its CMake cache
referenced a removed LLVM 22.1.5 OpenMP library. Fresh configuration discovered
LLVM 23.1.1 and resolved the issue without source changes.

For the Metal-enabled build, Homebrew LLVM 23.1.1 Objective-C++ emitted class
selector dispatch symbols that the installed linker could not resolve. Fresh
configuration with Homebrew LLVM for C/C++ and Apple's compiler for Objective-C++
resolved the link failure. The macOS build recipe in `manual.md` now records that
compiler selection. The existing Metal-enabled/CPU-default behavior is retained.

## Limits

No real waveform analysis, GPU calculation, frequency-band change, existing
workspace migration, commit, or publication was performed. Existing parent
launchers/configs/results were preserved. Test workspaces were temporary and
removed automatically. Formal analyses still need an explicitly rebuilt binary
and versioned input data; the recorded checkout SHA does not prove binary build
origin, and waveform contents are not hashed by the launcher.
