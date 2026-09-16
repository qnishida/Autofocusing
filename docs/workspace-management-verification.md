# Workspace management verification — 2026-09-16

The script changes separate experiment definitions from generated runs, retain
source/analysis provenance, and preserve the existing executable CLI. The
workflow and reproduction commands are in [Scripts/README.md](../Scripts/README.md).

## Checks performed

- `bash -n Scripts/run.sh Scripts/setup_workspace.sh`: passed.
- `python3 -B tests/workspace.py`: all 7 integration tests passed. Coverage includes
  non-overwriting setup, independent Git initialization, refusal to initialize
  inside another checkout, settings precedence, paths with spaces, source and
  analysis revisions/diffs, uninitialized and unborn repositories, timestamp
  collisions, concurrent runs, legacy no-argument use, missing inputs, child
  failure, Metal power alias precedence, and SIGINT/SIGTERM handling.
- Configured a fresh CPU build in `build-workspace-validation/` using the manual's
  Homebrew LLVM compiler discovery, with Metal/CUDA disabled. Built the
  `cal_ccf_clang` target successfully.
- `ctest --test-dir build-workspace-validation -R 'workspace_scripts|parent_launcher' --output-on-failure`:
  both tests passed. The parent-launcher test used the freshly built executable
  with an empty waveform archive, checked run metadata and relative paths, and
  retained checks of direct CLI compatibility and invalid backend/grid settings.
- `git diff --check`: passed.

An initial attempt to reuse `build-cpu-only/` failed because its CMake cache
referenced a removed LLVM 22.1.5 OpenMP library. Fresh configuration discovered
LLVM 23.1.1 and resolved the issue without source changes.

## Limits

No real waveform analysis, GPU calculation, frequency-band change, existing
workspace migration, commit, or publication was performed. Existing parent
launchers/configs/results were preserved. Test workspaces were temporary and
removed automatically. Formal analyses still need an explicitly rebuilt binary
and versioned input data; the recorded checkout SHA does not prove binary build
origin, and waveform contents are not hashed by the launcher.
