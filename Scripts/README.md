# Analysis workspaces and recorded runs

The public source checkout, analysis definitions and generated results have
separate lifecycles. Keep the source in `repo/`, track small analysis files in
an independent `analysis/` Git repository, and back up `results/` separately.
Do not initialize Git in the parent workspace.

```text
Autofocusing/
├── repo/                         # Public source repository
├── run.sh                        # Thin launcher created by setup
├── local_config.sh               # Machine settings, outside Git
├── analysis/                     # Optional independent Git repository
│   ├── README.md
│   ├── .gitignore
│   └── primary-microseisms/
│       ├── config.sh
│       └── README.md             # Objective, procedure and interpretation
└── results/
    └── primary-microseisms/
        └── 20260916T053025Z/      # UTC start timestamp; collisions get -01, -02, ...
            ├── manifest.json
            ├── effective_config.json
            ├── local_config.sh
            ├── experiment_config.sh
            ├── run.log
            ├── source.patch
            ├── analysis.patch
            └── *.dat
```

Branches organize code development. Experiment names organize scientific work;
one experiment can have many runs, including repeated runs of identical code
and settings. Each run records both source and analysis revisions. A timestamp
identifies the run, not the revision.

## Setup

Requirements: Bash, Git, Python 3.7 or newer, and a separately built Autofocusing
executable. Python helpers use only the standard library. These scripts support
macOS and Linux. See [the build instructions](../README.md#build).

From the source checkout:

```bash
bash Scripts/setup_workspace.sh --experiment primary-microseisms --init-git
```

The default workspace is the checkout's parent. Use `--workspace /path/to/work`
for another location. The checkout need not be named `repo`; the generated
launcher points to the actual source location relative to the workspace.

Setup creates directories and missing templates. It never overwrites existing
files, copies seismic data, builds the program, creates commits or configures
GitHub. `--init-git` initializes only `analysis/`; it refuses to initialize an
analysis directory inside another Git checkout. Omit it to defer Git setup.
The workspace itself must be outside the public source checkout.

Edit `local_config.sh` to set `HINET_ROOT`, `CMT_CATALOG` and, if necessary,
`AUTOFOCUSING_BIN`. Edit `analysis/primary-microseisms/config.sh` to set analysis
conditions. Commit analysis definitions independently when ready.

**The template retains the nominal 0.1–0.25 Hz band and ±0.165 s/km px/py grid.**
Naming an experiment `primary-microseisms` does not change its frequency band.
The intended 0.05–0.1 Hz analysis needs a separate change to the C++ frequency
handling. `AUTOFOCUSING_SLOWNESS_MAX=0.4` already widens each of px and py to
±0.4 s/km; it does not impose a radial 0.4 s/km limit.

## Run and settings

```bash
bash Scripts/run.sh --experiment primary-microseisms
# Or, using the new generated launcher:
bash ../run.sh --experiment primary-microseisms
```

Use `--workspace PATH` on either launcher to select another workspace. Both
launchers can be called from any current directory. A relative `--workspace`
is interpreted from the caller's directory. Other relative file paths are
resolved against the selected workspace, including those in experiment configs.
The executable runs with the source checkout as its working directory so it can
find `vel_Nishida2008.dat`.

Settings are applied in this order, with later assignments taking precedence:

1. Inherited environment.
2. Workspace `local_config.sh`, if present.
3. `analysis/<experiment>/config.sh` when `--experiment` is supplied (required).
4. Defaults for settings not supplied.

Configuration files are trusted Bash scripts. Ordinary assignments are exported
to the executable too. Prefer simple assignments and keep these files unchanged
while a run starts; their snapshots are taken after they are sourced. The runner
does not archive additional files sourced by those scripts. Store machine paths,
thread limits and device preferences in `local_config.sh`, and scientific
conditions in the experiment config. Experiment settings may also override
backend choices for controlled comparisons.

`--experiment` fixes the experiment name regardless of `PARAM_ID`. Names must
start with an ASCII letter or digit and contain only letters, digits, dots,
underscores or hyphens; `..` is forbidden. With no `--experiment`, the runner
uses only local/environment settings and `PARAM_ID` (default `tilt_horizontal`),
retaining the old no-argument entry point. Every new run uses a timestamp folder.

Defaults are horizontal components, start year 2004, CPU backend, slowness step
0.005 s/km and maximum 0.165 s/km. Inputs must be supplied. The executable defaults
to `bin/cal_ccf_clang`, then `bin/cal_ccf_gcc`, in the source checkout.
`RESULTS_ROOT` defaults to workspace `results/`. The existing C++ daily scan and
frequency band are unchanged; there is no new end-date or frequency option.

## What each run records

A run directory is allocated atomically, so parallel or repeated runs cannot
reuse it. UTC timestamps use `YYYYMMDDTHHMMSSZ`, with a numeric suffix on collision.
Input/configuration preflight errors return before allocating a run directory.

`manifest.json` records the run ID, timestamps, command, working directory,
allowlisted effective settings, exit code and status (`running`, `succeeded`,
`failed` or `interrupted`). It also records:

- Source and analysis repository paths, full commit SHAs, branches and dirty
  state. Detached HEAD has no branch; an unborn repository has no commit.
- Tracked differences from HEAD in `source.patch` and `analysis.patch`, status
  text files, and JSON lists of untracked files. Repositories that are not
  independently initialized are explicitly identified as such.
- SHA-256 hashes and paths of the binary, catalog and available velocity model.
  Waveforms are identified by their archive root, not by hashing the archive.
- Copies of configuration files used, and `effective_config.json` containing the
  settings passed by the launcher. These files contain local paths; review them
  before sharing results. Unrelated environment variables are not recorded.

The combined stdout/stderr stream is saved to `run.log` and shown on the terminal.
The runner forwards SIGINT/SIGTERM to the child process group, saves its final
status and returns its exit code (128 + signal for interruption). SIGKILL or
power loss can leave `status: running`; inspect the log and process state before
interpreting such a record as a live run. Results are retained on failure.

Recorded source SHA describes the checkout, **not a verified build origin**.
An existing binary may have been built from different source. The manifest
therefore sets `binary.source_revision_verified` to false and records its hash.
For formal analysis, commit source and analysis changes and rebuild the selected
binary before running. Untracked files are listed but their contents are not
archived, except the selected config snapshots; dirty runs are not guaranteed
reconstructible. Freeze/version the external waveform archive separately when
exact input reproducibility is required.

The C++ CLI is unchanged. The launcher passes the run ID in its existing
`<git-version>` directory argument and stores the actual Git revisions in the
manifest. Direct calls still use their supplied directory label and do not
create these run records; repeating a direct call can overwrite its data file.

## Existing workspaces

Setup preserves existing parent launchers, configs and results. An older parent
`run.sh` still implements the old commit-directory behavior until you migrate it.
Use `bash /path/to/repo/Scripts/run.sh --workspace /path/to/work --experiment NAME`
to select the new runner immediately. Move scientific settings from the old local
config into the experiment config as appropriate; experiment assignments win.

To migrate a customized parent launcher, review its changes and retain needed
analysis-specific steps in `analysis/`. Replace its common execution logic with
a call to `repo/Scripts/run.sh --workspace ...`, as in the new generated launcher.
Do not copy over it blindly. Newly generated launchers delegate to tracked code,
so later source updates need no launcher copy. Old output folders stay in place;
new executions are written to new timestamp folders.

Record interpretation and reusable findings in each experiment's README or
other tracked analysis documents, referencing the run IDs. GitHub synchronization
of `repo/` does not synchronize `analysis/` or `results/`. Configure an analysis
remote separately if desired, and back up large outputs independently.

## Verification and other tools

```bash
python3 -B tests/workspace.py
ctest --test-dir build -R 'workspace_scripts|parent_launcher' --output-on-failure
```

The standalone integration checks use temporary workspaces and a mock executable.
The CTest parent-launcher check also exercises the built executable with an empty
archive; neither test performs a scientific waveform analysis.

`compare_profiles.py` compares profiling logs. Global CMT conversion tools and
catalog provenance are documented in [GlobalCMT/README.md](GlobalCMT/README.md).
