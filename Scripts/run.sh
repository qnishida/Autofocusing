#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ -d "$SCRIPT_DIR/repo/src" ]]; then
  REPO_DIR="$SCRIPT_DIR/repo"
else
  # Tracked copy: repo/Scripts/run.sh. Deployed copy: ../run.sh.
  REPO_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)
fi
WORK_DIR=$(cd -- "$REPO_DIR/.." && pwd)
cd -- "$WORK_DIR"
if [[ -f "$WORK_DIR/local_config.sh" ]]; then
  source "$WORK_DIR/local_config.sh"
fi

HINET_ROOT=${HINET_ROOT:-/Volumes/Seismic_Data/hdf5/Hi-net_tilt}
COMPONENT_MODE=${COMPONENT_MODE:-horizontal}
PARAM_ID=${PARAM_ID:-tilt_horizontal}
START_YEAR=${START_YEAR:-2004}
RESULTS_ROOT=${RESULTS_ROOT:-$WORK_DIR/results}
GIT_VERSION=$(git -C "$REPO_DIR" describe --tags --always --dirty)

if [[ -z "${CMT_CATALOG:-}" || ! -f "$CMT_CATALOG" ]]; then
  echo "Set CMT_CATALOG to the converted Global CMT catalog in $WORK_DIR/local_config.sh or the environment." >&2
  exit 1
fi
if [[ ! -d "$HINET_ROOT" ]]; then
  echo "Waveform directory not found: $HINET_ROOT" >&2
  exit 1
fi
case "$COMPONENT_MODE" in
  horizontal|3c) ;;
  *) echo "COMPONENT_MODE must be horizontal or 3c" >&2; exit 1 ;;
esac
# Resolve supplied relative paths before changing to the repository for its model file.
CMT_CATALOG=$(cd -- "$(dirname -- "$CMT_CATALOG")" && pwd)/$(basename -- "$CMT_CATALOG")
HINET_ROOT=$(cd -- "$HINET_ROOT" && pwd)
[[ "$RESULTS_ROOT" = /* ]] || RESULTS_ROOT="$WORK_DIR/$RESULTS_ROOT"
if [[ -z "${AUTOFOCUSING_BIN:-}" ]]; then
  for candidate in "$REPO_DIR/bin/cal_ccf_clang" "$REPO_DIR/bin/cal_ccf_gcc"; do
    if [[ -x "$candidate" ]]; then
      AUTOFOCUSING_BIN=$candidate
      break
    fi
  done
fi
if [[ -z "${AUTOFOCUSING_BIN:-}" || ! -x "$AUTOFOCUSING_BIN" ]]; then
  echo "Build and install cal_ccf first (see repo/manual.md)." >&2
  exit 1
fi
[[ "$AUTOFOCUSING_BIN" = /* ]] || AUTOFOCUSING_BIN="$WORK_DIR/$AUTOFOCUSING_BIN"
printf 'Mode: %s; start year: %s; output: %s/%s/%s/\n' \
  "$COMPONENT_MODE" "$START_YEAR" "$RESULTS_ROOT" "$PARAM_ID" "$GIT_VERSION"
cd -- "$REPO_DIR"
exec "$AUTOFOCUSING_BIN" "$START_YEAR" "$PARAM_ID" "$GIT_VERSION" \
  "$HINET_ROOT" "$CMT_CATALOG" "$RESULTS_ROOT" "$COMPONENT_MODE"
