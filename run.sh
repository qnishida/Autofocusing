#!/usr/bin/env bash

PARAM_ID="param_set_A"
GIT_VERSION=$(git describe --tags --always)

if [[ -f local_config.sh ]]; then
  # Local machine paths. This file is intentionally not tracked.
  source local_config.sh
fi

if [[ -z "${HINET_ROOT:-}" || -z "${CMT_CATALOG:-}" ]]; then
  echo "Usage: HINET_ROOT=/path/to/hdf5 CMT_CATALOG=/path/to/moment_loc_76_24 ./run.sh" >&2
  exit 1
fi

OUTDIR="output/${PARAM_ID}/${GIT_VERSION}"
mkdir -p "$OUTDIR"

./bin/cal_ccf_gcc 2004 "$PARAM_ID" "$GIT_VERSION" "$HINET_ROOT" "$CMT_CATALOG"
