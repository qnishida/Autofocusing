#!/usr/bin/env bash
set -euo pipefail
usage() {
  printf 'Usage: bash %s --experiment NAME [--workspace PATH]\n' "$0"
  printf 'Example: bash %q --experiment primary-microseisms\n' "$0"
}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ -d "$SCRIPT_DIR/repo/src" ]]; then
  REPO_DIR="$SCRIPT_DIR/repo" # Compatibility with a copied parent launcher.
else
  REPO_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)
fi
WORK_DIR=$(dirname -- "$REPO_DIR")
EXPERIMENT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --workspace|--experiment)
      [[ $# -ge 2 && -n "$2" ]] || { echo "Missing value for $1" >&2; usage >&2; exit 2; }
      if [[ "$1" == --workspace ]]; then WORK_DIR=$2; else EXPERIMENT=$2; fi
      shift 2 ;;
    -h|--help)
      usage
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
if [[ -z "$EXPERIMENT" ]]; then
  echo "--experiment is required; no analysis was started." >&2
  usage >&2
  exit 2
fi
if [[ ! "$EXPERIMENT" =~ ^[a-zA-Z0-9][a-zA-Z0-9._-]*$ || "$EXPERIMENT" == *..* ]]; then
  echo "Invalid experiment name: $EXPERIMENT" >&2; exit 2
fi
WORK_DIR=$(cd -- "$WORK_DIR" && pwd)
readonly REPO_DIR WORK_DIR EXPERIMENT
cd -- "$WORK_DIR"
set -a
if [[ -f "$WORK_DIR/local_config.sh" ]]; then
  source "$WORK_DIR/local_config.sh"
fi
if [[ ! -f "$WORK_DIR/analysis/$EXPERIMENT/config.sh" ]]; then
  echo "Experiment config not found: analysis/$EXPERIMENT/config.sh" >&2; exit 1
fi
source "$WORK_DIR/analysis/$EXPERIMENT/config.sh"
set +a
cd -- "$WORK_DIR"
args=(run --repo "$REPO_DIR" --workspace "$WORK_DIR" --experiment "$EXPERIMENT")
exec python3 "$REPO_DIR/Scripts/workspace.py" "${args[@]}"
