# Paths are relative to the parent directory containing run.sh, repo/, results/.
export HINET_ROOT="/Volumes/Seismic_Data/hdf5/Hi-net_tilt"
export CMT_CATALOG="moment_loc_76_24"
export COMPONENT_MODE="horizontal"
export PARAM_ID="tilt_horizontal"
export START_YEAR="2004"
# Optional: export OMP_NUM_THREADS=4
# Optional: export RESULTS_ROOT="results"
# Apple Metal: slant stack uses float on the GPU; CPU remains the default.
# export AUTOFOCUSING_BACKEND="metal"
# export OMP_NUM_THREADS=16
# Additional horizontal fitting objectives (default off): bootstrap / grid / all.
# export AUTOFOCUSING_METAL_POWER="all"
# Finer/wider grid, in seconds/km (defaults: 0.005, 0.165).
# Use a distinct PARAM_ID for each backend/grid experiment to preserve results.
# export AUTOFOCUSING_SLOWNESS_STEP=0.0025
# export AUTOFOCUSING_SLOWNESS_MAX=0.25
# export PARAM_ID="tilt_horizontal_metal_dp0025_p025"
