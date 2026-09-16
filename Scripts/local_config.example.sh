# Machine settings. Relative paths are resolved against the workspace.
# Set these to your external inputs; no waveform data or catalog is installed.
export HINET_ROOT="/path/to/hdf5/Hi-net_tilt"
export CMT_CATALOG="moment_loc_76_24"
# Optional: export AUTOFOCUSING_BIN="repo/bin/cal_ccf_clang"
# Optional: export RESULTS_ROOT="results"
# Optional: export OMP_NUM_THREADS=4
# Optional: export AUTOFOCUSING_BACKEND="metal"  # cpu (default), metal, cuda
# Optional: export AUTOFOCUSING_GPU_POWER="bootstrap" # off, bootstrap, grid, all
# AUTOFOCUSING_METAL_POWER is supported when AUTOFOCUSING_GPU_POWER is unset.
