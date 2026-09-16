# Analysis settings override local_config.sh and inherited environment values.
export COMPONENT_MODE="horizontal"
# Inclusive dates; choose the intended interval before running.
export START_DATE=2004-01-01
export END_DATE=2004-01-07
export AUTOFOCUSING_SLOWNESS_STEP=0.005
export AUTOFOCUSING_SLOWNESS_MAX=0.165
# MAX sets each px/py grid extent in s/km (rounded down to a STEP multiple).
# Peak candidates are also masked to radial slowness <= that grid extent.
export AUTOFOCUSING_FREQ_MIN=0.1
export AUTOFOCUSING_FREQ_MAX=0.25
# Hz; for primary microseisms set MIN=0.05 and MAX=0.1 explicitly.
# Both endpoints round down to FFT bins; requested/effective bands are recorded.
