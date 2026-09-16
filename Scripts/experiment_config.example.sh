# Analysis settings override local_config.sh and inherited environment values.
export COMPONENT_MODE="horizontal"
export START_YEAR=2004
export AUTOFOCUSING_SLOWNESS_STEP=0.005
export AUTOFOCUSING_SLOWNESS_MAX=0.165
# MAX limits each of px and py in s/km, not the radial slowness.
# Frequency is currently fixed in C++ at nominal 0.1-0.25 Hz.
# Naming an experiment primary-microseisms does NOT change its frequency band.
