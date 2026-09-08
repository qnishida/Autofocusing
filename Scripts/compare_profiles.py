#!/usr/bin/env python3
"""Compare AUTOFOCUSING_PROFILE logs from identical input and thread settings."""
import re
import sys
from pathlib import Path

KEYS = ('fft_qc_s', 'pack_s', 'stack_s', 'fit_s', 'total_s')

def read(path):
    lines = Path(path).read_text().splitlines()
    profiles = [dict(re.findall(r'(\w+)=([\d.e+-]+)', line))
                for line in lines if line.startswith('#PROFILE segment=')]
    if not profiles:
        raise SystemExit(f'No segment profiles found: {path}')
    return profiles

if len(sys.argv) != 3:
    raise SystemExit('Usage: compare_profiles.py baseline.log optimized.log')
baseline, optimized = (read(path) for path in sys.argv[1:])
if [(p['segment'], p['windows']) for p in baseline] != [
    (p['segment'], p['windows']) for p in optimized
]:
    raise SystemExit('Segment/window counts differ; inspect results before comparing speed.')
print('stage baseline_seconds optimized_seconds speedup')
for key in KEYS:
    before = sum(float(p[key]) for p in baseline)
    after = sum(float(p[key]) for p in optimized)
    ratio = f'{before/after:.3f}' if after else 'n/a'
    print(f'{key} {before:.6f} {after:.6f} {ratio}')
print('accepted_windows', sum(int(p['windows']) for p in optimized))
