"""Compare deterministic pre/post I/O drivers on January 1–3, 2004."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('reference', type=Path)
p.add_argument('current', type=Path)
p.add_argument('hinet', type=Path)
p.add_argument('catalog', type=Path)
p.add_argument('output', type=Path, help='new directory for bounded input, logs and report')
p.add_argument('--backend', choices=('cpu', 'metal'), default='cpu')
p.add_argument('--repeats', type=int, default=2)
a = p.parse_args()
if a.repeats < 1:
    p.error('positive repeat count required')
out = a.output.resolve()
out.mkdir(parents=True, exist_ok=False)
for day in range(1, 4):
    source = a.hinet.resolve() / '2004' / f'010{day}' / f'200400{day}0000.h5'
    if not source.is_file():
        raise RuntimeError(f'Missing {source}')
    target = out / 'input' / '2004' / f'010{day}'
    target.mkdir(parents=True)
    (target / source.name).symlink_to(source)
report = {'platform': platform.platform(), 'backend': a.backend, 'threads': 16,
          'days': ['2004-01-01', '2004-01-02', '2004-01-03'],
          'bootstrap': '1837 + 104729 * replicate', 'fftw': 'ESTIMATE | UNALIGNED',
          'cache': 'not purged; alternating process order', 'runs': []}
expected = None
profiles = None
for trial in range(a.repeats):
    for mode in (('reference', 'current') if trial % 2 == 0 else ('current', 'reference')):
        exe = getattr(a, mode).resolve()
        tag = f'{mode}-{trial}'
        env = dict(os.environ, OMP_NUM_THREADS='16', AUTOFOCUSING_BACKEND=a.backend,
                   AUTOFOCUSING_PROFILE='1')
        for key in ('AUTOFOCUSING_SLOWNESS_STEP', 'AUTOFOCUSING_SLOWNESS_MAX', 'AUTOFOCUSING_VERIFY_METAL'):
            env.pop(key, None)
        command = [str(exe), '2004', 'io', tag, str(out / 'input'), str(a.catalog.resolve()),
                   str(out / 'results'), 'horizontal']
        print(f'Running {a.backend} {tag}', flush=True)
        start = time.monotonic()
        log = out / f'{tag}.log'
        with log.open('w') as stream:
            subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=900)
        wall = time.monotonic() - start
        files = list((out / 'results' / 'io' / tag).glob('*.dat'))
        if len(files) != 1:
            raise RuntimeError('Expected one result file')
        events = files[0].read_bytes()
        if len(events.splitlines()) != 20:
            raise RuntimeError('Expected 20 events')
        segments = [dict(x.split('=', 1) for x in line.split()[1:]) for line in log.read_text().splitlines()
                    if line.startswith('#PROFILE segment=')]
        if len(segments) != 12:
            raise RuntimeError('Expected 12 segments')
        windows = [(s['segment'], s['windows']) for s in segments]
        if expected is None:
            expected, profiles = events, windows
        if events != expected or windows != profiles:
            raise RuntimeError(f'Event bytes or accepted windows changed: {tag}')
        report['runs'].append({'mode': mode, 'trial': trial, 'wall_s': wall, 'events': 20,
                               'event_sha256': hashlib.sha256(events).hexdigest(),
                               'binary_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
                               'segments': segments})
        (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        print(f'PASS {tag}: 20 events byte-identical, {wall:.3f}s', flush=True)
report['median_s'] = {mode: statistics.median(r['wall_s'] for r in report['runs'] if r['mode'] == mode)
                      for mode in ('reference', 'current')}
report['speedup'] = report['median_s']['reference'] / report['median_s']['current']
report['within_5_percent_regression'] = report['median_s']['current'] <= 1.05 * report['median_s']['reference']
(out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
if not report['within_5_percent_regression']:
    raise RuntimeError('End-to-end regression exceeds 5%')
print(json.dumps(report['median_s']), flush=True)
