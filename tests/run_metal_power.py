"""Paired real-event qualification/timings for Metal or CUDA objectives.

The historical filename and default Metal backend are retained for compatibility.
Use a deterministic driver from build_io_probe.py, never a production binary.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time


def compare_power(rows, reference, tag, bootstrap=True):
    if not rows or len(rows) != len(reference):
        raise RuntimeError(f'{tag}: empty result or event count mismatch')
    max_relative = 0.
    for row, ref in zip(rows, reference):
        if len(row) != 38 or len(ref) != 38:
            raise RuntimeError('Expected 38 columns')
        for col, (value, expected) in enumerate(zip(row, ref), 1):
            if bootstrap and 21 <= col <= 33:
                x, y = float(value), float(expected)
                if not math.isfinite(x) or not math.isclose(x, y, rel_tol=1e-3, abs_tol=1e-30):
                    raise RuntimeError(f'{tag} column {col}: {x} != {y}')
                max_relative = max(max_relative, abs(x-y)/max(abs(y), 1e-30))
            elif col == 15:
                if not math.isclose(float(value), float(expected), rel_tol=0, abs_tol=1e-12):
                    raise RuntimeError(f'{tag}: convergence residual differs')
            elif value != expected:
                raise RuntimeError(f'{tag} column {col}: {value} != {expected}')
    return max_relative


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary', type=Path)
    p.add_argument('hinet', type=Path)
    p.add_argument('catalog', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--backend', choices=('metal', 'cuda'), default='metal')
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--modes', nargs='+', choices=('off', 'bootstrap', 'grid', 'all'),
                   default=['off', 'bootstrap', 'grid', 'all'])
    p.add_argument('--days', type=int, default=3, choices=(3, 5))
    p.add_argument('--repeats', type=int, default=5)
    p.add_argument('--compare-cpu', action='store_true', help='also pair CPU/off with GPU/off each repeat')
    p.add_argument('--expected-events', type=int, help='optional independently known count; never assumes Metal counts on CUDA')
    p.add_argument('--legacy-run', type=Path, help='reuse a completed legacy-0 run from this script')
    p.add_argument('--legacy', type=Path, help='optional frozen driver with this GPU backend')
    a = p.parse_args()
    if a.threads < 1 or a.repeats < 1 or len(set(a.modes)) != len(a.modes):
        p.error('positive threads/repeats and distinct modes required')
    if 'off' not in a.modes and not (a.legacy or a.legacy_run):
        p.error('include off as the power baseline, or supply a legacy baseline')
    if a.compare_cpu and 'off' not in a.modes:
        p.error('--compare-cpu requires off')
    if not a.binary.is_file() or not a.catalog.is_file():
        p.error('binary and catalog must exist')
    sources = [a.hinet.resolve()/'2004'/f'010{day}'/f'200400{day}0000.h5'
               for day in range(1, a.days+1)]
    for source in sources:
        if not source.is_file():
            raise RuntimeError(f'Missing {source}')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    for source in sources:
        target = out/'input/2004'/source.parent.name
        target.mkdir(parents=True)
        (target/source.name).symlink_to(source)
    report = {'platform': platform.platform(), 'backend': a.backend, 'threads': a.threads,
              'days': a.days, 'fftw': 'ESTIMATE | UNALIGNED',
              'seed': '1837 + 104729 * replicate',
              'baseline': 'legacy' if a.legacy or a.legacy_run else f'{a.backend}/off',
              'cache': 'not purged; alternate process order; includes initialization/packing/transfers/sync',
              'runs': []}
    reference = reference_choices = None
    if a.legacy_run:
        files = list((a.legacy_run/'results/power/legacy-0').glob('*.dat'))
        if len(files) != 1:
            raise RuntimeError('Expected legacy result')
        reference = [line.split() for line in files[0].read_text().splitlines()]
        reference_choices = [line for line in (a.legacy_run/'legacy-0.log').read_text().splitlines()
                             if line.startswith(('#deg max=', '#dp_Δ max='))]
        report['legacy_event_sha256'] = hashlib.sha256(files[0].read_bytes()).hexdigest()
    # off must establish a baseline even when the CLI lists it last.
    modes = (['off'] if 'off' in a.modes else []) + [m for m in a.modes if m != 'off']
    paired = (['cpu'] if a.compare_cpu else []) + modes
    jobs = ([('legacy', 0)] if a.legacy else []) + [
        (mode, r) for r in range(a.repeats) for mode in (paired if r % 2 == 0 else paired[::-1])]
    event_files = {}
    for mode, repeat in jobs:
        tag = f'{mode}-{repeat}'
        exe = (a.legacy if mode == 'legacy' else a.binary).resolve()
        backend = 'cpu' if mode == 'cpu' else a.backend
        power = 'off' if mode in ('legacy', 'cpu') else mode
        env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), AUTOFOCUSING_BACKEND=backend,
                   AUTOFOCUSING_GPU_POWER=power, AUTOFOCUSING_METAL_POWER=power, AUTOFOCUSING_PROFILE='1')
        for key in ('AUTOFOCUSING_VERIFY_METAL', 'AUTOFOCUSING_VERIFY_CUDA', 'AUTOFOCUSING_CUDA_PROFILE',
                    'AUTOFOCUSING_SLOWNESS_STEP', 'AUTOFOCUSING_SLOWNESS_MAX'):
            env.pop(key, None)
        print(f'Running {backend} {tag}', flush=True)
        start = time.monotonic()
        log = out/f'{tag}.log'
        with log.open('w') as stream:
            timing = (['/usr/bin/time', '-l'] if platform.system() == 'Darwin' else
                      ['/usr/bin/time', '-f', '#RSS_KIB %M'])
            subprocess.run(timing+[str(exe), '2004', 'power', tag, str(out/'input'), str(a.catalog.resolve()),
                                   str(out/'results'), 'horizontal'], env=env, stdout=stream,
                           stderr=subprocess.STDOUT, check=True, timeout=900)
        wall = time.monotonic()-start
        files = list((out/'results/power'/tag).glob('*.dat'))
        if len(files) != 1:
            raise RuntimeError('Expected one event output')
        event_files[(mode, repeat)] = files[0]
        contents = files[0].read_text()
        rows = [line.split() for line in contents.splitlines()]
        if not rows or (a.expected_events is not None and len(rows) != a.expected_events):
            raise RuntimeError(f'Unexpected event count: {len(rows)}')
        lines = log.read_text().splitlines()
        choices = [line for line in lines if line.startswith(('#deg max=', '#dp_Δ max='))]
        max_rel = 0.
        if mode != 'cpu':
            if reference is None:
                reference, reference_choices = rows, choices
            if choices != reference_choices:
                raise RuntimeError(f'{tag}: initial candidates differ')
            max_rel = compare_power(rows, reference, tag, mode in ('bootstrap', 'all'))
        if a.compare_cpu and ('cpu', repeat) in event_files and ('off', repeat) in event_files and mode in ('cpu', 'off'):
            subprocess.run([sys.executable, str(Path(__file__).with_name('compare_event_results.py')),
                            str(event_files[('cpu', repeat)]), str(event_files[('off', repeat)]), '--gpu'], check=True)
        stage, rss, segments = {}, None, []
        for line in lines:
            if 'maximum resident set size' in line:
                rss = int(line.split()[0])
            if line.startswith('#RSS_KIB '):
                rss = int(line.split()[1])*1024
            if line.startswith('#POWER_PROFILE '):
                fields = dict(x.split('=', 1) for x in line.split()[1:])
                stage[fields['stage']] = stage.get(fields['stage'], 0)+float(fields['total_s'])
            if line.startswith('#PROFILE segment='):
                segments.append(dict(x.split('=', 1) for x in line.split()[1:]))
        if not segments or not all(s in stage for s in ('grid', 'bootstrap')):
            raise RuntimeError('Missing segment/power profiles')
        windows = [(s['segment'], s['windows']) for s in segments]
        if report['runs'] and windows != report['runs'][0]['windows']:
            raise RuntimeError('Accepted segments/windows differ')
        report['runs'].append({'mode': mode, 'repeat': repeat, 'wall_s': wall, 'stages': stage,
                              'events': len(rows), 'max_rss_bytes': rss, 'windows': windows,
                              'stack_s': sum(float(s['stack_s']) for s in segments),
                              'choices': choices, 'max_bootstrap_relative': max_rel,
                              'binary_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
                              'event_sha256': hashlib.sha256(contents.encode()).hexdigest()})
        (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        print(f'PASS {tag}: {len(rows)} events, max_boot_rel={max_rel:.3g}, wall={wall:.3f}s', flush=True)
    summary = {}
    for mode in paired:
        runs = [r for r in report['runs'] if r['mode'] == mode]
        summary[mode] = {'median_s': statistics.median(r['wall_s'] for r in runs),
                         'min_s': min(r['wall_s'] for r in runs), 'max_s': max(r['wall_s'] for r in runs),
                         'stack_s': statistics.median(r['stack_s'] for r in runs),
                         'stages': {s: statistics.median(r['stages'][s] for r in runs) for s in ('grid', 'bootstrap')}}
    for mode in (m for m in modes if m != 'off' and 'off' in summary):
        stages = ('grid', 'bootstrap') if mode == 'all' else (mode,)
        before = sum(summary['off']['stages'][s] for s in stages)
        after = sum(summary[mode]['stages'][s] for s in stages)
        summary[mode].update(stage_speedup=before/after, stage_pass=after <= .9*before,
                             end_to_end_pass=summary[mode]['median_s'] <= 1.05*summary['off']['median_s'])
    if a.compare_cpu:
        summary['off'].update(cpu_speedup=summary['cpu']['median_s']/summary['off']['median_s'],
                              stack_speedup=summary['cpu']['stack_s']/summary['off']['stack_s'],
                              stage_pass=summary['off']['stack_s'] <= .9*summary['cpu']['stack_s'],
                              end_to_end_pass=summary['off']['median_s'] <= 1.05*summary['cpu']['median_s'])
    report['summary'] = summary
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
