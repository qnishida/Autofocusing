"""Compare pre-optimization CPU and current GPU fixed-seed event drivers.

Run from a checkout containing vel_Nishida2008.dat. Only use deterministic
drivers produced by build_io_probe.py; production random seeds are unsuitable.
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
import time


def compare(reference, current):
    if len(reference) != 20 or len(current) != 20:
        raise RuntimeError('Expected 20 events for the qualified three-day input')
    errors = {'beam': 0.0, 'bootstrap': 0.0}
    for index, (before, after) in enumerate(zip(reference, current)):
        if len(before) != 38 or len(after) != 38:
            raise RuntimeError('Expected 38 columns')
        for column, (a, b) in enumerate(zip(before, after), 1):
            category = ('beam' if column in (12, 13) else
                        'bootstrap' if 21 <= column <= 33 else None)
            if category or column == 15:
                x, y = float(a), float(b)
                rtol = {'beam': 1e-4, 'bootstrap': 1e-3}.get(category, 0)
                atol = 1e-12 if column == 15 else 1e-30
                ok = (math.isfinite(x) and math.isfinite(y)
                      and math.isclose(x, y, rel_tol=rtol, abs_tol=atol))
                if category:
                    errors[category] = max(errors[category], abs(x-y)/max(abs(x), 1e-30))
            else:
                ok = a == b
            if not ok:
                raise RuntimeError(f'Event {index}, column {column}: {a} != {b}')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('legacy', 'current', 'hinet', 'catalog', 'output'):
        parser.add_argument(key, type=Path)
    parser.add_argument('--backend', choices=('metal', 'cuda'), required=True)
    parser.add_argument('--power', choices=('off', 'bootstrap', 'grid', 'all'), default='bootstrap')
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if min(args.threads, args.repeats) < 1:
        parser.error('positive threads/repeats required')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    for day in range(1, 4):
        source = args.hinet.resolve()/'2004'/f'010{day}'/f'200400{day}0000.h5'
        if not source.is_file():
            raise RuntimeError(f'Missing {source}')
        target = out/'input/2004'/source.parent.name/source.name
        target.parent.mkdir(parents=True)
        target.symlink_to(source)
    report = dict(platform=platform.platform(), backend=args.backend, power=args.power,
                  threads=args.threads, days=3, fftw='ESTIMATE | UNALIGNED',
                  seed='1837 + 104729 * replicate', profile='disabled in both',
                  protocol='one excluded warmup per binary; alternating fresh processes; OS caches not purged',
                  bounds='beam columns 12/13 rtol=1e-4; bootstrap 21..33 rtol=1e-3; residual 15 atol=1e-12; other fields exact',
                  runs=[])
    jobs = [('legacy', -1), ('current', -1)] + [
        (mode, repeat) for repeat in range(args.repeats)
        for mode in (('legacy', 'current') if repeat % 2 == 0 else ('current', 'legacy'))]
    reference = windows = choices = None
    event_hashes = {}
    for mode, repeat in jobs:
        tag = f'{mode}-{repeat}'
        exe = getattr(args, mode).resolve()
        env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE',
                   AUTOFOCUSING_BACKEND='cpu' if mode == 'legacy' else args.backend,
                   AUTOFOCUSING_GPU_POWER='off' if mode == 'legacy' else args.power,
                   AUTOFOCUSING_METAL_POWER='off' if mode == 'legacy' else args.power)
        for key in ('AUTOFOCUSING_PROFILE', 'AUTOFOCUSING_VERIFY_METAL',
                    'AUTOFOCUSING_VERIFY_CUDA', 'AUTOFOCUSING_CUDA_PROFILE',
                    'AUTOFOCUSING_SLOWNESS_STEP', 'AUTOFOCUSING_SLOWNESS_MAX'):
            env.pop(key, None)
        timing = (['/usr/bin/time', '-l'] if platform.system() == 'Darwin'
                  else ['/usr/bin/time', '-f', '#RSS_KIB %M'])
        print(f'Running {tag}', flush=True)
        start = time.monotonic()
        log_path = out/f'{tag}.log'
        with log_path.open('w') as log:
            subprocess.run(timing+[str(exe), '2004', 'legacy-total', tag,
                                   str(out/'input'), str(args.catalog.resolve()),
                                   str(out/'results'), 'horizontal'], cwd=Path.cwd(),
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=900)
        elapsed = time.monotonic()-start
        files = list((out/'results/legacy-total'/tag).glob('*.dat'))
        if len(files) != 1:
            raise RuntimeError('Expected one event file')
        contents = files[0].read_bytes()
        rows = [line.split() for line in contents.decode().splitlines()]
        lines = log_path.read_text().splitlines()
        run_windows = [line for line in lines if line.startswith('#Segment ')]
        run_choices = [line for line in lines if line.startswith(('#deg max=', '#dp_Δ max='))]
        if len(run_windows) != 12 or len(run_choices) != 40:
            raise RuntimeError('Missing windows or initial candidate records')
        if reference is None:
            reference, windows, choices = rows, run_windows, run_choices
        errors = compare(reference, rows)
        if run_windows != windows or run_choices != choices:
            raise RuntimeError('Accepted windows or initial candidates changed')
        event_hash = hashlib.sha256(contents).hexdigest()
        if mode in event_hashes and event_hash != event_hashes[mode]:
            raise RuntimeError(f'{mode} output is not repeatable')
        event_hashes[mode] = event_hash
        rss = next((int(line.split()[0]) for line in lines
                    if 'maximum resident set size' in line), None)
        if rss is None:
            rss = next((int(line.split()[1])*1024 for line in lines
                        if line.startswith('#RSS_KIB ')), None)
        report['runs'].append(dict(mode=mode, repeat=repeat, wall_s=elapsed,
                                   events=len(rows), max_relative=errors, windows=run_windows,
                                   choices=run_choices, max_rss_bytes=rss,
                                   binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                                   event_sha256=event_hash, event_file=str(files[0].relative_to(out))))
        (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        print(f'PASS {tag}: {elapsed:.3f}s, errors={errors}', flush=True)
    report['summary'] = {}
    for mode in ('legacy', 'current'):
        values = [r['wall_s'] for r in report['runs'] if r['mode'] == mode and r['repeat'] >= 0]
        report['summary'][mode] = dict(median_s=statistics.median(values), min_s=min(values), max_s=max(values))
    before, after = (report['summary'][mode]['median_s'] for mode in ('legacy', 'current'))
    report.update(speedup=before/after, elapsed_reduction_percent=100*(1-after/before))
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['summary']), flush=True)


if __name__ == '__main__':
    main()
