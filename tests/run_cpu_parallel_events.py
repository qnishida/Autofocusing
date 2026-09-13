"""Paired same-backend pre/post CPU-change checks; fixed-seed drivers only."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time

from run_metal_power import compare_power


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('reference', 'current', 'hinet', 'catalog', 'output'):
        p.add_argument(name, type=Path)
    p.add_argument('--backend', choices=('cpu', 'cuda', 'metal'), default='cuda')
    p.add_argument('--power', choices=('off', 'bootstrap', 'grid', 'all'), default='bootstrap')
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--days', type=int, choices=(3, 5), default=3)
    p.add_argument('--repeats', type=int, default=5)
    p.add_argument('--warmups', type=int, choices=(0, 1), default=1,
                   help='use 0 for correctness-only qualification, 1 for timing')
    p.add_argument('--target', choices=('hessian', 'objective', 'rotation'))
    a = p.parse_args()
    if min(a.threads, a.repeats) < 1:
        p.error('positive threads and repeats required')
    for path in (a.reference, a.current, a.catalog):
        if not path.is_file():
            p.error(f'missing file: {path}')
    sources = [a.hinet.resolve()/'2004'/f'010{day}'/f'200400{day}0000.h5'
               for day in range(1, a.days+1)]
    if not all(path.is_file() for path in sources):
        p.error('missing daily input')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    for path in sources:
        link = out/'input'/'2004'/path.parent.name/path.name
        link.parent.mkdir(parents=True)
        link.symlink_to(path)
    report = dict(platform=platform.platform(), backend=a.backend, power=a.power,
                  threads=a.threads, days=a.days, target=a.target,
                  cache=f'{a.warmups} excluded warmup per binary; OS cache not purged; alternating fresh processes',
                  timing='inclusive function time; outer-parallel time is worker-call sum, not wall time',
                  fftw='ESTIMATE | UNALIGNED', seed='1837 + 104729 * replicate', runs=[])
    reference_rows = reference_windows = reference_choices = None
    jobs = ([('reference', -1), ('current', -1)] if a.warmups else []) + [
        (mode, repeat) for repeat in range(a.repeats)
        for mode in (('reference', 'current') if repeat % 2 == 0 else ('current', 'reference'))]
    for mode, repeat in jobs:
        exe = getattr(a, mode).resolve()
        tag = f'{mode}-{repeat}'
        env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), OMP_DYNAMIC='FALSE',
                   AUTOFOCUSING_BACKEND=a.backend, AUTOFOCUSING_GPU_POWER=a.power,
                   AUTOFOCUSING_METAL_POWER=a.power, AUTOFOCUSING_PROFILE='1')
        for key in ('AUTOFOCUSING_VERIFY_METAL', 'AUTOFOCUSING_VERIFY_CUDA',
                    'AUTOFOCUSING_CUDA_PROFILE', 'AUTOFOCUSING_SLOWNESS_STEP',
                    'AUTOFOCUSING_SLOWNESS_MAX'):
            env.pop(key, None)
        timing = (['/usr/bin/time', '-l'] if platform.system() == 'Darwin'
                  else ['/usr/bin/time', '-f', '#RSS_KIB %M'])
        start = time.monotonic()
        with (out/f'{tag}.log').open('w') as log:
            subprocess.run(timing+[str(exe), '2004', 'cpu-parallel', tag,
                                   str(out/'input'), str(a.catalog.resolve()),
                                   str(out/'results'), 'horizontal'],
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=900)
        wall = time.monotonic()-start
        files = list((out/'results'/'cpu-parallel'/tag).glob('*.dat'))
        if len(files) != 1:
            raise RuntimeError('Expected one event file')
        contents = files[0].read_bytes()
        rows = [line.split() for line in contents.decode().splitlines()]
        if len(rows) != (20 if a.days == 3 else 29):
            raise RuntimeError('Unexpected event count for qualified input')
        lines = (out/f'{tag}.log').read_text().splitlines()
        choices = [line for line in lines if line.startswith(('#deg max=', '#dp_Δ max='))]
        segments, loads, calls, powers, rss = [], [], {}, {}, None
        for line in lines:
            if line.startswith(('#PROFILE segment=', '#LOAD_PROFILE ', '#CPU_PROFILE ', '#POWER_PROFILE ')):
                fields = dict(x.split('=', 1) for x in line.split()[1:])
                if line.startswith('#PROFILE '):
                    segments.append(fields)
                elif line.startswith('#LOAD_PROFILE '):
                    loads.append(fields)
                elif line.startswith('#CPU_PROFILE '):
                    calls[fields['stage']+'/'+fields['context']] = {
                        'calls': int(fields['calls']), 'inclusive_s': float(fields['inclusive_s'])}
                else:
                    powers[fields['stage']] = powers.get(fields['stage'], 0)+float(fields['total_s'])
            elif line.startswith('#RSS_KIB '):
                rss = int(line.split()[1])*1024
            elif 'maximum resident set size' in line:
                rss = int(line.split()[0])
        if len(segments) != 4*a.days or len(loads) != a.days or len(calls) != 8:
            raise RuntimeError('Missing expected CPU/load/segment profiling')
        windows = [(s['segment'], s['windows']) for s in segments]
        if reference_rows is None:
            reference_rows, reference_windows, reference_choices = rows, windows, choices
        # No FP32 relaxation: this is a same-backend CPU implementation change.
        compare_power(rows, reference_rows, tag, bootstrap=False)
        if windows != reference_windows or choices != reference_choices:
            raise RuntimeError('Accepted windows or grid candidates changed')
        record = dict(mode=mode, repeat=repeat, wall_s=wall, events=len(rows),
                      event_sha256=hashlib.sha256(contents).hexdigest(),
                      binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                      cpu=calls, power=powers, max_rss_bytes=rss, windows=windows,
                      segment_s={key: sum(float(s[key]) for s in segments)
                                 for key in ('fft_qc_s', 'pack_s', 'stack_s', 'fit_s', 'total_s')},
                      load_s={key: sum(float(s[key]) for s in loads)
                              for key in ('init_s', 'read_decode_copy_s', 'filter_s', 'total_s')})
        report['runs'].append(record)
        (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        print(f'PASS {tag}: {len(rows)} events, wall={wall:.3f}s', flush=True)
    summary = {}
    for mode in ('reference', 'current'):
        runs = [r for r in report['runs'] if r['mode'] == mode and r['repeat'] >= 0]
        summary[mode] = dict(median_s=statistics.median(r['wall_s'] for r in runs),
                             min_s=min(r['wall_s'] for r in runs), max_s=max(r['wall_s'] for r in runs),
                             cpu={key: statistics.median(r['cpu'][key]['inclusive_s'] for r in runs)
                                  for key in runs[0]['cpu']},
                             power={key: statistics.median(r['power'][key] for r in runs)
                                    for key in runs[0]['power']})
    report['summary'] = summary
    report['whole_run_pass'] = summary['current']['median_s'] <= 1.05*summary['reference']['median_s']
    if a.target:
        key = a.target+'/serial_caller'
        report['target_stage_pass'] = summary['current']['cpu'][key] <= .9*summary['reference']['cpu'][key]
    report['all_event_bytes_identical'] = len({r['event_sha256'] for r in report['runs']}) == 1
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
