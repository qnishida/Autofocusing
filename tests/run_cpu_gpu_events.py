"""Compare a GPU matrix-audit run against an archived deterministic CPU run.

Uses identical executable, input hashes, start/end dates, FFTW plans, seeds and
thread count. GPU slant-stack verification also checks every accepted window.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

from run_spectral_matrix_events import digest


def compare(cpu, gpu):
    a = [s.split() for s in cpu.read_text().splitlines()]
    b = [s.split() for s in gpu.read_text().splitlines()]
    result = dict(cpu_events=len(a), gpu_events=len(b), checks={}, columns={}, mismatches=[])
    checks = result['checks']
    checks['event_counts_equal'] = len(a) == len(b) and bool(a)
    checks['columns_valid'] = all(len(r) == 38 for r in a + b)
    if not all(checks.values()):
        return result
    checks['all_output_finite'] = all(math.isfinite(float(v)) for r in a + b for i, v in enumerate(r) if i not in (0, 18))
    checks['event_order_and_times_equal'] = [(r[0], r[1], r[18]) for r in a] == [(r[0], r[1], r[18]) for r in b]
    checks['fields_within_declared_tolerances'] = True
    for col in range(38):
        values = []
        changed = 0
        for n, (x, y) in enumerate(zip(a, b), 1):
            before, after = x[col], y[col]
            changed += before != after
            ok = before == after
            if col not in (0, 18):
                u, v = float(before), float(after)
                values.append(abs(v - u) / max(abs(u), abs(v), 1e-30))
                if col in (11, 12):
                    ok = math.isclose(u, v, rel_tol=1e-4, abs_tol=1e-30)
                elif col == 14:
                    ok = math.isclose(u, v, rel_tol=0, abs_tol=1e-12)
            if not ok:
                checks['fields_within_declared_tolerances'] = False
                if len(result['mismatches']) < 30:
                    result['mismatches'].append(dict(row=n, column=col + 1, cpu=before, gpu=after))
        result['columns'][str(col + 1)] = dict(changed_rows=changed, max_scaled_difference=max(values, default=None))
    checks['bootstrap_fields_exact'] = all(x[20:33] == y[20:33] for x, y in zip(a, b))
    checks['matrix_fields_exact'] = all(x[33:] == y[33:] for x, y in zip(a, b))
    checks['threshold_selection_equal'] = True
    result['threshold_counts'] = {}
    for name, rows in [('cpu', a), ('gpu', b)]:
        selection = [float(r[12]) > 0 and float(r[11]) / float(r[12]) > (35 if int(r[1]) == 2 else 7) for r in rows]
        result['threshold_counts'][name] = {str(c): sum(flag and int(r[1]) == c for flag, r in zip(selection, rows)) for c in range(3)}
        if name == 'cpu':
            cpu_selection = selection
        else:
            checks['threshold_selection_equal'] = selection == cpu_selection
    result['tolerances'] = 'Columns 12/13: rtol=1e-4, atol=1e-30; column 15: atol=1e-12; all other columns exact output strings, including bootstrap and matrices.'
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cpu_run', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--backend', choices=['metal', 'cuda'], default='metal')
    a = p.parse_args()
    repo = Path(__file__).resolve().parents[1]
    cpu = a.cpu_run.resolve()
    baseline = json.loads((cpu / 'report.json').read_text())
    if baseline['settings']['AUTOFOCUSING_BACKEND'] != 'cpu' or not all(baseline['checks'].values()):
        raise RuntimeError('Expected a passing CPU audit')
    binary = Path(baseline['command'][0])
    catalog = Path(baseline['command'][5])
    if digest(binary) != baseline['binary_sha256'] or digest(catalog) != baseline['catalog_sha256']:
        raise RuntimeError('Binary or source catalog differs from the CPU audit')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    for s in baseline['sources']:
        source = Path(s['path'])
        if digest(source) != s['sha256']:
            raise RuntimeError(f'Input differs from CPU audit: {source}')
        link = out / 'input' / source.parent.parent.name / source.parent.name / source.name
        link.parent.mkdir(parents=True)
        link.symlink_to(source)
    env = {k: v for k, v in os.environ.items() if not k.startswith(('AUTOFOCUSING_', 'OMP_', 'KMP_', 'MATRIX_'))}
    settings = dict(baseline['settings'])
    settings.update(AUTOFOCUSING_BACKEND=a.backend, AUTOFOCUSING_GPU_POWER='off')
    settings['AUTOFOCUSING_VERIFY_' + a.backend.upper()] = '1'
    env.update(settings)
    env['MATRIX_REFERENCE_OUTPUT'] = str(out / 'reference-events.dat')
    cmd = [str(binary), baseline['start'][:4], 'gpu-audit', 'corrected', str(out / 'input'),
           str(catalog), str(out / 'results'), '3c', baseline['start'], baseline['end']]
    report = dict(cpu_run=str(cpu), cpu_report_sha256=digest(cpu / 'report.json'),
                  binary_sha256=digest(binary), script_sha256=digest(Path(__file__)), settings=settings,
                  sources=baseline['sources'], catalog_sha256=digest(catalog), command=cmd, checks={})

    def save():
        (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')

    save()
    with (out / 'run.log').open('w') as log:
        process = subprocess.run(cmd, cwd=repo, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=7200)
    report['exit_code'] = process.returncode
    save()
    process.check_returncode()
    event_file, = (out / 'results/gpu-audit/corrected').glob('*.dat')
    shutil.copyfile(event_file, out / 'corrected-events.dat')
    log = (out / 'run.log').read_text()
    cpu_log = (cpu / 'run.log').read_text()
    checks = report['checks']
    checks['gpu_backend_confirmed'] = '#SlantStack backend=' + a.backend + ' (float)' in log
    segments = [s for s in log.splitlines() if s.startswith('#Segment ')]
    checks['segments_and_windows_equal'] = segments == baseline['segments']
    markers = ('#deg max=', '#dp_Δ max=')
    checks['initial_candidates_equal'] = [s for s in log.splitlines() if s.startswith(markers)] == [s for s in cpu_log.splitlines() if s.startswith(markers)]
    samples = [tuple(map(float, m)) for m in re.findall(r'#(?:Metal|Cuda)Check max_scaled_error=(\S+) relative_l2=(\S+) peak_mismatches=(\d+)', log)]
    expected = sum(int(s.split('accepted_windows=')[1]) for s in baseline['segments'])
    checks['every_accepted_window_verified'] = len(samples) == expected and bool(samples)
    report['window_verification'] = dict(count=len(samples), max_scaled_error=max((s[0] for s in samples), default=None),
                                       max_relative_l2=max((s[1] for s in samples), default=None),
                                       total_peak_mismatches=sum(int(s[2]) for s in samples))
    checks['inputs_unchanged'] = all(digest(Path(s['path'])) == s['sha256'] for s in baseline['sources'])
    checks['source_catalog_unchanged'] = digest(catalog) == baseline['catalog_sha256']
    report['comparisons'] = {}
    for label in ['reference', 'corrected']:
        if digest(cpu / (label + '-events.dat')) != baseline['event_hashes'][label]:
            raise RuntimeError('CPU catalog changed')
        comparison = compare(cpu / (label + '-events.dat'), out / (label + '-events.dat'))
        report['comparisons'][label] = comparison
        checks[label + '_comparison_passed'] = all(comparison['checks'].values())
    report['event_hashes'] = {v: digest(out / (v + '-events.dat')) for v in ['reference', 'corrected']}
    save()
    print(json.dumps(dict(checks=checks, window_verification=report['window_verification'],
                         comparisons=report['comparisons']), indent=2))
    if not all(checks.values()):
        raise RuntimeError('GPU comparison failed; see saved diagnostics')


if __name__ == '__main__':
    main()
