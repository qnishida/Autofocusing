"""Run a bounded 3c matrix audit, recording both matrices and input hashes."""
import argparse
import datetime
import json
import math
import os
from pathlib import Path
import re
import subprocess

from run_spectral_matrix_events import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('binary', 'hinet', 'catalog', 'output'):
        p.add_argument(name, type=Path)
    p.add_argument('--start', required=True)
    p.add_argument('--days', type=int, required=True)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--expected', type=Path,
                   help='saved independent reference-events.dat/corrected-events.dat pair')
    a = p.parse_args()
    if min(a.days,a.threads)<1:
        p.error('days and threads must be positive')
    repo = Path(__file__).resolve().parents[1]
    start = datetime.date.fromisoformat(a.start)
    end = start+datetime.timedelta(days=a.days-1)
    out = a.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    env = {k:v for k,v in os.environ.items() if not k.startswith(('AUTOFOCUSING_','OMP_','KMP_','MATRIX_'))}
    settings = dict(AUTOFOCUSING_BACKEND='cpu',AUTOFOCUSING_GPU_POWER='off',
                    AUTOFOCUSING_FREQ_MIN='0.1',AUTOFOCUSING_FREQ_MAX='0.25',
                    AUTOFOCUSING_SLOWNESS_MAX='0.165',AUTOFOCUSING_SLOWNESS_STEP='0.005',
                    AUTOFOCUSING_PROFILE='1',OMP_NUM_THREADS=str(a.threads),OMP_DYNAMIC='FALSE')
    env.update(settings)
    env['MATRIX_REFERENCE_OUTPUT']=str(out/'reference-events.dat')
    report=dict(start=str(start),end=str(end),component_mode='3c',settings=settings,
                sources=[],checks={},binary_sha256=digest(a.binary.resolve()),
                catalog_sha256=digest(a.catalog.resolve()),
                seed='1837 + 104729 * replicate',fftw='ESTIMATE | UNALIGNED')
    def save():
        (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    for i in range(a.days):
        day=start+datetime.timedelta(days=i)
        source=a.hinet.resolve()/day.strftime('%Y/%m%d/%Y%j0000.h5')
        before=source.stat();sha=digest(source);after=source.stat()
        if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
            raise RuntimeError(f'Input is changing: {source}')
        report['sources'].append(dict(path=str(source),size=after.st_size,mtime_ns=after.st_mtime_ns,sha256=sha))
        link=out/'input'/day.strftime('%Y/%m%d')/source.name
        link.parent.mkdir(parents=True)
        link.symlink_to(source)
    report['frequency_band']=json.loads(subprocess.check_output([str(a.binary.resolve()),'--frequency-info'],env=env,text=True))
    assert (report['frequency_band']['min_bin'],report['frequency_band']['max_bin'])==(102,256)
    command=[str(a.binary.resolve()),str(start.year),'matrix-audit','corrected',str(out/'input'),
             str(a.catalog.resolve()),str(out/'results'),'3c',str(start),str(end)]
    report['command']=command;save()
    print(f'Running {start} through {end}',flush=True)
    with (out/'run.log').open('w') as log:
        process=subprocess.run(command,cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT,
                               timeout=max(3600,300*a.days))
    report['exit_code']=process.returncode;save();process.check_returncode()
    corrected,=(out/'results/matrix-audit/corrected').glob('*.dat')
    (out/'corrected-events.dat').write_bytes(corrected.read_bytes())
    old=[s.split() for s in (out/'reference-events.dat').read_text().splitlines()]
    new=[s.split() for s in (out/'corrected-events.dat').read_text().splitlines()]
    checks=report['checks']
    checks['columns_valid']=all(len(r)==38 for r in old+new)
    if not checks['columns_valid']:
        save();raise RuntimeError('Unexpected event columns')
    checks['non_matrix_columns_exact']=[r[:33] for r in old]==[r[:33] for r in new]
    checks['all_output_finite']=all(math.isfinite(float(v)) for r in old+new for i,v in enumerate(r) if i not in (0,18))
    log=(out/'run.log').read_text()
    segments=[s for s in log.splitlines() if s.startswith('#Segment ')]
    report['segments']=segments
    report['loaded_station_counts']=[int(s) for s in re.findall(r'^#Read data: load_h5 \S+ (\d+) ',log,re.M)]
    checks['all_days_loaded']=len(report['loaded_station_counts'])==a.days and min(report['loaded_station_counts'])>0
    checks['grid_settings_verified']='#SlantStack step_s_per_km=0.005 max_s_per_km=0.165 grid=67x67' in log
    report['events']=len(new)
    report['threshold_counts']={str(c):sum(int(r[1])==c and float(r[12])>0 and float(r[11])/float(r[12])>t for r in new) for c,t in [(0,7),(1,7),(2,35)]}
    max_error=0.
    for r in new:
        scalar=float(r[32])/(155/1024);diagonal=float(r[33+int(r[1])])
        max_error=max(max_error,abs(scalar-diagonal)/max(abs(scalar),abs(diagonal),1e-30))
    report['diagonal_scalar_max_scaled_error']=max_error
    checks['diagonal_matches_scalar']=max_error<2e-5
    if a.expected:
        for label in ['reference','corrected']:
            checks[label+'_matches_independent_run']=(out/(label+'-events.dat')).read_bytes()==(a.expected/(label+'-events.dat')).read_bytes()
    checks['inputs_unchanged']=all(digest(Path(s['path']))==s['sha256'] for s in report['sources'])
    checks['catalog_unchanged']=digest(a.catalog.resolve())==report['catalog_sha256']
    report['event_hashes']={label:digest(out/(label+'-events.dat')) for label in ['reference','corrected']}
    save()
    print(json.dumps(dict(events=report['events'],threshold_counts=report['threshold_counts'],checks=checks),indent=2),flush=True)
    if not all(checks.values()):
        raise RuntimeError('Audit check failed')


if __name__=='__main__':
    main()
