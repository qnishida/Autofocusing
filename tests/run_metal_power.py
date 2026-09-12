"""Paired real-event qualification and timings for opt-in Metal objectives."""
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

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('binary',type=Path)
p.add_argument('hinet',type=Path)
p.add_argument('catalog',type=Path)
p.add_argument('output',type=Path)
p.add_argument('--days',type=int,default=3,choices=(3,5))
p.add_argument('--repeats',type=int,default=5)
p.add_argument('--legacy-run',type=Path,help='reuse a completed legacy-0 run from this script')
p.add_argument('--legacy',type=Path,help='optional frozen f745119 driver')
a=p.parse_args()
out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
for day in range(1,a.days+1):
    source=a.hinet.resolve()/'2004'/f'010{day}'/f'200400{day}0000.h5'
    if not source.is_file(): raise RuntimeError(f'Missing {source}')
    target=out/'input/2004'/f'010{day}';target.mkdir(parents=True)
    (target/source.name).symlink_to(source)
report={'platform':platform.platform(),'baseline_commit':'f745119','threads':16,'days':a.days,
        'fftw':'ESTIMATE | UNALIGNED','seed':'1837 + 104729 * replicate',
        'cache':'OS cache not purged; rotate mode order; includes initialization/packing/sync',
        'runs':[]}
reference=None; reference_choices=None
if a.legacy_run:
    files=list((a.legacy_run/'results/power/legacy-0').glob('*.dat'))
    if len(files)!=1: raise RuntimeError('Expected legacy result')
    reference=[line.split() for line in files[0].read_text().splitlines()]
    reference_choices=[line for line in (a.legacy_run/'legacy-0.log').read_text().splitlines() if line.startswith(('#deg max=','#dp_Δ max='))]
    report['legacy_event_sha256']=hashlib.sha256(files[0].read_bytes()).hexdigest()
modes=['off','bootstrap','grid','all']
jobs=([('legacy',0)] if a.legacy else [])+[(mode,r) for r in range(a.repeats) for mode in (modes if r%2==0 else list(reversed(modes)))]
for mode,repeat in jobs:
    tag=f'{mode}-{repeat}';exe=(a.legacy if mode=='legacy' else a.binary).resolve()
    env=dict(os.environ,OMP_NUM_THREADS='16',AUTOFOCUSING_BACKEND='metal',AUTOFOCUSING_METAL_POWER='off' if mode=='legacy' else mode,AUTOFOCUSING_PROFILE='1')
    for key in ('AUTOFOCUSING_VERIFY_METAL','AUTOFOCUSING_SLOWNESS_STEP','AUTOFOCUSING_SLOWNESS_MAX'): env.pop(key,None)
    print(f'Running {tag}',flush=True)
    start=time.monotonic();log=out/f'{tag}.log'
    with log.open('w') as stream:
        timing=['/usr/bin/time','-l'] if platform.system()=='Darwin' else []
        subprocess.run(timing+[str(exe),'2004','power',tag,str(out/'input'),str(a.catalog.resolve()),str(out/'results'),'horizontal'],env=env,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=900)
    wall=time.monotonic()-start
    files=list((out/'results/power'/tag).glob('*.dat'))
    if len(files)!=1: raise RuntimeError('Expected one event output')
    text=files[0].read_text();rows=[line.split() for line in text.splitlines()]
    if len(rows)!=(20 if a.days==3 else 29): raise RuntimeError(f'Unexpected event count: {len(rows)}')
    choices=[line for line in log.read_text().splitlines() if line.startswith(('#deg max=','#dp_Δ max='))]
    if reference is None: reference,reference_choices=rows,choices
    if choices!=reference_choices or len(rows)!=len(reference): raise RuntimeError('Initial candidates/event counts differ')
    max_rel=0.
    for row,ref in zip(rows,reference):
        if len(row)!=38: raise RuntimeError('Expected 38 columns')
        for col,(value,expected) in enumerate(zip(row,ref),1):
            if 21<=col<=33 and mode in ('bootstrap','all'):
                x,y=float(value),float(expected)
                if not math.isfinite(x) or not math.isclose(x,y,rel_tol=1e-3,abs_tol=1e-30):
                    raise RuntimeError(f'{tag} column {col}: {x} != {y}')
                max_rel=max(max_rel,abs(x-y)/max(abs(y),1e-30))
            elif col==15:
                if not math.isclose(float(value),float(expected),rel_tol=0,abs_tol=1e-12): raise RuntimeError('Convergence residual differs')
            elif value!=expected: raise RuntimeError(f'{tag} column {col}: {value} != {expected}')
    stage={}
    rss=None
    for line in log.read_text().splitlines():
        if 'maximum resident set size' in line: rss=int(line.split()[0])
        if line.startswith('#POWER_PROFILE '):
            fields=dict(x.split('=',1) for x in line.split()[1:])
            stage[fields['stage']]=stage.get(fields['stage'],0)+float(fields['total_s'])
    report['runs'].append({'mode':mode,'repeat':repeat,'wall_s':wall,'stages':stage,'events':len(rows),'max_rss_bytes_macos':rss,
                          'choices':choices,'max_bootstrap_relative':max_rel,
                          'event_sha256':hashlib.sha256(text.encode()).hexdigest()})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {tag}: {len(rows)} events, max_boot_rel={max_rel:.3g}, wall={wall:.3f}s stages={stage}',flush=True)
summary={}
for mode in modes:
    r=[x for x in report['runs'] if x['mode']==mode]
    summary[mode]={'median_s':statistics.median(x['wall_s'] for x in r),'min_s':min(x['wall_s'] for x in r),'max_s':max(x['wall_s'] for x in r),
                   'stages':{stage:statistics.median(x['stages'][stage] for x in r) for stage in ('grid','bootstrap')}}
for mode in modes[1:]:
    stage='bootstrap' if mode=='bootstrap' else 'grid'
    summary[mode]['stage_speedup']=summary['off']['stages'][stage]/summary[mode]['stages'][stage]
    summary[mode]['end_to_end_pass']=summary[mode]['median_s']<=1.05*summary['off']['median_s']
report['summary']=summary
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(summary),flush=True)
