"""Compare peak-search candidates by window, component and original slowness seed."""
import argparse,collections,hashlib,json,math
from pathlib import Path


def load(path):
    report=json.loads((path/'report.json').read_text())
    records={k:[] for k in ['PeakSeed','ConvergenceResult','ConvergenceEvent','ConvergenceStep']}
    for line in (path/'run.log').read_text().splitlines():
        for kind in records:
            if line.startswith('#'+kind+' '):records[kind].append(json.loads(line.split(' ',1)[1]))
    seeds={r['id']:r for r in records['PeakSeed']}
    results={r['id']:r for r in records['ConvergenceResult']}
    rows=[s.split() for s in (path/'corrected-events.dat').read_text().splitlines()]
    events={r['id']:row for r,row in zip(records['ConvergenceEvent'],rows)}
    assert len(rows)==len(records['ConvergenceEvent']) and seeds.keys()==results.keys()
    keyed={}
    for ident,s in seeds.items():
        key=(s['start'],s['end'],s['component'],s['px'],s['py'])
        assert key not in keyed
        keyed[key]={'seed':s,'result':results[ident],'row':events.get(ident)}
    return report,keyed,records


def selected(r):return r is not None and float(r[12])>0 and float(r[11])/float(r[12])>(35 if r[1]=='2' else 7)


def summarize(items):
    rows=[v['row'] for v in items if v['row'] is not None]
    return {'candidates':len(items),'outcomes':dict(collections.Counter(v['result']['reason'] for v in items)),
            'events':len(rows),'selected_counts':{c:sum(int(r[1])==i and selected(r) for r in rows) for i,c in enumerate(['R','T','U'])}}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('baseline',type=Path);p.add_argument('fixed',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    ra,old,oldlog=load(a.baseline);rb,new,newlog=load(a.fixed)
    checks={'baseline_run_passed':all(ra['checks'].values()),'fixed_run_passed':all(rb['checks'].values()),
            'settings_equal':ra['settings']==rb['settings'],'segments_and_windows_equal':ra['segments']==rb['segments'],
            'input_hashes_equal':[s['sha256'] for s in ra['sources']]==[s['sha256'] for s in rb['sources']],
            'catalog_hash_equal':ra['catalog_sha256']==rb['catalog_sha256']}
    # Preserve a failed run verdict while allowing diagnostic comparison.
    # Do not relax the runner's finite-output requirement.
    required=['settings_equal','segments_and_windows_equal','input_hashes_equal','catalog_hash_equal']
    assert all(checks[k] for k in required),checks
    assert ra.get('exit_code')==0 and rb.get('exit_code')==0
    common=sorted(old.keys()&new.keys());added=sorted(new.keys()-old.keys());lost=sorted(old.keys()-new.keys())
    result={'checks':checks,'baseline':summarize(list(old.values())),'fixed':summarize(list(new.values())),
            'common_candidates':len(common),'added':summarize([new[k] for k in added]),'lost':summarize([old[k] for k in lost]),
            'added_candidates':[new[k] for k in added],'lost_candidates':[old[k] for k in lost],
            'common_result_changes':[],'common_catalog_changes':[]}
    for k in common:
        x,y=old[k],new[k]
        if {i:v for i,v in x['result'].items() if i!='id'}!={i:v for i,v in y['result'].items() if i!='id'}:
            result['common_result_changes'].append({'key':k,'baseline':x['result'],'fixed':y['result']})
        if x['row']!=y['row']:result['common_catalog_changes'].append({'key':k,'baseline':x['row'],'fixed':y['row']})
    result['nonfinite_events']={name:[{'key':key,'row':value['row'],'selected':selected(value['row'])} for key,value in items.items() if value['row'] is not None and any(not math.isfinite(float(v)) for i,v in enumerate(value['row']) if i not in (0,18))] for name,items in [('baseline',old),('fixed',new)]}
    result['common_event_count']=sum(old[k]['row'] is not None and new[k]['row'] is not None for k in common)
    result['all_existing_candidates_retained']=not lost
    result['all_common_event_rows_exact']=not result['common_catalog_changes']
    result['added_seed_radius_s_per_km']=[math.hypot(k[3],k[4]) for k in added]
    files=[a.baseline/'report.json',a.fixed/'report.json',a.baseline/'run.log',a.fixed/'run.log',Path(__file__)]
    result['sha256']={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['added_candidates','lost_candidates','sha256']},indent=2))


if __name__=='__main__':main()
