"""Summarize paired matrix outputs after detected-component max/MAD selection."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics


def distribution(values, unit):
    if not values:
        return dict(count=0,unit=unit)
    ab=sorted(abs(x) for x in values)
    p=.95*(len(ab)-1);i=int(p)
    return dict(count=len(values),unit=unit,mean_signed=statistics.mean(values),
                median_signed=statistics.median(values),median_absolute=statistics.median(ab),
                p95_absolute=ab[i]+(p-i)*(ab[min(i+1,len(ab)-1)]-ab[i]),
                maximum_absolute=ab[-1])


def summarize(rows):
    result=dict(selected_events=len(rows),ratios={},fractions={},diagonal_nonpositive={})
    for name,c in [('R',0),('T',1),('U',2)]:
        result['diagonal_nonpositive'][name]={version:sum(r[version][c]<=0 for r in rows)
                                             for version in ['old','corrected']}
    for name,a,b in [('R/U',0,2),('T/U',1,2),('T/R',1,0)]:
        eligible=[r for r in rows if all(r[v][c]>0 for v in ['old','corrected'] for c in [a,b])]
        old=[r['old'][a]/r['old'][b] for r in eligible]
        new=[r['corrected'][a]/r['corrected'][b] for r in eligible]
        entry=distribution([100*(x/y-1) for x,y in zip(old,new)],'percent')
        entry.update(excluded_nonpositive=len(rows)-len(eligible),
                     median_old_ratio=statistics.median(old) if old else None,
                     median_corrected_ratio=statistics.median(new) if new else None,
                     event_ids=[r['id'] for r in eligible])
        result['ratios'][name]=entry
    common=[r for r in rows if all(x>0 for v in ['old','corrected'] for x in r[v])]
    result['all_components_positive_events']=len(common)
    for name,c in [('R',0),('T',1),('U',2)]:
        old=[r['old'][c]/sum(r['old']) for r in common]
        new=[r['corrected'][c]/sum(r['corrected']) for r in common]
        result['fractions'][name]=dict(change=distribution([100*(x-y) for x,y in zip(old,new)],'percentage points'),
                                      mean_old_fraction=statistics.mean(old) if old else None,
                                      mean_corrected_fraction=statistics.mean(new) if new else None)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',nargs='+',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    a=parser.parse_args()
    rows=[]; sources=[]
    for run in a.runs:
        buffers=[(run/(v+'-events.dat')).read_bytes() for v in ['reference','corrected']]
        old,new=[[line.split() for line in data.decode().splitlines()] for data in buffers]
        if len(old)!=len(new) or any(len(r)!=38 for r in old+new):
            raise RuntimeError('Event count or column mismatch')
        sources.append(dict(run=str(run),event_hashes={v:hashlib.sha256(data).hexdigest()
                                                      for v,data in zip(['reference','corrected'],buffers)}))
        for n,(x,y) in enumerate(zip(old,new),1):
            if x[:33]!=y[:33]:
                raise RuntimeError('Non-matrix fields changed; selection may be inconsistent')
            old_diag=list(map(float,x[33:36]));new_diag=list(map(float,y[33:36]))
            if not all(math.isfinite(v) for v in old_diag+new_diag):
                raise RuntimeError('Nonfinite diagonal entry')
            mad=float(x[12]);peak=float(x[11]);component=int(x[1]);threshold=35 if component==2 else 7
            score=peak/mad if mad>0 else None
            rows.append(dict(id=f'{len(sources)}:{n}',start=x[0],end=x[18],component=component,
                             max_over_mad=score,threshold=threshold,selected=score is not None and score>threshold,
                             old=old_diag,corrected=new_diag))
    selected=[r for r in rows if r['selected']]
    result=dict(sources=sources,thresholds={'R':7,'T':7,'U_or_V':35},comparison='strictly greater than threshold; selected using output columns 12/13, identical in both versions',
                interpretation='Ratios use same fitted event. Positive-pair selection is not a significance test of each component; no clipping or absolute-value replacement of powers.',
                total_events=len(rows),time_range=[min(r['start'] for r in rows),max(r['end'] for r in rows)] if rows else None,
                by_component={name:summarize([r for r in selected if r['component']==c]) for name,c in [('R',0),('T',1),('U',2)]},
                pooled=summarize(selected),
                maximum_score_by_component={name:max((r['max_over_mad'] for r in rows if r['component']==c and r['max_over_mad'] is not None),default=None) for name,c in [('R',0),('T',1),('U',2)]},
                events=rows)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(total_events=len(rows),selected={k:v['selected_events'] for k,v in result['by_component'].items()},
                         pooled_ratios=result['pooled']['ratios']),indent=2))


if __name__=='__main__':
    main()
