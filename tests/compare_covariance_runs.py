"""Compare saved float/double covariance runs without relaxing audit checks."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from compare_peak_search_runs import load, selected


def distribution(values):
    values = sorted(values)
    if not values:
        return None
    def quantile(q):
        index = (len(values)-1)*q
        low, high = math.floor(index), math.ceil(index)
        return values[low]*(high-index)+values[high]*(index-low) if low != high else values[low]
    return {'count':len(values), 'median':quantile(.5), 'p95':quantile(.95), 'max':values[-1]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'fixed', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    ra, old, old_log = load(args.baseline)
    rb, new, new_log = load(args.fixed)
    checks = {
        'corrected_audit_passed':all(rb['checks'].values()),
        'same_settings':ra['settings'] == rb['settings'],
        'same_input_hashes':[s['sha256'] for s in ra['sources']] == [s['sha256'] for s in rb['sources']],
        'same_catalog_hash':ra['catalog_sha256'] == rb['catalog_sha256'],
        'same_segments':ra['segments'] == rb['segments'],
        'same_candidates':old.keys() == new.keys(),
        'same_optimizer_records':old_log == new_log,
        'same_events':True, 'non_covariance_columns_exact':True,
        'all_output_finite':True, 'same_selection':True,
    }
    changes, selected_changes, recovered = [], [], []
    negative_diagonal_rows = {'baseline':[], 'fixed':[]}
    changed_rows = sign_changes = events = 0
    for key in sorted(old.keys() & new.keys()):
        a, b = old[key]['row'], new[key]['row']
        checks['same_events'] &= (a is None) == (b is None)
        if a is None or b is None:
            continue
        events += 1
        checks['non_covariance_columns_exact'] &= a[:20]+a[30:] == b[:20]+b[30:]
        checks['same_selection'] &= selected(a) == selected(b)
        checks['all_output_finite'] &= all(math.isfinite(float(v)) for i,v in enumerate(b) if i not in (0,18))
        changed_rows += a[20:30] != b[20:30]
        for label, row in [('baseline',a), ('fixed',b)]:
            if any(float(row[i]) < 0 for i in (20,24,27,29)):
                negative_diagonal_rows[label].append({'key':key, 'selected':selected(row)})
        if any(not math.isfinite(float(v)) for v in a[20:30]):
            recovered.append({'key':key, 'selected':selected(b), 'row':b})
            continue
        for av, bv in zip(a[20:30], b[20:30]):
            x, y = float(av), float(bv)
            error = abs(y-x)/max(abs(x),1e-300)
            changes.append(error)
            if selected(b):
                selected_changes.append(error)
            sign_changes += x*y < 0
    result = {'checks':checks, 'baseline_checks':ra['checks'], 'events':events,
              'changed_covariance_rows':changed_rows, 'recovered_nonfinite_rows':recovered,
              'relative_covariance_change_all_finite_entries':distribution(changes),
              'relative_covariance_change_selected_entries':distribution(selected_changes),
              'finite_covariance_sign_changes':sign_changes, 'selected_counts':rb['threshold_counts']}
    result['negative_diagonal_rows'] = negative_diagonal_rows
    files = [args.baseline/'report.json', args.fixed/'report.json',
             args.baseline/'corrected-events.dat', args.fixed/'corrected-events.dat', Path(__file__)]
    result['sha256'] = {str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'sha256'},indent=2))
    if not all(checks.values()):
        raise RuntimeError('Covariance comparison failed')


if __name__ == '__main__':
    main()
