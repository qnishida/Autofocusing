"""Match a historical 38-column catalog to archived paired matrix-audit runs."""
import argparse
from collections import Counter, defaultdict
import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path

from summarize_matrix_estimates import describe


def selected(row):
    return float(row[12]) > 0 and float(row[11]) / float(row[12]) > (35 if int(row[1]) == 2 else 7)


def distance(a, b):
    return math.hypot(float(a[5]) - float(b[5]), float(a[6]) - float(b[6]))


def metadata(item):
    ident, row = item
    return dict(id=ident, start=row[0], end=row[18], component=int(row[1]),
                px=float(row[5]), py=float(row[6]), selected=selected(row),
                max_over_mad=float(row[11]) / float(row[12]) if float(row[12]) > 0 else None)


def summarize(pairs):
    def relative(value):
        values = [(value(a), value(b)) for _, a, _, b in pairs]
        values = [(a, b) for a, b in values if a > 0 and b > 0]
        differences = [100 * (a / b - 1) for a, b in values]
        return dict(signed_percent=describe(differences), absolute_percent=describe([abs(x) for x in differences]))

    return dict(count=len(pairs),
                slowness_distance_s_per_km=describe([distance(a, b) for _, a, _, b in pairs]),
                source_distance_difference_deg=describe([abs(float(a[2]) - float(b[2])) for _, a, _, b in pairs]),
                accepted_window_count_equal=sum(a[13] == b[13] for _, a, _, b in pairs),
                peak_frequency_equal=sum(a[15] == b[15] for _, a, _, b in pairs),
                array_center_equal=sum(a[16:18] == b[16:18] for _, a, _, b in pairs),
                max_over_mad=relative(lambda a: float(a[11]) / float(a[12]) if float(a[12]) > 0 else 0),
                max_power=relative(lambda a: float(a[11])),
                mad=relative(lambda a: float(a[12])),
                MS_max=relative(lambda a: float(a[32])),
                fitted_component_original_matrix=relative(lambda a: float(a[33 + int(a[1])])),
                fit_shift_dp_absolute_difference=describe([abs(float(a[8]) - float(b[8])) for _, a, _, b in pairs]),
                fit_shift_dD_absolute_difference=describe([abs(float(a[9]) - float(b[9])) for _, a, _, b in pairs]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('historical', type=Path)
    p.add_argument('expanded', type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--max-slowness-distance', type=float, default=.005)
    args = p.parse_args()
    if args.max_slowness_distance <= 0:
        p.error('slowness distance must be positive')
    data = json.loads(args.expanded.read_text())
    days = set()
    current = []
    for i, source in enumerate(data['sources'], 1):
        run = Path(source['run'])
        report = json.loads((run / 'report.json').read_text())
        start = datetime.date.fromisoformat(report['start'])
        end = datetime.date.fromisoformat(report['end'])
        dates = {(start + datetime.timedelta(days=j)).strftime('%Y%m%d') for j in range((end - start).days + 1)}
        if days & dates:
            raise RuntimeError('Repeated input dates')
        days.update(dates)
        raw = (run / 'reference-events.dat').read_bytes()
        if hashlib.sha256(raw).hexdigest() != source['event_hashes']['reference']:
            raise RuntimeError('Archived event hash mismatch')
        current.extend((f'{i}:{n}', line.split()) for n, line in enumerate(raw.decode().splitlines(), 1))
    raw = args.historical.read_bytes()
    historical = [(n, line.split()) for n, line in enumerate(raw.decode().splitlines(), 1) if line[:8] in days]
    if not all(len(row) == 38 for _, row in historical + current):
        raise RuntimeError('Expected 38 columns')
    groups = []
    for items in [historical, current]:
        group = defaultdict(list)
        for item in items:
            row = item[1]
            group[(row[0], row[18], row[1])].append(item)
        groups.append(group)
    old_groups, new_groups = groups
    pairs, unmatched_old, unmatched_new = [], [], []
    for key in sorted(old_groups.keys() | new_groups.keys()):
        a, b = old_groups[key], new_groups[key]
        if max(len(a), len(b)) > 8:
            raise RuntimeError('Group too large for bounded exhaustive matching')
        if len(a) <= len(b):
            candidates = (list(enumerate(order)) for order in itertools.permutations(range(len(b)), len(a)))
        else:
            candidates = ([(j, i) for i, j in enumerate(order)] for order in itertools.permutations(range(len(a)), len(b)))

        def cost(candidate):
            ds = [distance(a[i][1], b[j][1]) for i, j in candidate]
            return sum(d > args.max_slowness_distance for d in ds), sum(ds)

        assignment = min(candidates, key=cost)
        accepted = [(i, j) for i, j in assignment if distance(a[i][1], b[j][1]) <= args.max_slowness_distance]
        pairs.extend((a[i][0], a[i][1], b[j][0], b[j][1]) for i, j in accepted)
        unmatched_old.extend(metadata(item) for i, item in enumerate(a) if i not in {x[0] for x in accepted})
        unmatched_new.extend(metadata(item) for j, item in enumerate(b) if j not in {x[1] for x in accepted})
    result = dict(historical=str(args.historical), historical_sha256=hashlib.sha256(raw).hexdigest(),
                  expanded=str(args.expanded), expanded_sha256=hashlib.sha256(args.expanded.read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  matching='Exact start/end/component; one-to-one assignment maximizing matches within slowness tolerance, then minimizing total slowness distance. No power values used in matching.',
                  max_slowness_distance_s_per_km=args.max_slowness_distance,
                  comparison='Historical versus current ORIGINAL matrix formula, before the candidate self-term correction. Positive-valued comparisons use 100*(historical/current-1).',
                  input_identity='Current input files are stable by recorded hashes. Historical waveform hashes/source settings are unavailable; file identity is not established.',
                  days=sorted(days), counts={}, all_matched=summarize(pairs), current_selected={},
                  unmatched_historical=unmatched_old, unmatched_current=unmatched_new,
                  pairs=[dict(historical=metadata((i, a)), current=metadata((j, b)), slowness_distance_s_per_km=distance(a, b)) for i, a, j, b in pairs])
    for label, items in [('historical', historical), ('current', current)]:
        result['counts'][label] = dict(total=len(items), all_by_component=dict(Counter(row[1] for _, row in items)),
                                     selected_by_component=dict(Counter(row[1] for _, row in items if selected(row))))
    for c, label in enumerate(['R', 'T', 'U']):
        result['current_selected'][label] = summarize([pair for pair in pairs if selected(pair[3]) and int(pair[3][1]) == c])
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(counts=result['counts'], matched=len(pairs), unmatched_historical=len(unmatched_old),
                          unmatched_current=len(unmatched_new), matched_current_selected={k:v['count'] for k,v in result['current_selected'].items()}), indent=2))


if __name__ == '__main__':
    main()
