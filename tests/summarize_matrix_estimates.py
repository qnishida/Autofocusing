"""Describe catalog estimate distributions and paired signed implementation shifts.

These are across-event distributions, not repeated-estimation uncertainty or
bias against a known physical truth. Input is analyze_matrix_ratios.py JSON.
"""
import argparse
import hashlib
import json
from pathlib import Path
import statistics


def describe(values):
    if not values:
        return {'count': 0}
    ordered = sorted(values)

    def quantile(p):
        position = p * (len(ordered) - 1)
        i = int(position)
        return ordered[i] + (position - i) * (ordered[min(i + 1, len(ordered) - 1)] - ordered[i])

    return dict(count=len(values), mean=statistics.mean(values),
                sample_sd=statistics.stdev(values) if len(values) > 1 else None,
                q25=quantile(.25), median=quantile(.5), q75=quantile(.75),
                minimum=min(values), maximum=max(values))


def paired(original, corrected):
    delta = [x - y for x, y in zip(original, corrected)]
    return dict(original=describe(original), corrected=describe(corrected),
                original_minus_corrected=describe(delta),
                original_higher=sum(x > 0 for x in delta),
                original_lower=sum(x < 0 for x in delta),
                equal=sum(x == 0 for x in delta))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    data = json.loads(raw)
    result = dict(input=str(args.input), input_sha256=hashlib.sha256(raw).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  interpretation=__doc__.strip(),
                  spread='Across selected catalog rows; sample SD uses n-1; quartiles use linear interpolation.',
                  shift='Original minus corrected on the same fitted event; corrected is not known physical truth.',
                  by_component={})
    for c, name in enumerate(['R', 'T', 'U']):
        rows = [e for e in data['events'] if e['selected'] and e['component'] == c]
        cohort = dict(selected_events=len(rows), unique_start_dates=len({e['start'][:8] for e in rows}),
                      ratios={}, fractions_percent={})
        for label, a, b in [('R/U', 0, 2), ('T/U', 1, 2), ('T/R', 1, 0)]:
            eligible = [e for e in rows if all(e[v][i] > 0 for v in ['old', 'corrected'] for i in [a, b])]
            stats = paired([e['old'][a] / e['old'][b] for e in eligible],
                           [e['corrected'][a] / e['corrected'][b] for e in eligible])
            stats.update(excluded_nonpositive=len(rows) - len(eligible),
                         event_ids=[e['id'] for e in eligible], unit='dimensionless ratio')
            cohort['ratios'][label] = stats
        common = [e for e in rows if all(x > 0 for v in ['old', 'corrected'] for x in e[v])]
        for a, label in enumerate(['R', 'T', 'U']):
            stats = paired([100 * e['old'][a] / sum(e['old']) for e in common],
                           [100 * e['corrected'][a] / sum(e['corrected']) for e in common])
            stats.update(event_ids=[e['id'] for e in common],
                         estimate_unit='percent of diagonal sum', shift_unit='percentage points')
            cohort['fractions_percent'][label] = stats
        result['by_component'][name] = cohort
    args.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
