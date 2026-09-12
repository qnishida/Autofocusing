#!/usr/bin/env python3
"""Compare the frozen 51f60a7 loader with the common CPU I/O implementation."""
import argparse
import json
import platform
from pathlib import Path
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path, help='built test_io executable')
    parser.add_argument('input', type=Path, help='one real daily HDF5 file')
    parser.add_argument('output', type=Path, help='JSON report destination')
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 1 or not args.input.is_file():
        parser.error('provide an existing HDF5 file and positive repeat count')
    binary, source = str(args.binary.resolve()), str(args.input.resolve())
    # Validate every stored sample before timing; also warms the OS cache.
    subprocess.run([binary, '--real', source], check=True)
    records = []
    summaries = []
    for threads in (1, 4, 8, 16):
        groups = {mode: [] for mode in ('reference', 'current')}
        # Alternate process order to reduce systematic thermal/order bias.
        for trial in range(args.repeats):
            for mode in (('reference', 'current') if trial % 2 == 0 else ('current', 'reference')):
                run = subprocess.run([binary, '--benchmark', source, mode, str(threads), '1'],
                                     check=True, text=True, capture_output=True)
                rows = [json.loads(line) for line in run.stdout.splitlines() if line.startswith('{')]
                if len(rows) != 1:
                    raise RuntimeError(f'Expected one timing record: {run.stdout}')
                row = rows[0]
                row['trial'] = trial
                records.append(row)
                groups[mode].append(row['wall_s'])
        summary = {'threads': threads}
        for mode, values in groups.items():
            summary[mode] = {'median_s': statistics.median(values), 'min_s': min(values), 'max_s': max(values)}
        summary['speedup'] = summary['reference']['median_s'] / summary['current']['median_s']
        summaries.append(summary)
        print(json.dumps(summary), flush=True)
    report = {'platform': platform.platform(), 'input': source, 'baseline': '51f60a7',
              'cache': 'warmed by full equivalence check; not a cold USB throughput measurement',
              'rss_unit': 'bytes on macOS, KiB on Linux; per-process high water',
              'records': records, 'summary': summaries}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
