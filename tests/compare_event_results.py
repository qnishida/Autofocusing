"""Compare event output while excluding the existing time-seeded bootstrap fields."""
import math
from pathlib import Path
import sys

gpu = len(sys.argv) == 4 and sys.argv[3] in ('--metal', '--gpu', '--cuda')
if len(sys.argv) != 3 and not gpu:
    raise SystemExit('Usage: compare_event_results.py baseline.dat optimized.dat [--gpu|--metal|--cuda]')
rows = [[line.split() for line in Path(path).read_text().splitlines()]
        for path in sys.argv[1:3]]
if not rows[0] or len(rows[0]) != len(rows[1]):
    raise SystemExit('Empty result or event count mismatch')
for index, (before, after) in enumerate(zip(*rows)):
    if len(before) != 38 or len(after) != 38:
        raise SystemExit('Expected 38 output columns')
    for column, (a,b) in enumerate(zip(before,after), start=1):
        # Columns 21–32: covariance, bootstrap beam power and its spread.
        if 21 <= column <= 32:
            continue
        if column == 15:  # Near-zero relative optimizer convergence residual.
            ok = math.isclose(float(a), float(b), rel_tol=0, abs_tol=1e-12)
        elif gpu and column in (12, 13):  # Float beam maximum and MAD.
            ok = math.isclose(float(a), float(b), rel_tol=1e-4, abs_tol=1e-30)
        else:
            # Source locations, powers, times and flags must match at output precision.
            ok = a == b
        if not ok:
            raise SystemExit(f'Event {index}, column {column}: {a} != {b}')
print(f'PASS: {len(rows[0])} event identities and non-bootstrap outputs match '
      '(optimizer residual tolerance 1e-12' +
      ('; GPU beam maximum/MAD rtol=1e-4)' if gpu else ')'))
