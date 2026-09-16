"""Exercise inclusive date ranges and sparse archives with the real executable."""
import datetime as dt
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

EXECUTABLE = Path(sys.argv.pop(1)).resolve()
BUILD = Path(sys.argv.pop(1)).resolve()


class DateRangeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='date-range-', dir=BUILD)
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        (self.work / 'input').mkdir()
        (self.work / 'catalog').write_text(
            '104 001 00 00 0 140 35 10 1e27\n125 001 00 00 0 140 35 10 1e27\n')
        self.env = {k: v for k, v in os.environ.items() if not k.startswith('AUTOFOCUSING_')}
        self.env.update(AUTOFOCUSING_BACKEND='cpu', AUTOFOCUSING_GPU_POWER='off', OMP_NUM_THREADS='1')

    def run_dates(self, start, end, year=None):
        args = [str(EXECUTABLE), year or start[:4], 'test', 'run', 'input', 'catalog',
                'results', 'horizontal', start, end]
        return subprocess.run(args, cwd=self.work, env=self.env, capture_output=True, text=True)

    def day(self, date):
        value = dt.date.fromisoformat(date)
        path = self.work / 'input' / str(value.year) / value.strftime('%m%d')
        path.mkdir(parents=True)
        return path

    def test_leap_day_and_empty_interval(self):
        result = self.run_dates('2004-02-28', '2004-03-01')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('#ScanSummary days=3 missing=3 empty=0 multiple=0 loaded=0', result.stderr)
        self.assertIn('No input files loaded', result.stderr)
        self.assertEqual(len(list((self.work / 'results/test/run').glob('2004_*.dat'))), 1)

    def test_year_boundary_and_single_day(self):
        for start, end, days in [('2024-12-31', '2025-01-01', 2),
                                 ('2025-01-01', '2025-01-01', 1)]:
            result = self.run_dates(start, end)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(f'#ScanSummary days={days} missing={days}', result.stderr)

    def test_explicit_range_is_not_capped(self):
        start, end = '2004-01-01', '2025-01-01'
        result = self.run_dates(start, end)
        days = (dt.date.fromisoformat(end) - dt.date.fromisoformat(start)).days + 1
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f'#ScanSummary days={days} missing={days}', result.stderr)

    def test_missing_empty_and_multiple_days_are_skipped(self):
        self.day('2004-01-02')
        multiple = self.day('2004-01-03')
        (multiple / 'one.h5').write_bytes(b'invalid')
        (multiple / 'two.h5').write_bytes(b'invalid')
        # A malformed singleton demonstrates that the scanner reaches the later
        # existing file after the gaps. Corrupt input must fail, not be skipped.
        selected = self.day('2004-01-04') / 'selected.h5'
        selected.write_bytes(b'invalid')
        result = self.run_dates('2004-01-01', '2004-01-04')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('#Startinput/2004/0104/selected.h5', result.stderr)
        self.assertIn('Cannot open HDF5: input/2004/0104/selected.h5', result.stderr)
        self.assertNotIn('#Startinput/2004/0103/', result.stderr)
        result = self.run_dates('2004-01-01', '2004-01-03')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('#ScanSummary days=3 missing=1 empty=1 multiple=1 loaded=0', result.stderr)

    def test_outside_files_are_not_opened_and_start_is_inclusive(self):
        for day in ('2004-02-27', '2004-03-02'):
            (self.day(day) / 'outside.h5').write_bytes(b'invalid')
        result = self.run_dates('2004-02-28', '2004-03-01')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('#Startinput', result.stderr)
        selected = self.day('2004-02-28') / 'first.h5'
        selected.write_bytes(b'invalid')
        result = self.run_dates('2004-02-28', '2004-03-01')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('#Startinput/2004/0228/first.h5', result.stderr)

    def test_invalid_dates_fail_before_creating_output(self):
        for start, end, year in [('2004-02-30', '2004-03-01', '2004'),
                                 ('2005-02-29', '2005-03-01', '2005'),
                                 ('2004-2-01', '2004-03-01', '2004'),
                                 ('2004-01-02', '2004-01-01', '2004'),
                                 ('2004-01-01', '2004-01-02', '2005'),
                                 ('1399-01-01', '2004-01-01', '2004')]:
            result = self.run_dates(start, end, year)
            self.assertNotEqual(result.returncode, 0, (start, end))
            self.assertIn('DATE', result.stderr)
            self.assertFalse((self.work / 'results').exists())


if __name__ == '__main__':
    unittest.main()
