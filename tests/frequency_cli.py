"""Frequency setting, validation and run-record checks against the real binary."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

REPO, BIN, BUILD = (Path(p).resolve() for p in sys.argv[1:4])
del sys.argv[1:4]


class FrequencyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='frequency-', dir=BUILD)
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.env = {k: v for k, v in os.environ.items() if not k.startswith('AUTOFOCUSING_')}
        self.env.update(AUTOFOCUSING_BIN=str(BIN), AUTOFOCUSING_BACKEND='cpu', OMP_NUM_THREADS='1')
        (self.work / 'input').mkdir()
        (self.work / 'catalog').write_text('104 001 00 00 0 140 35 10 1e27\n')
        (self.work / 'local_config.sh').write_text('HINET_ROOT=input\nCMT_CATALOG=catalog\n')
        self.config = self.work / 'analysis/trial/config.sh'
        self.config.parent.mkdir(parents=True)
        self.base = 'START_DATE=2004-01-01\nEND_DATE=2004-01-01\n'
        self.config.write_text(self.base)

    def command(self, args):
        return subprocess.run(args, cwd=self.work, env=self.env, capture_output=True, text=True)

    def run_analysis(self):
        return self.command(['bash', str(REPO / 'Scripts/run.sh'),
                             '--workspace', str(self.work), '--experiment', 'trial'])

    def test_default_and_primary_metadata(self):
        default = self.command([str(BIN), '--frequency-info'])
        self.assertEqual(default.returncode, 0, default.stderr)
        info = json.loads(default.stdout)
        self.assertEqual((info['min_bin'], info['max_bin'], info['spectrum_bins']), (102, 256, 266))
        self.config.write_text(self.base + 'AUTOFOCUSING_FREQ_MIN=0.05\nAUTOFOCUSING_FREQ_MAX=0.1\n')
        result = self.run_analysis()
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest_file, = (self.work / 'results/trial').glob('*/manifest.json')
        manifest = json.loads(manifest_file.read_text())
        band = manifest['frequency_band']
        self.assertEqual(band['requested_min_hz'], .05)
        self.assertEqual(band['requested_max_hz'], .1)
        self.assertEqual((band['min_bin'], band['max_bin']), (51, 102))
        self.assertEqual((band['min_hz'], band['max_hz']), (51/1024, 102/1024))
        self.assertTrue((manifest_file.parent / '2004_2048_0.049805-0.099609.dat').exists())
        line, = (line for line in (manifest_file.parent / 'run.log').read_text().splitlines()
                 if line.startswith('#FrequencyBand '))
        self.assertEqual(json.loads(line[len('#FrequencyBand '):]), band)

    def test_invalid_settings_do_not_start(self):
        for lo, hi in [('0', '.1'), ('-.05', '.1'), ('nan', '.1'),
                       ('.05', 'inf'), ('', '.1'), ('junk', '.1'),
                       ('.05', '.1junk'), ('.1', '.05'), ('.1', '.1'),
                       ('.05', '1'), ('.05', '1.1'), ('.0001', '.1')]:
            self.config.write_text(self.base +
                                  f"AUTOFOCUSING_FREQ_MIN='{lo}'\nAUTOFOCUSING_FREQ_MAX='{hi}'\n")
            result = self.run_analysis()
            self.assertNotEqual(result.returncode, 0, (lo, hi))
            self.assertIn('AUTOFOCUSING_FREQ_', result.stderr)
            self.assertFalse((self.work / 'results').exists())

    def test_near_nyquist_and_single_bin(self):
        for lo, hi, first, last in [('.8', '.9999', 819, 1023),
                                    ('.05', '.0501', 51, 51)]:
            self.env.update(AUTOFOCUSING_FREQ_MIN=lo, AUTOFOCUSING_FREQ_MAX=hi)
            result = self.command([str(BIN), '--frequency-info'])
            self.assertEqual(result.returncode, 0, result.stderr)
            band = json.loads(result.stdout)
            self.assertEqual((band['min_bin'], band['max_bin']), (first, last))
            self.assertLessEqual(band['spectrum_bins'], 1024)


if __name__ == '__main__':
    unittest.main()
