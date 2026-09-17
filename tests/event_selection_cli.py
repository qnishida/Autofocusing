"""Check early-selection configuration against the real executable and runner."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

REPO, BIN, BUILD = (Path(p).resolve() for p in sys.argv[1:4])
del sys.argv[1:4]


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='selection-', dir=BUILD)
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.env = {k: v for k, v in os.environ.items() if not k.startswith('AUTOFOCUSING_')}
        self.env.update(AUTOFOCUSING_BIN=str(BIN), AUTOFOCUSING_BACKEND='cpu', OMP_NUM_THREADS='1')
        (self.work / 'input').mkdir()
        (self.work / 'catalog').write_text('104 001 00 00 0 140 35 10 1e27\n')
        (self.work / 'local_config.sh').write_text('HINET_ROOT=input\nCMT_CATALOG=catalog\n')
        self.config = self.work / 'analysis/trial/config.sh'
        self.config.parent.mkdir(parents=True)
        self.config.write_text('START_DATE=2004-01-01\nEND_DATE=2004-01-01\n')

    def command(self, command):
        return subprocess.run(command, cwd=self.work, env=self.env, capture_output=True, text=True)

    def test_default_and_selected_configuration(self):
        result = self.command([str(BIN), '--event-selection-info'])
        self.assertEqual(result.returncode, 0, result.stderr)
        info = json.loads(result.stdout)
        self.assertEqual(info['mode'], 'all')
        self.assertEqual(info['minimum_max_mad'], dict(R=7, T=7, U=35))
        self.env.update(AUTOFOCUSING_EVENT_SELECTION='selected', AUTOFOCUSING_MIN_MAX_MAD_R='9',
                        AUTOFOCUSING_MIN_MAX_MAD_T='8', AUTOFOCUSING_MIN_MAX_MAD_U='40')
        info = json.loads(self.command([str(BIN), '--event-selection-info']).stdout)
        self.assertEqual(info['mode'], 'selected')
        self.assertEqual(info['minimum_max_mad'], dict(R=9, T=8, U=40))
        result = self.command(['bash', str(REPO / 'Scripts/run.sh'), '--workspace',
                               str(self.work), '--experiment', 'trial'])
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest_file, = (self.work / 'results/trial').glob('*/manifest.json')
        manifest = json.loads(manifest_file.read_text())
        self.assertEqual(manifest['event_selection'], info)
        self.assertEqual(manifest['settings']['AUTOFOCUSING_MIN_MAX_MAD_T'], '8')
        self.assertIn('#EventSelection ' + json.dumps(info, separators=(',', ':')),
                      (manifest_file.parent / 'run.log').read_text())

    def test_invalid_settings_fail_before_start(self):
        cases = [('AUTOFOCUSING_EVENT_SELECTION', v) for v in ('', 'Selected', 'on', 'select')]
        cases += [('AUTOFOCUSING_MIN_MAX_MAD_' + c, v)
                  for c in ('R', 'T', 'U') for v in ('', '0', '-1', 'nan', 'inf', '7junk')]
        for key, value in cases:
            self.env[key] = value
            result = self.command([str(BIN), '--event-selection-info'])
            self.assertNotEqual(result.returncode, 0, (key, value))
            self.assertIn(key, result.stderr)
            result = self.command(['bash', str(REPO / 'Scripts/run.sh'), '--workspace',
                                   str(self.work), '--experiment', 'trial'])
            self.assertNotEqual(result.returncode, 0, (key, value))
            self.assertIn(key, result.stderr)
            self.assertFalse((self.work / 'results').exists())
            del self.env[key]


if __name__ == '__main__':
    unittest.main()
