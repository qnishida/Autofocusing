"""Exercise qualification ordering, CPU comparison and rejection without a GPU."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('power_runner', HERE/'run_metal_power.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class QualificationTests(unittest.TestCase):
    def test_tolerances_and_nonfinite(self):
        ref = [['1']*38]
        row = [ref[0].copy()]
        row[0][20] = '1.0001'
        runner.compare_power(row, ref, 'bootstrap')
        with self.assertRaises(RuntimeError):
            runner.compare_power(row, ref, 'off', bootstrap=False)
        row[0][20] = 'nan'
        with self.assertRaises(RuntimeError):
            runner.compare_power(row, ref, 'bootstrap')
        row[0][20] = '1'
        row[0][5] = '2'
        with self.assertRaises(RuntimeError):
            runner.compare_power(row, ref, 'all')

    def test_runner_pairs_off_first_and_reports(self):
        with tempfile.TemporaryDirectory(prefix='gpu-tools-') as temp:
            root = Path(temp)
            for day in range(1, 4):
                file = root/'input/2004'/f'010{day}'/f'200400{day}0000.h5'
                file.parent.mkdir(parents=True)
                file.touch()
            (root/'catalog').touch()
            fake = root/'driver'
            fake.write_text('#!' + sys.executable + '\n' + '''
import os, sys
from pathlib import Path
tag=sys.argv[3]
out=Path(sys.argv[6])/sys.argv[2]/tag
out.mkdir(parents=True)
row=['1']*38
if os.environ['AUTOFOCUSING_GPU_POWER'] in ('bootstrap','all'): row[20]='1.00001'
(out/'events.dat').write_text(' '.join(row)+'\\n')
print('#deg max=1')
print('#POWER_PROFILE stage=grid total_s=1')
print('#POWER_PROFILE stage=bootstrap total_s=2')
print('#PROFILE segment=0 windows=1 stack_s=1')
''')
            fake.chmod(0o755)
            command = [sys.executable, str(HERE/'run_metal_power.py'), str(fake),
                       str(root/'input'), str(root/'catalog'), str(root/'output'),
                       '--backend', 'cuda', '--threads', '2', '--repeats', '2',
                       '--compare-cpu', '--modes', 'all', 'off', 'bootstrap', 'grid']
            env = dict(os.environ, AUTOFOCUSING_GPU_POWER='invalid',
                       AUTOFOCUSING_VERIFY_CUDA='1', AUTOFOCUSING_CUDA_PROFILE='1')
            result = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            report = json.loads((root/'output/report.json').read_text())
            self.assertEqual(report['backend'], 'cuda')
            self.assertEqual(report['threads'], 2)
            self.assertEqual([r['mode'] for r in report['runs'][:2]], ['cpu', 'off'])
            self.assertEqual(report['runs'][-1]['mode'], 'cpu')
            self.assertEqual(len(report['runs']), 10)
            self.assertFalse(report['summary']['all']['stage_pass'])
            self.assertTrue(all(r['max_rss_bytes'] is not None for r in report['runs']))


if __name__ == '__main__':
    unittest.main()
