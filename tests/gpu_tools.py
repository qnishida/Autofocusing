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
    def test_cpu_runner_three_components_across_year_boundary(self):
        with tempfile.TemporaryDirectory(prefix='cpu-bootstrap-') as temp:
            root = Path(temp)
            for year, day in ((2004, '1231'), (2005, '0101'), (2005, '0102')):
                file = root/'input'/str(year)/day/'waveforms.h5'
                file.parent.mkdir(parents=True)
                file.touch()
            (root/'catalog').touch()
            fake = root/'driver'
            fake.write_text('#!' + sys.executable + '\n' + '''
import sys
from pathlib import Path
assert sys.argv[1] == '2004'
assert sys.argv[7:] == ['3c', '2004-12-31', '2005-01-02']
for year, day in ((2004, '1231'), (2005, '0101'), (2005, '0102')):
    assert (Path(sys.argv[4])/str(year)/day/'waveforms.h5').is_file()
out = Path(sys.argv[6])/sys.argv[2]/sys.argv[3]
out.mkdir(parents=True)
(out/'events.dat').write_text(' '.join(['1']*38)+'\\n')
print('#deg max=1')
print('#POWER_PROFILE stage=bootstrap total_s=1')
for day in range(3):
    print('#LOAD_PROFILE init_s=1 read_decode_copy_s=1 filter_s=1 total_s=3')
    for segment in range(4):
        print(f'#PROFILE segment={segment} windows=1 fft_qc_s=1 pack_s=1 stack_s=1 fit_s=1 total_s=4')
for stage in ('objective', 'hessian', 'rotation', 'fitting'):
    for context in ('serial_caller', 'outer_parallel'):
        print(f'#CPU_PROFILE stage={stage} context={context} calls=1 inclusive_s=1')
''')
            fake.chmod(0o755)
            command = [sys.executable, str(HERE/'run_cpu_parallel_events.py'), str(fake), str(fake),
                       str(root/'input'), str(root/'catalog'), str(root/'output'),
                       '--backend', 'cpu', '--power', 'off', '--components', '3c',
                       '--start-date', '2004-12-31', '--days', '3', '--expected-events', '1',
                       '--repeats', '2', '--warmups', '0', '--target', 'bootstrap']
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            report = json.loads((root/'output/report.json').read_text())
            self.assertEqual(report['components'], '3c')
            self.assertEqual(report['end_date'], '2005-01-02')
            self.assertTrue(report['all_event_bytes_identical'])
            self.assertFalse(report['target_stage_pass'])
            self.assertEqual([r['mode'] for r in report['runs']],
                             ['reference', 'current', 'current', 'reference'])

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
