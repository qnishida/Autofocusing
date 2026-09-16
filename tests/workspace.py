"""Integration checks for workspace setup and recorded runs, without seismic data."""
import concurrent.futures
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest

REPO = Path(__file__).resolve().parents[1]


class WorkspaceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='autofocusing workspace ')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.repo = self.base / 'source checkout'
        self.repo.mkdir()
        shutil.copytree(REPO / 'Scripts', self.repo / 'Scripts')
        (self.repo / 'src').mkdir()
        (self.repo / 'vel_Nishida2008.dat').write_text('model fixture\n')
        self.workspace = self.base / 'analysis workspace'
        self.env = {k: v for k, v in os.environ.items()
                    if not k.startswith(('AUTOFOCUSING_', 'OMP_', 'GIT_')) and
                    k not in ('PARAM_ID', 'START_YEAR', 'HINET_ROOT', 'CMT_CATALOG',
                              'COMPONENT_MODE', 'RESULTS_ROOT', 'EXPERIMENT', 'WORK_DIR', 'REPO_DIR')}
        self.setup()
        (self.workspace / 'input data').mkdir()
        (self.workspace / 'catalog').write_text('catalog fixture\n')
        binary = self.workspace / 'mock program'
        binary.write_text('''#!/usr/bin/env python3
import os, pathlib, sys, time
root = pathlib.Path(sys.argv[6]) / sys.argv[2] / sys.argv[3]
(root / 'result.dat').write_text(os.environ['AUTOFOCUSING_SLOWNESS_MAX'])
print('mock stdout', flush=True)
print('mock stderr', file=sys.stderr, flush=True)
if os.environ.get('MOCK_WAIT'):
    print('ready', flush=True)
    time.sleep(60)
sys.exit(int(os.environ.get('MOCK_EXIT', '0')))
''')
        binary.chmod(0o755)
        self.local = self.workspace / 'local_config.sh'
        self.local.write_text("HINET_ROOT='input data'\nCMT_CATALOG=catalog\n"
                              "AUTOFOCUSING_BIN='mock program'\nAUTOFOCUSING_SLOWNESS_MAX=0.2\n")
        self.config = self.workspace / 'analysis/trial/config.sh'

    def command(self, args, **kwargs):
        return subprocess.run(args, cwd=self.base, env=self.env,
                              capture_output=True, text=True, **kwargs)

    def setup(self, *extra):
        result = self.command(['bash', str(self.repo / 'Scripts/setup_workspace.sh'),
                               '--workspace', str(self.workspace), '--experiment', 'trial', *extra])
        self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def run_command(self, experiment=True):
        args = ['bash', str(self.workspace / 'run.sh')]
        return args + (['--experiment', 'trial'] if experiment else [])

    def outputs(self, name='trial'):
        return sorted((self.workspace / 'results' / name).glob('*/manifest.json'))

    def initialize(self, path):
        for args in [('init',), ('add', '.'),
                     ('-c', 'user.name=Test', '-c', 'user.email=test@example.invalid',
                      '-c', 'commit.gpgsign=false', 'commit', '-m', 'fixture')]:
            result = self.command(['git', '-C', str(path), *args])
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_setup_preserves_and_initializes(self):
        originals = {p: p.read_bytes() for p in (self.local, self.config, self.workspace / 'run.sh')}
        self.setup('--init-git')
        self.setup('--init-git')
        self.assertTrue((self.workspace / 'analysis/.git').is_dir())
        self.assertFalse((self.workspace / '.git').exists())
        for path, content in originals.items():
            self.assertEqual(path.read_bytes(), content)
        result = self.command(['bash', str(self.repo / 'Scripts/setup_workspace.sh'),
                               '--workspace', str(self.repo), '--experiment', 'trial', '--init-git'])
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.command(['git', 'init', str(self.workspace)]).returncode, 0)
        # A fresh nested analysis directory must not inherit the parent's Git root.
        nested = self.workspace / 'nested'
        result = self.command(['bash', str(self.repo / 'Scripts/setup_workspace.sh'),
                               '--workspace', str(nested), '--experiment', 'trial', '--init-git'])
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((nested / 'analysis').exists())

    def test_success_precedence_and_git_provenance(self):
        self.env['AUTOFOCUSING_SLOWNESS_MAX'] = '0.3'
        self.env['UNRELATED_SECRET'] = 'must-not-be-recorded'
        self.initialize(self.repo)
        self.initialize(self.workspace / 'analysis')
        self.config.write_text(self.config.read_text() + 'AUTOFOCUSING_SLOWNESS_MAX=0.4\n')
        (self.workspace / 'analysis/notes.txt').write_text('untracked notes')
        result = self.command(self.run_command())
        self.assertEqual(result.returncode, 0, result.stderr)
        output = self.outputs()[0].parent
        manifest = json.loads((output / 'manifest.json').read_text())
        self.assertEqual(manifest['status'], 'succeeded')
        self.assertEqual(manifest['settings']['AUTOFOCUSING_SLOWNESS_MAX'], '0.4')
        self.assertEqual(manifest['settings']['AUTOFOCUSING_BACKEND'], 'cpu')
        self.assertEqual((output / 'result.dat').read_text(), '0.4')
        self.assertEqual(len(manifest['source']['commit']), 40)
        self.assertFalse(manifest['source']['dirty'])
        self.assertTrue(manifest['analysis']['dirty'])
        self.assertIn('AUTOFOCUSING_SLOWNESS_MAX=0.4', (output / 'analysis.patch').read_text())
        self.assertIn('notes.txt', json.loads((output / 'analysis-untracked.json').read_text()))
        self.assertEqual((output / 'experiment_config.sh').read_bytes(), self.config.read_bytes())
        self.assertIn('mock stdout', (output / 'run.log').read_text())
        self.assertIn('mock stderr', result.stdout)
        self.assertNotIn('UNRELATED_SECRET', (output / 'manifest.json').read_text())
        self.assertFalse(manifest['binary']['source_revision_verified'])
        self.assertTrue(manifest['finished_at'])

    def test_repeat_and_parallel_runs(self):
        for _ in range(2):
            self.assertEqual(self.command(self.run_command()).returncode, 0)
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            results = list(pool.map(lambda _: self.command(self.run_command()), range(4)))
        self.assertTrue(all(result.returncode == 0 for result in results))
        self.assertEqual(len(self.outputs()), 6)
        for path in self.outputs():
            manifest = json.loads(path.read_text())
            self.assertFalse(manifest['analysis']['git_initialized'])
            self.assertEqual(manifest['run_id'], path.parent.name)

    def test_timestamp_collision(self):
        spec = importlib.util.spec_from_file_location('workspace_helper', REPO / 'Scripts/workspace.py')
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        fixed = helper.utc_now()
        helper.utc_now = lambda: fixed
        first = helper.allocate_run(self.base / 'collision')
        second = helper.allocate_run(self.base / 'collision')
        self.assertEqual(second.name, first.name + '-01')

    def test_unborn_repository_and_power_alias(self):
        self.setup('--init-git')
        self.env.update(AUTOFOCUSING_BACKEND='metal', AUTOFOCUSING_METAL_POWER='bootstrap')
        result = self.command(self.run_command())
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = json.loads(self.outputs()[0].read_text())
        self.assertTrue(manifest['analysis']['git_initialized'])
        self.assertIsNone(manifest['analysis']['commit'])
        self.assertTrue(manifest['analysis']['dirty'])
        self.assertEqual(manifest['settings']['AUTOFOCUSING_GPU_POWER'], 'bootstrap')
        self.env['AUTOFOCUSING_GPU_POWER'] = 'off'
        before = set(self.outputs())
        result = self.command(self.run_command())
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = json.loads((set(self.outputs()) - before).pop().read_text())
        self.assertEqual(manifest['settings']['AUTOFOCUSING_GPU_POWER'], 'off')

    def test_legacy_failure_and_missing_inputs(self):
        self.env.update(PARAM_ID='legacy', MOCK_EXIT='7')
        result = self.command(self.run_command(experiment=False))
        self.assertEqual(result.returncode, 7, result.stderr)
        manifest = json.loads(self.outputs('legacy')[0].read_text())
        self.assertEqual(manifest['status'], 'failed')
        self.assertEqual(manifest['exit_code'], 7)
        self.assertEqual(manifest['settings']['AUTOFOCUSING_SLOWNESS_MAX'], '0.2')
        (self.workspace / 'catalog').unlink()
        result = self.command(self.run_command())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('CMT_CATALOG', result.stderr)
        self.assertEqual(self.outputs(), [])
        for name in ('../escape', 'nested/name', '..', '-option'):
            result = self.command(['bash', str(self.workspace / 'run.sh'), '--experiment', name])
            self.assertNotEqual(result.returncode, 0)

    def test_interruption(self):
        self.env['MOCK_WAIT'] = '1'
        for sig in (signal.SIGTERM, signal.SIGINT):
            before = set(self.outputs())
            with tempfile.TemporaryFile() as terminal:
                process = subprocess.Popen(self.run_command(), cwd=self.base, env=self.env,
                                           stdout=terminal, stderr=terminal)
                try:
                    deadline = time.monotonic() + 10
                    output = None
                    while time.monotonic() < deadline:
                        paths = set(self.outputs()) - before
                        if paths:
                            output = next(iter(paths)).parent
                            log = output / 'run.log'
                            if log.exists() and 'ready' in log.read_text():
                                break
                        time.sleep(0.02)
                    else:
                        self.fail('Mock process did not start')
                    process.send_signal(sig)
                    self.assertEqual(process.wait(timeout=10), 128 + sig)
                    manifest = json.loads((output / 'manifest.json').read_text())
                    self.assertEqual(manifest['status'], 'interrupted')
                    self.assertEqual(manifest['signal'], sig)
                    self.assertTrue(manifest['finished_at'])
                finally:
                    if process.poll() is None:
                        process.kill()
                        process.wait()


if __name__ == '__main__':
    unittest.main()
