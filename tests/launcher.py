"""Exercise the deployed parent launcher from another cwd using an empty archive."""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import tempfile

repo, executable, build = (Path(p).resolve() for p in sys.argv[1:])
with tempfile.TemporaryDirectory(prefix='launcher-', dir=build) as tmp:
    workspace = Path(tmp)
    (workspace / 'repo').symlink_to(repo, target_is_directory=True)
    shutil.copy2(repo / 'Scripts/run.sh', workspace / 'run.sh')
    (workspace / 'input').mkdir()
    (workspace / 'catalog').write_text('104 001 00 00 0 140 35 10 1e27\n')
    (workspace / 'local_config.sh').write_text(
        'HINET_ROOT=input\nCMT_CATALOG=catalog\nCOMPONENT_MODE=horizontal\n'
        'PARAM_ID=launcher_test\nSTART_YEAR=2004\nRESULTS_ROOT=results\n')
    env = os.environ.copy()
    env.pop('START_DATE', None)
    env.pop('END_DATE', None)
    env['AUTOFOCUSING_BIN'] = str(executable)
    env['AUTOFOCUSING_BACKEND'] = 'cpu'
    env['AUTOFOCUSING_GPU_POWER'] = 'off'
    for name in ('launcher_test', 'relative_bin_test'):
        config = workspace / 'analysis' / name / 'config.sh'
        config.parent.mkdir(parents=True)
        config.write_text('COMPONENT_MODE=horizontal\nSTART_DATE=2004-02-28\nEND_DATE=2004-03-01\n')
    result = subprocess.run(['bash', str(workspace / 'run.sh')], cwd=tempfile.gettempdir(),
                            env=env, capture_output=True, text=True)
    assert result.returncode == 2 and 'Example:' in result.stderr
    assert not (workspace / 'results').exists()
    result = subprocess.run(['bash', str(workspace / 'run.sh'), '--experiment', 'launcher_test'],
                            cwd=tempfile.gettempdir(), env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    files = list((workspace / 'results/launcher_test').glob('*/*.dat'))
    assert len(files) == 1 and files[0].stat().st_size == 0, result.stdout
    manifest = json.loads((files[0].parent / 'manifest.json').read_text())
    assert manifest['status'] == 'succeeded' and manifest['exit_code'] == 0
    assert manifest['run_id'] == files[0].parent.name
    assert manifest['source']['commit'] and not manifest['binary']['source_revision_verified']
    assert manifest['settings']['START_DATE'] == '2004-02-28'
    assert manifest['settings']['END_DATE'] == '2004-03-01'
    assert manifest['command'][-2:] == ['2004-02-28', '2004-03-01']
    assert '#ScanSummary days=3 missing=3 empty=0 multiple=0 loaded=0' in result.stdout
    assert (files[0].parent / 'run.log').is_file()
    # A copied, editable launcher also resolves a configured executable path
    # relative to its workspace, even when invoked from another cwd.
    assert not (workspace / 'run.sh').is_symlink()
    (workspace / 'bin').mkdir()
    (workspace / 'bin/cal_ccf').symlink_to(executable)
    with (workspace / 'local_config.sh').open('a') as config:
        config.write('AUTOFOCUSING_BIN=bin/cal_ccf\nPARAM_ID=relative_bin_test\n')
    result = subprocess.run(['bash', str(workspace / 'run.sh'), '--experiment', 'relative_bin_test'], cwd=tempfile.gettempdir(),
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert len(list((workspace / 'results/relative_bin_test').glob('*/*.dat'))) == 1
    legacy = subprocess.run([str(executable), '2004', 'legacy', 'test', 'input', 'catalog'],
                            cwd=workspace, capture_output=True, text=True)
    assert legacy.returncode == 0, legacy.stderr
    assert len(list((workspace / 'output/legacy/test').glob('*.dat'))) == 1
    invalid = subprocess.run([str(executable), '2004', 'invalid', 'test', 'input',
                              'catalog', 'results', 'bad-mode'],
                             cwd=workspace, capture_output=True, text=True)
    assert invalid.returncode != 0 and 'Component mode' in invalid.stderr
    for key, value, message in [
        ('AUTOFOCUSING_BACKEND', 'invalid', 'AUTOFOCUSING_BACKEND'),
        ('AUTOFOCUSING_GPU_POWER', 'invalid', 'AUTOFOCUSING_GPU_POWER'),
        ('AUTOFOCUSING_SLOWNESS_STEP', '0', 'positive finite'),
        ('AUTOFOCUSING_SLOWNESS_MAX', 'nan', 'positive finite'),
        ('AUTOFOCUSING_SLOWNESS_STEP', '.005junk', 'positive finite'),
        ('AUTOFOCUSING_SLOWNESS_STEP', '1e-10', 'half-width'),
    ]:
        invalid_env = env.copy()
        invalid_env[key] = value
        invalid = subprocess.run([str(executable)], env=invalid_env,
                                 capture_output=True, text=True)
        assert invalid.returncode != 0 and message in invalid.stderr, invalid.stderr
    (workspace / 'catalog').unlink()
    result = subprocess.run(['bash', str(workspace / 'run.sh'), '--experiment', 'launcher_test'], cwd=tempfile.gettempdir(),
                            env=env, capture_output=True, text=True)
    assert result.returncode != 0 and 'CMT_CATALOG' in result.stderr
print('PASS: parent launcher resolves paths from another cwd; missing catalog fails clearly')
