"""Exercise the deployed parent launcher from another cwd using an empty archive."""
from pathlib import Path
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
    env['AUTOFOCUSING_BIN'] = str(executable)
    env['AUTOFOCUSING_BACKEND'] = 'cpu'
    env['AUTOFOCUSING_GPU_POWER'] = 'off'
    result = subprocess.run(['bash', str(workspace / 'run.sh')], cwd=tempfile.gettempdir(),
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    files = list((workspace / 'results/launcher_test').glob('*/*.dat'))
    assert len(files) == 1 and files[0].stat().st_size == 0, result.stdout
    # A copied, editable launcher also resolves a configured executable path
    # relative to its workspace, even when invoked from another cwd.
    assert not (workspace / 'run.sh').is_symlink()
    (workspace / 'bin').mkdir()
    (workspace / 'bin/cal_ccf').symlink_to(executable)
    with (workspace / 'local_config.sh').open('a') as config:
        config.write('AUTOFOCUSING_BIN=bin/cal_ccf\nPARAM_ID=relative_bin_test\n')
    result = subprocess.run(['bash', str(workspace / 'run.sh')], cwd=tempfile.gettempdir(),
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
    result = subprocess.run(['bash', str(workspace / 'run.sh')], cwd=tempfile.gettempdir(),
                            env=env, capture_output=True, text=True)
    assert result.returncode != 0 and 'CMT_CATALOG' in result.stderr
print('PASS: parent launcher resolves paths from another cwd; missing catalog fails clearly')
