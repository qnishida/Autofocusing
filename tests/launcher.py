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
    result = subprocess.run(['bash', str(workspace / 'run.sh')], cwd=tempfile.gettempdir(),
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    files = list((workspace / 'results/launcher_test').glob('*/*.dat'))
    assert len(files) == 1 and files[0].stat().st_size == 0, result.stdout
    legacy = subprocess.run([str(executable), '2004', 'legacy', 'test', 'input', 'catalog'],
                            cwd=workspace, capture_output=True, text=True)
    assert legacy.returncode == 0, legacy.stderr
    assert len(list((workspace / 'output/legacy/test').glob('*.dat'))) == 1
    invalid = subprocess.run([str(executable), '2004', 'invalid', 'test', 'input',
                              'catalog', 'results', 'bad-mode'],
                             cwd=workspace, capture_output=True, text=True)
    assert invalid.returncode != 0 and 'Component mode' in invalid.stderr
    (workspace / 'catalog').unlink()
    result = subprocess.run(['bash', str(workspace / 'run.sh')], cwd=tempfile.gettempdir(),
                            env=env, capture_output=True, text=True)
    assert result.returncode != 0 and 'CMT_CATALOG' in result.stderr
print('PASS: parent launcher resolves paths from another cwd; missing catalog fails clearly')
