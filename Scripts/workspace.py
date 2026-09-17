#!/usr/bin/env python3
"""Workspace setup and recorded runs; Python standard library only."""
import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys


DEFAULTS = {
    'COMPONENT_MODE': 'horizontal',
    'AUTOFOCUSING_BACKEND': 'cpu', 'AUTOFOCUSING_SLOWNESS_STEP': '0.005',
    'AUTOFOCUSING_SLOWNESS_MAX': '0.165',
    'AUTOFOCUSING_FREQ_MIN': '0.1', 'AUTOFOCUSING_FREQ_MAX': '0.25',
    'AUTOFOCUSING_EVENT_SELECTION': 'all',
    'AUTOFOCUSING_MIN_MAX_MAD_R': '7', 'AUTOFOCUSING_MIN_MAX_MAD_T': '7',
    'AUTOFOCUSING_MIN_MAX_MAD_U': '35',
}
OPTIONAL_SETTINGS = (
    'OMP_NUM_THREADS', 'OMP_DYNAMIC', 'OMP_PROC_BIND', 'OMP_PLACES',
    'OMP_SCHEDULE', 'OMP_THREAD_LIMIT', 'OMP_STACKSIZE', 'OMP_WAIT_POLICY',
    'CUDA_VISIBLE_DEVICES', 'AUTOFOCUSING_GPU_POWER', 'AUTOFOCUSING_METAL_POWER',
    'AUTOFOCUSING_VERIFY_METAL', 'AUTOFOCUSING_VERIFY_CUDA',
    'AUTOFOCUSING_PROFILE', 'AUTOFOCUSING_CUDA_PROFILE',
)


def experiment_name(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', value) or '..' in value:
        raise ValueError('Invalid experiment name: ' + value)
    return value


def git(path, *args):
    result = subprocess.run(['git', '-C', str(path), *args],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.returncode, result.stdout.decode('utf-8', errors='replace')


def git_root(path):
    code, output = git(path, 'rev-parse', '--show-toplevel')
    return Path(output.strip()).resolve() if code == 0 else None


def create_file(path, content, executable=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open('x') as stream:
            stream.write(content)
    except FileExistsError:
        print('Preserved:', path)
        return
    if executable:
        path.chmod(0o755)
    print('Created:', path)


def setup(args):
    repo, workspace = args.repo, args.workspace
    name = experiment_name(args.experiment)
    if workspace == repo or repo in workspace.parents:
        raise ValueError('Choose a workspace outside the public source repository')
    analysis = workspace / 'analysis'
    workspace.mkdir(parents=True, exist_ok=True)
    # Validate before populating a directory that belongs to another checkout.
    if args.init_git:
        existing_root = git_root(analysis if analysis.exists() else workspace)
        if existing_root is not None and existing_root != analysis.resolve():
            raise ValueError('analysis/ is inside another Git repository; refusing --init-git')
    analysis.mkdir(exist_ok=True)
    (workspace / 'results').mkdir(exist_ok=True)
    create_file(analysis / 'README.md', '''# Analysis workspace

Track experiment settings, procedures, plotting scripts and interpretation here.
This repository is separate from the public source checkout. Each experiment
has its own directory; generated results live in ../results/ and need a separate
backup. Do not commit external waveform data, catalogs or machine-local settings.

Record meaningful configuration changes before formal runs. The runner saves
both source and analysis revisions, configuration snapshots and logs with each run.
''')
    create_file(analysis / '.gitignore',
                '# Generated and machine-local files\nresults/\nlocal_config.sh\n__pycache__/\n.DS_Store\n')
    create_file(analysis / name / 'config.sh',
                (repo / 'Scripts/experiment_config.example.sh').read_text())
    create_file(analysis / name / 'README.md', f'''# {name}

Describe the scientific objective, input selection, procedure and interpretation.
Refer to output runs by their timestamp ID and record comparisons here.

The template retains the current nominal 0.1–0.25 Hz frequency band and
px/py range ±0.165 s/km. An experiment name does not change the calculation.
For primary microseisms, set AUTOFOCUSING_FREQ_MIN=0.05 and
AUTOFOCUSING_FREQ_MAX=0.1 in config.sh; choose dates and slowness separately.
''')
    create_file(workspace / 'local_config.sh',
                (repo / 'Scripts/local_config.example.sh').read_text())
    # Relative to the launcher, so moving the complete workspace still works.
    relative_repo = os.path.relpath(repo, workspace)
    wrapper = '''#!/usr/bin/env bash
set -euo pipefail
WORKSPACE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SOURCE_REPO="$WORKSPACE/"''' + shlex.quote(relative_repo) + '''
exec bash "$SOURCE_REPO/Scripts/run.sh" --workspace "$WORKSPACE" "$@"
'''
    create_file(workspace / 'run.sh', wrapper, executable=True)
    if args.init_git and not (analysis / '.git').exists():
        subprocess.run(['git', 'init', str(analysis)], check=True)
    print('Edit local_config.sh and analysis/' + name + '/config.sh before running.')
    print('Run: bash ' + shlex.quote(str(repo / 'Scripts/run.sh')) +
          ' --workspace ' + shlex.quote(str(workspace)) + ' --experiment ' + name)
    print('Existing parent launchers are preserved; see Scripts/README.md for migration.')


def utc_now():
    return dt.datetime.now(dt.timezone.utc)


def write_json(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True) + '\n')
    temporary.replace(path)


def checksum(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path):
    return {'path': str(path), 'sha256': checksum(path), 'size_bytes': path.stat().st_size}


def provenance(path, output, label):
    root = git_root(path)
    if root != path.resolve():
        return {'path': str(path), 'git_initialized': False,
                'enclosing_repository': str(root) if root else None}
    _, sha = git(path, 'rev-parse', '--verify', 'HEAD')
    sha = sha.strip() or None
    _, branch = git(path, 'symbolic-ref', '--quiet', '--short', 'HEAD')
    _, status = git(path, 'status', '--porcelain=v1', '--untracked-files=all')
    if sha:
        _, diff = git(path, 'diff', '--binary', '--no-ext-diff', '--no-textconv', 'HEAD', '--')
    else:
        _, staged = git(path, 'diff', '--cached', '--binary', '--no-ext-diff', '--no-textconv')
        _, unstaged = git(path, 'diff', '--binary', '--no-ext-diff', '--no-textconv')
        diff = staged + unstaged
    _, untracked = git(path, 'ls-files', '--others', '--exclude-standard', '-z')
    (output / (label + '.patch')).write_text(diff)
    write_json(output / (label + '-untracked.json'), untracked.rstrip('\0').split('\0') if untracked else [])
    (output / (label + '-status.txt')).write_text(status)
    return {'path': str(path), 'git_initialized': True, 'commit': sha,
            'branch': branch.strip() or None, 'dirty': bool(status),
            'patch': label + '.patch', 'status_file': label + '-status.txt',
            'untracked_files': label + '-untracked.json'}


def resolve_path(workspace, value):
    path = Path(value).expanduser()
    return (workspace / path).resolve() if not path.is_absolute() else path.resolve()


def allocate_run(root):
    root.mkdir(parents=True, exist_ok=True)
    stamp = utc_now().strftime('%Y%m%dT%H%M%SZ')
    for suffix in range(1000000):
        candidate = root / (stamp if suffix == 0 else f'{stamp}-{suffix:02d}')
        try:
            candidate.mkdir()
            return candidate
        except FileExistsError:
            continue
    raise RuntimeError('Could not allocate a unique run directory')


def run(args):
    workspace, repo = args.workspace, args.repo
    name = experiment_name(args.experiment)
    settings = {key: os.environ.get(key) or default for key, default in DEFAULTS.items()}
    # An explicitly empty backend/grid setting is invalid to C++, not a default.
    for key in DEFAULTS:
        if key.startswith('AUTOFOCUSING_') and key in os.environ:
            settings[key] = os.environ[key]
    settings.update({key: os.environ[key] for key in OPTIONAL_SETTINGS if key in os.environ})
    if 'AUTOFOCUSING_GPU_POWER' not in settings:
        settings['AUTOFOCUSING_GPU_POWER'] = (
            settings.get('AUTOFOCUSING_METAL_POWER', 'off')
            if settings['AUTOFOCUSING_BACKEND'] == 'metal' else 'off')
    if settings['COMPONENT_MODE'] not in ('horizontal', '3c'):
        raise ValueError('COMPONENT_MODE must be horizontal or 3c')
    if settings['AUTOFOCUSING_EVENT_SELECTION'] not in ('all', 'selected'):
        raise ValueError('AUTOFOCUSING_EVENT_SELECTION must be all or selected')
    selection = {'mode': settings['AUTOFOCUSING_EVENT_SELECTION'],
                 'score': 'initial_grid_max_over_mad', 'comparison': '>',
                 'minimum_max_mad': {}}
    for component in ('R', 'T', 'U'):
        key = 'AUTOFOCUSING_MIN_MAX_MAD_' + component
        try:
            value = float(settings[key])
        except ValueError:
            raise ValueError(key + ' must be a positive finite number') from None
        if not math.isfinite(value) or value <= 0:
            raise ValueError(key + ' must be a positive finite number')
        selection['minimum_max_mad'][component] = value
    dates = {}
    for key in ('START_DATE', 'END_DATE'):
        value = os.environ.get(key, '')
        if not re.fullmatch(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', value):
            raise ValueError('Set ' + key + ' in experiment config.sh using YYYY-MM-DD')
        try:
            dates[key] = dt.date.fromisoformat(value)
        except ValueError:
            raise ValueError(key + ' is not a valid calendar date: ' + value) from None
        if dates[key].year < 1400:
            raise ValueError(key + ' must be in the supported year range 1400–9999')
        settings[key] = value
    if dates['END_DATE'] < dates['START_DATE']:
        raise ValueError('END_DATE must be on or after START_DATE')
    settings['START_YEAR'] = str(dates['START_DATE'].year)
    for key, is_directory in [('HINET_ROOT', True), ('CMT_CATALOG', False)]:
        value = os.environ.get(key)
        if not value:
            raise ValueError('Set ' + key + ' in local_config.sh or the environment')
        path = resolve_path(workspace, value)
        if not (path.is_dir() if is_directory else path.is_file()):
            raise ValueError(key + ' not found: ' + str(path))
        settings[key] = str(path)
    binary = os.environ.get('AUTOFOCUSING_BIN')
    if binary:
        binary = resolve_path(workspace, binary)
    else:
        binary = next((p for p in (repo / 'bin/cal_ccf_clang', repo / 'bin/cal_ccf_gcc')
                       if p.is_file() and os.access(p, os.X_OK)), None)
    if binary is None or not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError('Build and install cal_ccf first, or set AUTOFOCUSING_BIN')
    results = resolve_path(workspace, os.environ.get('RESULTS_ROOT') or 'results')
    settings.update(AUTOFOCUSING_BIN=str(binary), RESULTS_ROOT=str(results), PARAM_ID=name)
    # Every recorded setting is passed to the child, including normalized defaults.
    env = dict(os.environ, **settings)
    info = subprocess.run([str(binary), '--frequency-info'], cwd=repo, env=env,
                          capture_output=True, text=True, timeout=30)
    if info.returncode != 0:
        raise ValueError('Cannot read frequency configuration; use a rebuilt cal_ccf with '
                         '--frequency-info support.\n' + info.stderr.strip())
    try:
        frequency_band = json.loads(info.stdout)
        if not isinstance(frequency_band, dict) or not all(key in frequency_band for key in
                ('requested_min_hz', 'requested_max_hz', 'min_hz', 'max_hz', 'df_hz', 'min_bin', 'max_bin')):
            raise ValueError('Missing frequency fields')
    except (ValueError, TypeError):
        raise ValueError('Invalid --frequency-info response; rebuild the selected cal_ccf') from None
    # Legacy binaries still support the default full-catalog behavior. Require
    # an explicit capability response before enabling early selection.
    if selection['mode'] == 'selected':
        info = subprocess.run([str(binary), '--event-selection-info'], cwd=repo, env=env,
                              capture_output=True, text=True, timeout=30)
        if info.returncode != 0:
            raise ValueError('Early event selection requires a rebuilt cal_ccf with '
                             '--event-selection-info support.\n' + info.stderr.strip())
        try:
            reported_selection = json.loads(info.stdout)
        except ValueError:
            raise ValueError('Invalid --event-selection-info response; rebuild cal_ccf') from None
        if reported_selection != selection:
            raise ValueError('The executable did not confirm the requested event selection settings')
    output = allocate_run(results / name)
    command = [str(binary), settings['START_YEAR'], name, output.name,
               settings['HINET_ROOT'], settings['CMT_CATALOG'], str(results), settings['COMPONENT_MODE'],
               settings['START_DATE'], settings['END_DATE']]
    manifest = {'schema_version': 1, 'experiment': name, 'run_id': output.name,
                'started_at': utc_now().isoformat(), 'finished_at': None,
                'status': 'running', 'exit_code': None, 'command': command,
                'cwd': str(repo), 'workspace': str(workspace), 'settings': settings,
                'frequency_band': frequency_band, 'event_selection': selection}
    manifest_path = output / 'manifest.json'
    write_json(manifest_path, manifest)
    child = None
    received_signal = None

    def stop(signum, _frame):
        nonlocal received_signal
        received_signal = signum
        if child is not None and child.poll() is None:
            try:
                os.killpg(child.pid, signum)
            except ProcessLookupError:
                pass

    previous_handlers = {sig: signal.signal(sig, stop) for sig in (signal.SIGINT, signal.SIGTERM)}
    code = 1
    try:
        manifest['source'] = provenance(repo, output, 'source')
        analysis = workspace / 'analysis'
        manifest['analysis'] = provenance(analysis, output, 'analysis') if analysis.is_dir() else {
            'path': str(analysis), 'git_initialized': False}
        manifest['binary'] = file_identity(binary)
        manifest['binary']['source_revision_verified'] = False
        manifest['catalog'] = file_identity(Path(settings['CMT_CATALOG']))
        manifest['waveforms'] = {'root': settings['HINET_ROOT'], 'content_hashed': False}
        model = repo / 'vel_Nishida2008.dat'
        if model.is_file():
            manifest['velocity_model'] = file_identity(model)
        manifest['config_snapshots'] = []
        configs = [('local_config.sh', workspace / 'local_config.sh'),
                   ('experiment_config.sh', analysis / name / 'config.sh')]
        for filename, source in configs:
            if source.is_file():
                (output / filename).write_bytes(source.read_bytes())
                manifest['config_snapshots'].append({'source': str(source), 'snapshot': filename})
        write_json(output / 'effective_config.json', settings)
        write_json(manifest_path, manifest)
        print('Output: ' + str(output), flush=True)
        with (output / 'run.log').open('wb') as log:
            if received_signal is None:
                child = subprocess.Popen(command, cwd=repo, env=env, stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT, start_new_session=True)
                if received_signal is not None:
                    stop(received_signal, None)
                for block in iter(lambda: child.stdout.read1(65536), b''):
                    log.write(block)
                    log.flush()
                    try:
                        sys.stdout.buffer.write(block)
                        sys.stdout.buffer.flush()
                    except BrokenPipeError:
                        # Closing a terminal consumer must not lose the saved log.
                        pass
                child.stdout.close()
                code = child.wait()
                if code < 0:
                    code = 128 - code
        if received_signal is not None:
            code = 128 + received_signal
        manifest['status'] = 'interrupted' if received_signal else ('succeeded' if code == 0 else 'failed')
    except Exception as error:
        code = 1
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        manifest['status'] = 'failed'
        manifest['error'] = str(error)
        print('Run failed: ' + str(error), file=sys.stderr)
    finally:
        manifest.update(exit_code=code, finished_at=utc_now().isoformat(),
                        signal=received_signal)
        write_json(manifest_path, manifest)
        for sig, handler in previous_handlers.items():
            signal.signal(sig, handler)
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for command in ('setup', 'run'):
        sub = commands.add_parser(command)
        sub.add_argument('--repo', type=Path, required=True)
        sub.add_argument('--workspace', type=Path)
        sub.add_argument('--experiment', required=True)
        if command == 'setup':
            sub.add_argument('--init-git', action='store_true')
    args = parser.parse_args()
    args.repo = args.repo.resolve()
    args.workspace = args.workspace.resolve() if args.workspace else args.repo.parent
    try:
        return setup(args) if args.command == 'setup' else run(args)
    except (ValueError, OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
