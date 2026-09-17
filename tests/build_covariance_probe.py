"""Build a deterministic covariance probe from the saved peak-search audit.

Requires build/peak-search-probe/fixed/cal_ccf.cpp and its documented audit
dependencies. Replaces only the bootstrap function with current production
code plus its covariance helper, retaining the existing deterministic seed.
"""
import argparse
import hashlib
import json
import shlex
import subprocess
from pathlib import Path


def bootstrap(source):
    start = source.rindex('static int est_dist_boot(')
    end = source.index('\n/**', start)
    return source[start:end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    args = parser.parse_args()
    build = args.build.resolve()
    repo = Path(__file__).resolve().parents[1]
    prior = build/'peak-search-probe/fixed/cal_ccf.cpp'
    base = prior.read_text()
    current = (repo/'src/cal_ccf.cpp').read_text()
    old = subprocess.check_output(['git', 'show', '485d761:src/cal_ccf.cpp'], cwd=repo, text=True)
    random_seed = 'static_cast<unsigned long>(time(0) * i)'
    fixed_seed = 'static_cast<unsigned long>(1837u + 104729u * i)'
    assert bootstrap(base) == bootstrap(old).replace(random_seed, fixed_seed)
    assert base[base.rindex('int search_max('):] == current[current.rindex('int search_max('):]
    helper_start = current.index('static void set_parameter_covariance(')
    helper = current[helper_start:current.index('\n/**', helper_start)]
    source = base.replace(bootstrap(base), helper+'\n'+bootstrap(current).replace(random_seed, fixed_seed))

    out = build/'covariance-double-probe'
    out.mkdir(exist_ok=True)
    generated = out/'cal_ccf.cpp'
    generated.write_text(source)
    entry = next(e for e in json.loads((build/'compile_commands.json').read_text())
                 if e['file'].endswith('/cal_ccf.cpp'))
    command = shlex.split(entry['command'])
    original_object = Path(command[command.index('-o')+1])
    obj = out/'cal_ccf.cpp.o'
    command[command.index('-o')+1] = str(obj)
    command[command.index('-c')+1] = str(generated)
    command += ['-I'+str(repo/'src')]
    subprocess.run(command, cwd=entry['directory'], check=True)
    link = shlex.split((Path(entry['directory'])/original_object.parent/'link.txt').read_text())
    objects = {'cal_ccf.cpp.o':str(obj), 'station_info.cpp.o':str(build/'io-fixed-current/station_info.cpp.o')}
    link = [objects.get(Path(s).name, s) for s in link]
    exe = out/'cal_ccf_covariance'
    link[link.index('-o')+1] = str(exe)
    subprocess.run(link, cwd=entry['directory'], check=True)
    files = [prior, repo/'src/cal_ccf.cpp', generated, exe, Path(__file__)]
    (out/'build.json').write_text(json.dumps({
        'compile':command, 'link':link,
        'sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    }, indent=2)+'\n')
    print(exe)


if __name__ == '__main__':
    main()
