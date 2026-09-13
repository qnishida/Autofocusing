"""Build a test-only deterministic event driver using an existing CMake build.

Run from a CPU, Metal or CUDA checkout after configuring CMAKE_EXPORT_COMPILE_COMMANDS.
--reference selects the pre-I/O commit for cal_ccf.cpp and station_info.cpp;
other objects/headers must match that backend. Production sources are untouched.
"""
import argparse
import json
from pathlib import Path
import re
import shlex
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('build', type=Path)
p.add_argument('--reference', help='51f60a7 for CPU, 27fccd3 for Metal')
a = p.parse_args()
repo = Path(__file__).resolve().parents[1]
build = a.build.resolve()
out = build / ('io-fixed-reference' if a.reference else 'io-fixed-current')
out.mkdir(exist_ok=True)
entries = json.loads((build / 'compile_commands.json').read_text())
driver = next(e for e in entries if e['file'].endswith('/cal_ccf.cpp'))
original = Path(shlex.split(driver['command'])[shlex.split(driver['command']).index('-o') + 1])
objects = {}
for name in ('cal_ccf.cpp', 'station_info.cpp'):
    source = (subprocess.check_output(['git', 'show', f'{a.reference}:src/{name}'], cwd=repo, text=True)
              if a.reference else (repo / 'src' / name).read_text())
    if name == 'cal_ccf.cpp':
        if source.count('time(0) * i') != 1:
            raise RuntimeError('Unexpected bootstrap seed sites')
        source = source.replace('time(0) * i', '1837u + 104729u * i')
    else:
        source, count = re.subn(r'(?m)^(\s*plan\s*=\s*fftw_plan_r2r_1d\([^\n]*)FFTW_MEASURE',
                                r'\1FFTW_ESTIMATE | FFTW_UNALIGNED', source)
        if count != 2:
            raise RuntimeError('Unexpected FFT planning sites')
    generated = out / name
    generated.write_text(source)
    entry = next(e for e in entries if e['file'].endswith('/' + name)
                 and str(original.parent) in e['command'])
    command = shlex.split(entry['command'])
    obj = out / (name + '.o')
    command[command.index('-o') + 1] = str(obj)
    command[command.index('-c') + 1] = str(generated)
    command += ['-I' + str(repo / 'src')]
    subprocess.run(command, cwd=entry['directory'], check=True)
    objects[name + '.o'] = str(obj)
link = Path(driver['directory']) / original.parent / 'link.txt'
command = [objects.get(Path(arg).name, arg) for arg in shlex.split(link.read_text())]
exe = out / 'cal_ccf_io_probe'
command[command.index('-o') + 1] = str(exe)
subprocess.run(command, cwd=driver['directory'], check=True)
print(exe)
