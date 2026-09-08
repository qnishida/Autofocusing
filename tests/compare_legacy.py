"""Compare current 3c spectra/QC with the pre-change station loader at a Git ref.
Usage: python3 tests/compare_legacy.py build-clang 149995f
Requires compile_commands.json and the horizontal_components test fixture.
"""
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys

repo = Path(__file__).resolve().parents[1]
build = Path(sys.argv[1]).resolve()
ref = sys.argv[2]
scratch = build / 'legacy-comparison'
scratch.mkdir(exist_ok=True)
for name in ('station_info.cpp', 'station_info.h'):
    (scratch / name).write_bytes(subprocess.check_output(
        ['git', 'show', f'{ref}:src/{name}'], cwd=repo))
entries = json.loads((build / 'compile_commands.json').read_text())
entry = next(e for e in entries if e['file'].endswith('/station_info.cpp')
             and 'test_horizontal.dir' in e['command'])
base = shlex.split(entry['command'])

def compile_file(source, target):
    cmd = base.copy()
    cmd[cmd.index('-o') + 1] = str(target)
    cmd[cmd.index('-c') + 1] = str(source)
    cmd += ['-I' + str(repo / 'src')]
    subprocess.run(cmd, cwd=entry['directory'], check=True)

compile_file(repo / 'tests/spectrum_dump.cpp', scratch / 'dump.o')
compile_file(scratch / 'station_info.cpp', scratch / 'baseline.o')
link = shlex.split((build / 'CMakeFiles/test_horizontal.dir/link.txt').read_text())
outputs = []
for name, station_obj in [('baseline', scratch / 'baseline.o'),
                          ('current', build / 'CMakeFiles/test_horizontal.dir/src/station_info.cpp.o')]:
    cmd = link.copy()
    for i, token in enumerate(cmd):
        if token.endswith('tests/horizontal.cpp.o'):
            cmd[i] = str(scratch / 'dump.o')
        elif token.endswith('src/station_info.cpp.o'):
            cmd[i] = str(station_obj)
    executable = scratch / name
    cmd[cmd.index('-o') + 1] = str(executable)
    subprocess.run(cmd, cwd=build, check=True)
    outputs.append(subprocess.check_output([str(executable), str(build / 'horizontal-fixture.h5')],
                                          stderr=subprocess.DEVNULL))
baseline, current = ([float(v) for v in out.split()] for out in outputs)
if len(baseline) != len(current) or not all(
    math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-18) for a,b in zip(baseline,current)
):
    raise SystemExit('FAIL: baseline and current 3c spectra differ beyond numerical tolerance')
print(f'PASS: 3c spectra and band powers match {ref} (rtol=1e-10, atol=1e-18; {len(current)} values)')
