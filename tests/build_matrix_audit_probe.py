"""Build a deterministic test driver emitting original and corrected matrices.

First generate io-fixed-reference and io-fixed-current with build_io_probe.py.
The normal output contains the corrected matrix; MATRIX_REFERENCE_OUTPUT names
a second event file containing the original matrix at identical fitted points.
"""
import argparse
import json
from pathlib import Path
import shlex
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('build', type=Path)
a = p.parse_args()
build = a.build.resolve()
repo = Path(__file__).resolve().parents[1]
out = build/'matrix-audit'
out.mkdir(exist_ok=True)
current = (build/'io-fixed-current/cal_ccf.cpp').read_text()
reference = (build/'io-fixed-reference/cal_ccf.cpp').read_text()
expected = (repo/'src/cal_ccf.cpp').read_text().replace('time(0) * i', '1837u + 104729u * i')
if current != expected:
    raise RuntimeError('Regenerate the deterministic current driver after source changes')

def function_span(source):
    start = source.rindex('static cmatrix cal_S_matrix(')
    body = source.index('{', start)
    depth = 1
    end = body+1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return start, end

start, end = function_span(reference)
original_function = reference[start:end].replace('cal_S_matrix(', 'cal_S_matrix_reference(', 1)
start, end = function_span(current)
body = current[start:end]
if body.count('return (S_RTU);') != 1:
    raise RuntimeError('Unexpected matrix return')
body = body.replace('return (S_RTU);',
    'matrix_audit_original = cal_S_matrix_reference(prm, buf_specENURT, w_spec, dx, dy, num_ss);\n'
    '  return (S_RTU);')
current = current[:start] + original_function + '\n\n' + body + current[end:]
prototype = current.index('static cmatrix cal_S_matrix(')
stop = current.index(';', prototype)+1
declaration = current[prototype:stop].replace('cal_S_matrix(', 'cal_S_matrix_reference(', 1)
current = current[:prototype] + declaration + '\nstatic cmatrix matrix_audit_original(3, 3);\n' \
    'static std::ofstream matrix_audit_stream;\nstatic bool matrix_audit_writing = false;\n' + current[prototype:]
main = 'int main(int argc, char *argv[]) try {'
assert current.count(main) == 1
current = current.replace(main, main + '''
  if (!(argc == 2 && std::string(argv[1]) == "--frequency-info")) {
    const char *path = std::getenv("MATRIX_REFERENCE_OUTPUT");
    if (!path) throw std::runtime_error("Set MATRIX_REFERENCE_OUTPUT for this test driver");
    matrix_audit_stream.open(path);
    if (!matrix_audit_stream) throw std::runtime_error("Cannot open matrix reference output");
  }
''', 1)
signature = 'cmatrix &S_matrix, int flag2) {'
assert current.count(signature) == 1
current = current.replace(signature, signature + '''
  if (!matrix_audit_writing) {
    matrix_audit_writing = true;
    output_result(matrix_audit_stream, prm, icmp, max, mad, num_ss,
                  t_s, t_e, matrix_audit_original, flag2);
    matrix_audit_writing = false;
  }
''', 1)
generated = out/'cal_ccf.cpp'
generated.write_text(current)
entries = json.loads((build/'compile_commands.json').read_text())
driver = next(e for e in entries if e['file'].endswith('/cal_ccf.cpp'))
command = shlex.split(driver['command'])
original = Path(command[command.index('-o')+1])
obj = out/'cal_ccf.cpp.o'
command[command.index('-o')+1] = str(obj)
command[command.index('-c')+1] = str(generated)
command += ['-I'+str(repo/'src')]
subprocess.run(command, cwd=driver['directory'], check=True)
link = Path(driver['directory'])/original.parent/'link.txt'
objects = {'cal_ccf.cpp.o': str(obj),
           'station_info.cpp.o': str(build/'io-fixed-current/station_info.cpp.o')}
command = [objects.get(Path(arg).name, arg) for arg in shlex.split(link.read_text())]
exe = out/'cal_ccf_matrix_audit'
command[command.index('-o')+1] = str(exe)
subprocess.run(command, cwd=driver['directory'], check=True)
print(exe)
