"""Build temporary CPU ablations without changing the production kernel.
Usage: python3 tests/ablate_slant_stack.py build-clang [threads=4] [repeats=5]
Requires compile_commands.json and the built test_slant_stack target.
"""
import json
from pathlib import Path
import shlex
import subprocess
import sys

repo = Path(__file__).resolve().parents[1]
build = Path(sys.argv[1]).resolve()
threads = int(sys.argv[2]) if len(sys.argv) > 2 else 4
repeats = int(sys.argv[3]) if len(sys.argv) > 3 else 5
scratch = build / 'performance/ablation'
scratch.mkdir(parents=True, exist_ok=True)
original = (repo / 'src/slant_stack.cpp').read_text()

def repeated_phase(source):
    start = source.index('        for (std::size_t ist = 0; ist < station_count; ++ist) {')
    end = source.index('        const double ex =', start)
    return source[:start] + '''        for (std::size_t ist = 0; ist < station_count; ++ist) {
          for (int component = 0; component < components; ++component) {
            const double tau = -(px * dx[ist] + py * dy[ist]);
            const double phase = tau * 2. * M_PI * grid.df;
            const double dc = std::cos(phase), ds = std::sin(phase);
            double cp0 = std::cos(phase * grid.first_bin);
            double sp0 = std::sin(phase * grid.first_bin);
            const auto *input = spectra + component * component_stride + ist * nf;
            auto *sum = component == 0 ? sumE : component == 1 ? sumN : sumU;
            for (int k = 0; k < nf; ++k) {
              const double cp1 = cp0 * dc - sp0 * ds;
              const double sp1 = sp0 * dc + cp0 * ds;
              sum[k] += input[k] * std::complex<double>(cp0, sp0);
              cp0 = cp1;
              sp0 = sp1;
            }
          }
        }
''' + source[end:]

def full_grid_buffer(source):
    start = source.index('        const double ex =')
    end = source.index('\n      }\n    }\n  }\n}', start)
    reduction = source[start:end]
    source = source[:start] + source[end:]
    source = source.replace('#pragma omp parallel\n',
        '  std::vector<std::complex<double>> grid_sums(components * cells * nf);\n#pragma omp parallel\n')
    source = source.replace('    std::vector<std::complex<double>> sums(components * nf);\n', '')
    source = source.replace('        std::fill(sums.begin(), sums.end(), std::complex<double>(0., 0.));\n', '')
    source = source.replace('        auto *sumE = sums.data();', '''        const std::size_t offset = (std::size_t(ipx + grid.half_width)*width + ipy + grid.half_width)*components*nf;
        auto *sumE = grid_sums.data() + offset;''')
    pos = source.rfind('\n}')
    source = source[:pos] + '''
  // Restore a full-grid complex buffer followed by separate serial power reduction.
  for (int ipx = -grid.half_width; ipx <= grid.half_width; ++ipx) {
    for (int ipy = -grid.half_width; ipy <= grid.half_width; ++ipy) {
      const double px = grid.px0 + ipx * grid.step;
      const double py = grid.py0 + ipy * grid.step;
      const std::size_t offset = (std::size_t(ipx + grid.half_width)*width + ipy + grid.half_width)*components*nf;
      const auto *sumE = grid_sums.data() + offset;
      const auto *sumN = sumE + nf;
      const auto *sumU = sumN + nf;
''' + reduction + '\n    }\n  }\n' + source[pos:]
    return source

entries = json.loads((build / 'compile_commands.json').read_text())
entry = next(e for e in entries if e['file'].endswith('/slant_stack.cpp')
             and 'test_slant_stack.dir/src' in e['command'])
link_template = shlex.split((build / 'CMakeFiles/test_slant_stack.dir/link.txt').read_text())
variants = {'repeat_phase': repeated_phase(original),
            'grid_buffer': full_grid_buffer(original),
            'both_disabled': full_grid_buffer(repeated_phase(original))}
executables = {'current': build / 'test_slant_stack'}
for name, code in variants.items():
    source = scratch / (name + '.cpp')
    obj = scratch / (name + '.o')
    executable = scratch / name
    source.write_text(code)
    command = shlex.split(entry['command'])
    command[command.index('-c') + 1] = str(source)
    command[command.index('-o') + 1] = str(obj)
    command += ['-I' + str(repo / 'src')]
    subprocess.run(command, cwd=entry['directory'], check=True)
    link = [str(obj) if x.endswith('dir/src/slant_stack.cpp.o') else x for x in link_template]
    link[link.index('-o') + 1] = str(executable)
    subprocess.run(link, cwd=build, check=True)
    executables[name] = executable
for name, executable in executables.items():
    result = subprocess.check_output([str(executable), '--benchmark', '650', '33',
                                      str(threads), str(repeats)], text=True)
    (scratch / f'{name}-{threads}.log').write_text(result)
    print(name, result, flush=True)
