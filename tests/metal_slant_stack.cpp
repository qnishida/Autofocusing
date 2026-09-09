#include "../src/slant_stack.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <vector>
using Clock = std::chrono::steady_clock;
void run(int count, int half, int nf, bool horizontal, bool coherent,
         bool offset, int repeats, double step = .005) {
  const int capacity = count + 3, components = horizontal ? 2 : 3;
  const size_t cells = size_t(2 * half + 1) * (2 * half + 1);
  SlantStackGrid g{half, step,         offset ? .007 : 0, offset ? -.012 : 0,
                   102,  102 + nf - 1, 1. / 1024,         horizontal};
  std::vector<std::complex<double>> spectra(3 * capacity * nf);
  std::vector<double> dx(capacity), dy(capacity), cpu(3 * cells, 1e-20),
      gpu(cpu);
  std::mt19937 rng(1837);
  std::uniform_real_distribution<double> noise(-1, 1);
  for (int j = 0; j < capacity; ++j) {
    dx[j] = 1500 * noise(rng);
    dy[j] = 1500 * noise(rng);
    for (int k = 0; k < nf; ++k)
      for (int c = 0; c < 3; ++c) {
        std::complex<double> v(noise(rng), noise(rng));
        if (coherent)
          v += 10. * std::polar(1., 2 * M_PI * g.df * (g.first_bin + k) *
                                        (.02 * dx[j] - .015 * dy[j]));
        spectra[(c * capacity + j) * nf + k] = v * 1e-8;
      }
  }
  double tc = 0, tg = 0;
  // Warm both paths; compilation is reported separately at startup.
  slant_stack_cpu(spectra.data(), capacity, count, dx.data(), dy.data(), g,
                  cpu.data());
  slant_stack_metal(spectra.data(), capacity, count, dx.data(), dy.data(), g,
                    gpu.data());
  for (int r = 0; r < repeats; ++r) {
    auto t = Clock::now();
    slant_stack_cpu(spectra.data(), capacity, count, dx.data(), dy.data(), g,
                    cpu.data());
    tc += std::chrono::duration<double>(Clock::now() - t).count();
    t = Clock::now();
    slant_stack_metal(spectra.data(), capacity, count, dx.data(), dy.data(), g,
                      gpu.data());
    tg += std::chrono::duration<double>(Clock::now() - t).count();
  }
  double scale = 0, error = 0, sq = 0, ref = 0;
  for (size_t i = 0; i < cpu.size(); ++i) {
    if (!std::isfinite(gpu[i]))
      throw std::runtime_error("nonfinite GPU output");
    double d = gpu[i] - cpu[i];
    scale = std::max(scale, std::abs(cpu[i]));
    error = std::max(error, std::abs(d));
    sq += d * d;
    ref += cpu[i] * cpu[i];
    if (horizontal && i >= 2 * cells && gpu[i] != 1e-20)
      throw std::runtime_error("horizontal changed U");
  }
  auto previous = gpu;
  slant_stack_metal(spectra.data(), capacity, 1, dx.data(), dy.data(), g,
                    gpu.data());
  if (previous != gpu)
    throw std::runtime_error("one-station guard");
  int peak_errors = 0;
  for (int c = 0; c < components; ++c) {
    auto start = c * cells;
    if (std::max_element(cpu.begin() + start, cpu.begin() + start + cells) -
            cpu.begin() !=
        std::max_element(gpu.begin() + start, gpu.begin() + start + cells) -
            gpu.begin())
      ++peak_errors;
  }
  std::cout << std::setprecision(9)
            << "mode=" << (horizontal ? "horizontal" : "3c")
            << " stations=" << count << " grid=" << 2 * half + 1
            << " step=" << step << " bins=" << nf << " coherent=" << coherent
            << " offset=" << offset << " cpu_s=" << tc / repeats
            << " metal_s=" << tg / repeats << " speedup=" << tc / tg
            << " max_scaled_error=" << error / scale
            << " relative_l2=" << std::sqrt(sq / ref)
            << " peak_mismatches=" << peak_errors << std::endl;
  if (error > scale * 5e-4 + 1e-30 || std::sqrt(sq / ref) > 5e-4 ||
      (coherent && peak_errors))
    throw std::runtime_error("GPU accuracy threshold exceeded");
}
int main(int argc, char **argv) try {
  omp_set_dynamic(0);
  omp_set_num_threads(argc >= 5 ? std::stoi(argv[3]) : 16);
  auto start = Clock::now();
  std::cout << "device=" << slant_stack_metal_device() << " init_s="
            << std::chrono::duration<double>(Clock::now() - start).count()
            << std::endl;
  if (argc == 5 || argc == 6) {
    for (bool h : {true, false})
      run(std::stoi(argv[1]), std::stoi(argv[2]), 155, h, false, false,
          std::stoi(argv[4]), argc == 6 ? std::stod(argv[5]) : .005);
  } else {
    for (bool h : {true, false})
      for (bool coherent : {false, true})
        for (bool offset : {false, true})
          for (int bins : {1, 31, 155, 257})
            run(17, 6, bins, h, coherent, offset, 2);
  }
} catch (const std::exception &e) {
  std::cerr << e.what() << '\n';
  return 1;
}
