#include "slant_stack.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
#ifndef AUTOFOCUSING_HAVE_METAL
std::string slant_stack_metal_device() {
  throw std::runtime_error("Metal backend is not compiled; use an Apple build "
                           "with DELTAP_ENABLE_METAL=ON");
}
void slant_stack_metal(const std::complex<double> *, std::size_t, std::size_t,
                       const double *, const double *, const SlantStackGrid &,
                       double *) {
  slant_stack_metal_device();
}
#endif
namespace {
bool use_metal() {
  const char *value = std::getenv("AUTOFOCUSING_BACKEND");
  if (!value || std::string(value) == "cpu")
    return false;
  if (std::string(value) == "metal")
    return true;
  throw std::invalid_argument("AUTOFOCUSING_BACKEND must be cpu or metal");
}
} // namespace
std::string slant_stack_backend_name() {
  return use_metal() ? "metal (float), device=" + slant_stack_metal_device()
                     : "cpu (double)";
}
void slant_stack(const std::complex<double> *s, std::size_t capacity,
                 std::size_t count, const double *dx, const double *dy,
                 const SlantStackGrid &g, double *rtu) {
  if (use_metal()) {
    if (std::getenv("AUTOFOCUSING_VERIFY_METAL")) {
      const std::size_t cells =
          std::size_t(2 * g.half_width + 1) * (2 * g.half_width + 1);
      std::vector<double> cpu(3 * cells, 0.), gpu(3 * cells, 0.);
      slant_stack_cpu(s, capacity, count, dx, dy, g, cpu.data());
      slant_stack_metal(s, capacity, count, dx, dy, g, gpu.data());
      double error = 0, scale = 0, sq = 0, ref = 0;
      for (std::size_t i = 0; i < cpu.size(); ++i) {
        if (!std::isfinite(cpu[i]) || !std::isfinite(gpu[i]))
          throw std::runtime_error("Nonfinite slant-stack verification result");
        double d = gpu[i] - cpu[i];
        error = std::max(error, std::abs(d));
        scale = std::max(scale, std::abs(cpu[i]));
        sq += d * d;
        ref += cpu[i] * cpu[i];
      }
      int peaks = 0;
      for (int c = 0; c < (g.horizontal_only ? 2 : 3); ++c)
        if (std::max_element(cpu.begin() + c * cells,
                             cpu.begin() + (c + 1) * cells) -
                cpu.begin() !=
            std::max_element(gpu.begin() + c * cells,
                             gpu.begin() + (c + 1) * cells) -
                gpu.begin())
          ++peaks;
      std::cerr << "#MetalCheck max_scaled_error="
                << (scale ? error / scale : 0)
                << " relative_l2=" << (ref ? std::sqrt(sq / ref) : 0)
                << " peak_mismatches=" << peaks << '\n';
      if (error > 5e-4 * scale + 1e-30 || (ref && sq > ref * 2.5e-7))
        throw std::runtime_error(
            "Metal slant-stack error exceeds verification tolerance");
      for (std::size_t i = 0; i < gpu.size(); ++i)
        rtu[i] += gpu[i];
    } else
      slant_stack_metal(s, capacity, count, dx, dy, g, rtu);
  } else
    slant_stack_cpu(s, capacity, count, dx, dy, g, rtu);
}
