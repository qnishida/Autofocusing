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
#ifndef AUTOFOCUSING_HAVE_CUDA
std::string slant_stack_cuda_device() {
  throw std::runtime_error("CUDA backend is not compiled; use DELTAP_ENABLE_CUDA=ON");
}
void slant_stack_cuda(const std::complex<double> *, std::size_t, std::size_t,
                     const double *, const double *, const SlantStackGrid &, double *) {
  slant_stack_cuda_device();
}
#endif
SlantStackBackend slant_stack_backend() {
  const char *value = std::getenv("AUTOFOCUSING_BACKEND");
  if (!value || std::string(value) == "cpu") return SlantStackBackend::Cpu;
  if (std::string(value) == "metal") return SlantStackBackend::Metal;
  if (std::string(value) == "cuda") return SlantStackBackend::Cuda;
  throw std::invalid_argument("AUTOFOCUSING_BACKEND must be cpu, metal or cuda");
}
std::string slant_stack_backend_name() {
  switch (slant_stack_backend()) {
  case SlantStackBackend::Metal: return "metal (float), device=" + slant_stack_metal_device();
  case SlantStackBackend::Cuda: return "cuda (float), device=" + slant_stack_cuda_device();
  case SlantStackBackend::Cpu: return "cpu (double)";
  }
  throw std::logic_error("Invalid slant-stack backend");
}
void slant_stack(const std::complex<double> *s, std::size_t capacity,
                 std::size_t count, const double *dx, const double *dy,
                 const SlantStackGrid &g, double *rtu) {
  const auto selected = slant_stack_backend();
  if (selected != SlantStackBackend::Cpu) {
    const bool cuda = selected == SlantStackBackend::Cuda;
    const auto gpu_stack = cuda ? slant_stack_cuda : slant_stack_metal;
    const char *verify = std::getenv(cuda ? "AUTOFOCUSING_VERIFY_CUDA" : "AUTOFOCUSING_VERIFY_METAL");
    if (verify && (!cuda || std::string(verify) == "1")) {
      if (g.half_width < 0 || g.half_width > 16383)
        throw std::invalid_argument("Unsupported GPU verification grid");
      const std::size_t cells =
          std::size_t(2 * g.half_width + 1) * (2 * g.half_width + 1);
      std::vector<double> cpu(3 * cells, 0.), gpu(3 * cells, 0.);
      gpu_stack(s, capacity, count, dx, dy, g, gpu.data());
      slant_stack_cpu(s, capacity, count, dx, dy, g, cpu.data());
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
      std::cerr << (cuda ? "#CudaCheck" : "#MetalCheck") << " max_scaled_error="
                << (scale ? error / scale : 0)
                << " relative_l2=" << (ref ? std::sqrt(sq / ref) : 0)
                << " peak_mismatches=" << peaks << '\n';
      if (error > 5e-4 * scale + 1e-30 || (ref && sq > ref * 2.5e-7))
        throw std::runtime_error(
            "GPU slant-stack error exceeds verification tolerance");
      for (std::size_t i = 0; i < gpu.size(); ++i)
        rtu[i] += gpu[i];
    } else
      gpu_stack(s, capacity, count, dx, dy, g, rtu);
  } else
    slant_stack_cpu(s, capacity, count, dx, dy, g, rtu);
}
