#include "power.h"
#include "slant_stack.h"
#include <cstdlib>
#include <stdexcept>

#ifndef AUTOFOCUSING_HAVE_CUDA
std::vector<double> cuda_power_batch(const PowerData &, const std::vector<PowerPoint> &,
    const std::vector<double> &, bool, bool, bool) {
  throw std::runtime_error("CUDA power backend is not compiled; use DELTAP_ENABLE_CUDA=ON");
}
#endif

std::string gpu_power_mode() {
  const char *name="AUTOFOCUSING_GPU_POWER";
  const char *value=std::getenv(name);
  if (!value && slant_stack_backend()==SlantStackBackend::Metal) {
    name="AUTOFOCUSING_METAL_POWER";
    value=std::getenv(name);
  }
  const std::string mode=value?value:"off";
  if (mode!="off" && mode!="bootstrap" && mode!="grid" && mode!="all")
    throw std::invalid_argument(std::string(name)+" must be off, bootstrap, grid or all");
  return mode;
}
bool gpu_power_enabled(const char *stage, bool horizontal_only) {
  const auto mode=gpu_power_mode();
  return horizontal_only && slant_stack_backend()!=SlantStackBackend::Cpu &&
         (mode==stage || mode=="all");
}
std::vector<double> gpu_power_batch(const PowerData &data,
    const std::vector<PowerPoint> &points, const std::vector<double> &weights,
    bool shared_weights, bool shared_geometry, bool subtract_bias) {
  switch (slant_stack_backend()) {
  case SlantStackBackend::Metal:
    return metal_power_batch(data,points,weights,shared_weights,shared_geometry,subtract_bias);
  case SlantStackBackend::Cuda:
    return cuda_power_batch(data,points,weights,shared_weights,shared_geometry,subtract_bias);
  case SlantStackBackend::Cpu:
    throw std::logic_error("GPU power batch requested with CPU backend");
  }
  throw std::logic_error("Invalid GPU backend");
}
