#include "metal_power.h"
#ifndef AUTOFOCUSING_HAVE_METAL
#include <stdexcept>
std::vector<double> metal_power_batch(const PowerData &, const std::vector<PowerPoint> &,
    const std::vector<double> &, bool, bool, bool) {
  throw std::runtime_error("Metal power backend is not compiled");
}
#endif
