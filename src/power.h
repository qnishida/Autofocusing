#pragma once
#include <complex>
#include <string>
#include <vector>

// Internal batched objective interface; spectra are [window][station][bin].
struct PowerPoint { double p, theta, delta, curvature; };
struct PowerData {
  unsigned windows, stations, bins, first_bin;
  double df;
  std::vector<std::complex<double>> spectra;
  std::vector<double> xy; // interleaved x/y, km
  // CUDA 3c Bootstrap needs FP64 for cancellation-sensitive corrected powers.
  bool double_precision = false;
};
// weights: [batch][window][station], or one shared weight array.
std::vector<double> metal_power_batch(const PowerData &data,
    const std::vector<PowerPoint> &points, const std::vector<double> &weights,
    bool shared_weights, bool shared_geometry, bool subtract_bias);
std::vector<double> cuda_power_batch(const PowerData &data,
    const std::vector<PowerPoint> &points, const std::vector<double> &weights,
    bool shared_weights, bool shared_geometry, bool subtract_bias);
std::vector<double> gpu_power_batch(const PowerData &data,
    const std::vector<PowerPoint> &points, const std::vector<double> &weights,
    bool shared_weights, bool shared_geometry, bool subtract_bias);
std::string gpu_power_mode();
bool gpu_power_enabled(const char *stage, bool horizontal_only);
