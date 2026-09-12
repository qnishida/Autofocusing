#pragma once
#include <complex>
#include <vector>

// Internal batched objective interface; spectra are [window][station][bin].
struct PowerPoint { double p, theta, delta, curvature; };
struct PowerData {
  unsigned windows, stations, bins, first_bin;
  double df;
  std::vector<std::complex<double>> spectra;
  std::vector<double> xy; // interleaved x/y, km
};
// weights: [batch][window][station], or one shared weight array.
std::vector<double> metal_power_batch(const PowerData &data,
    const std::vector<PowerPoint> &points, const std::vector<double> &weights,
    bool shared_weights, bool shared_geometry, bool subtract_bias);
