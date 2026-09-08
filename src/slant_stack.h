#ifndef AUTOFOCUSING_SLANT_STACK_H
#define AUTOFOCUSING_SLANT_STACK_H

#include <complex>
#include <cstddef>

struct SlantStackGrid {
  int half_width;
  double step, px0, py0;
  int first_bin, last_bin;
  double df;
  bool horizontal_only;
};

// Input: component-major ENU, then station, then contiguous frequency bins.
// station_capacity is the allocated station dimension; station_count is active.
// Output: component-major RTU, then px, then py, accumulated into existing data.
// All arithmetic remains double; each grid cell has one OpenMP owner.
void slant_stack_cpu(const std::complex<double> *spectra,
                     std::size_t station_capacity, std::size_t station_count,
                     const double *dx, const double *dy,
                     const SlantStackGrid &grid, double *rtu);
#endif
