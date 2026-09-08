#include "slant_stack.h"
#include <algorithm>
#include <cmath>
#include <vector>

void slant_stack_cpu(const std::complex<double> *spectra,
                     std::size_t station_capacity, std::size_t station_count,
                     const double *dx, const double *dy,
                     const SlantStackGrid &grid, double *rtu) {
  const int nf = grid.last_bin - grid.first_bin + 1;
  const int width = 2 * grid.half_width + 1;
  const std::size_t cells = std::size_t(width) * width;
  const std::size_t component_stride = station_capacity * nf;
  if (station_count < 2 || nf <= 0) return;
  const int components = grid.horizontal_only ? 2 : 3;
  std::vector<double> msE(nf, 0.), msN(nf, 0.), msZ(nf, 0.), msEN(nf, 0.);
  // Preserve the station and frequency summation order of the original path.
  for (std::size_t ist = 0; ist < station_count; ++ist) {
    const auto *e = spectra + ist * nf;
    const auto *n = e + component_stride;
    const auto *u = n + component_stride;
    for (int k = 0; k < nf; ++k) {
      msE[k] += std::norm(e[k]);
      msN[k] += std::norm(n[k]);
      if (!grid.horizontal_only) msZ[k] += std::norm(u[k]);
      msEN[k] += std::real(std::conj(e[k]) * n[k]);
    }
  }
  const double factor = 1. / (station_count * (station_count - 1)) / nf;
#pragma omp parallel
  {
    // Small scratch space private to each worker, reused for all its grid cells.
    std::vector<std::complex<double>> sums(components * nf);
// Chunking reduces scheduling traffic while balancing heterogeneous CPU cores.
#pragma omp for collapse(2) schedule(dynamic, 8)
    for (int ipx = -grid.half_width; ipx <= grid.half_width; ++ipx) {
      for (int ipy = -grid.half_width; ipy <= grid.half_width; ++ipy) {
        std::fill(sums.begin(), sums.end(), std::complex<double>(0., 0.));
        auto *sumE = sums.data();
        auto *sumN = sumE + nf;
        auto *sumU = sumN + nf;
        const double px = grid.px0 + ipx * grid.step;
        const double py = grid.py0 + ipy * grid.step;
        for (std::size_t ist = 0; ist < station_count; ++ist) {
          const double tau = -(px * dx[ist] + py * dy[ist]);
          const double phase = tau * 2. * M_PI * grid.df;
          const double dc = std::cos(phase), ds = std::sin(phase);
          double cp0 = std::cos(phase * grid.first_bin);
          double sp0 = std::sin(phase * grid.first_bin);
          const auto *e = spectra + ist * nf;
          const auto *n = e + component_stride;
          const auto *u = n + component_stride;
          for (int k = 0; k < nf; ++k) {
            const double cp1 = cp0 * dc - sp0 * ds;
            const double sp1 = sp0 * dc + cp0 * ds;
            const std::complex<double> rotation(cp0, sp0);
            sumE[k] += e[k] * rotation;
            sumN[k] += n[k] * rotation;
            if (!grid.horizontal_only) sumU[k] += u[k] * rotation;
            cp0 = cp1;
            sp0 = sp1;
          }
        }
        const double ex = (ipx != 0 || ipy != 0) ? px / std::sqrt(px*px+py*py) : std::sqrt(2.)/2.;
        const double ey = (ipx != 0 || ipy != 0) ? py / std::sqrt(px*px+py*py) : std::sqrt(2.)/2.;
        const std::size_t cell = std::size_t(ipx + grid.half_width)*width + ipy + grid.half_width;
        // Fuse power calculation into the grid owner; no full-grid complex array.
        for (int k = 0; k < nf; ++k) {
          rtu[cell] += (std::norm(ex*sumE[k] + ey*sumN[k]) -
              ex*ex*msE[k] - 2*ex*ey*msEN[k] - ey*ey*msN[k]) * factor;
          rtu[cells+cell] += (std::norm(ey*sumE[k] - ex*sumN[k]) -
              ey*ey*msE[k] + 2*ex*ey*msEN[k] - ex*ex*msN[k]) * factor;
          if (!grid.horizontal_only)
            rtu[2*cells+cell] += (std::norm(sumU[k]) - msZ[k]) * factor;
        }
      }
    }
  }
}
