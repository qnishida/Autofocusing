// Exercise the same spectrum and fitting functions used by the executable.
#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main
#include <array>

static void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
static FrequencyBand configure(const char *low, const char *high) {
  setenv("AUTOFOCUSING_FREQ_MIN", low, 1);
  setenv("AUTOFOCUSING_FREQ_MAX", high, 1);
  const auto band = frequency_band_from_environment();
  STATION::df = band.df; STATION::len = FrequencyBand::fft_length;
  STATION::dt_msec = FrequencyBand::dt_msec; STATION::npts = 172800;
  STATION::if1 = band.first; STATION::if2 = band.last;
  STATION::nfreq = band.spectrum_bins;
  return band;
}
static void selected_band_objective() {
  configure("0.05", "0.1");
  array3c spec(boost::extents[1][1][range3c(STATION::if1, STATION::nfreq)]);
  array2d weights(boost::extents[1][1]); weights[0][0] = 1;
  std::fill_n(spec.data(), spec.num_elements(), std::complex<double>(0., 0.));
  spec[0][0][64] = 2.; // 0.0625 Hz, inside the requested band.
  dvector x(1, 0.), y(1, 0.);
  PARAM point{}; point.p = .3; point.Δ = M_PI / 2;
  const double expected = cal_S(point, spec, weights, x, y, 1, 0);
  require(expected > 0, "in-band power must contribute to objective");
  // Far stronger power above the selected band must affect neither the
  // objective nor the reported peak frequency.
  spec[0][0][150] = 1e8;
  require(cal_S(point, spec, weights, x, y, 1, 0) == expected, "out-of-band fitting power");
  require(est_fmax(point, spec, weights, x, y, 1) == .0625, "out-of-band peak frequency");
  spec[0][0][64] = 0.;
  spec[0][0][STATION::if1] = 3.;
  require(est_fmax(point, spec, weights, x, y, 1) == STATION::if1 * STATION::df,
          "inclusive lower frequency endpoint");
  spec[0][0][STATION::if1] = 0.;
  spec[0][0][STATION::if2] = 3.;
  require(est_fmax(point, spec, weights, x, y, 1) == STATION::if2 * STATION::df,
          "inclusive upper frequency endpoint");
  spec[0][0][STATION::if2] = 0.;
  require(est_fmax(point, spec, weights, x, y, 1) == STATION::if1 * STATION::df,
          "zero-power spectrum must not report DC outside the band");
}
static std::array<double, 3> qc_power(const char *high) {
  configure("0.05", high);
  STATION::init_Freq(); STATION::horizontal_only = true;
  STATION station;
  station.set_station("Hi-net", "TEST", 36., 138., 0.);
  std::vector<float> samples(STATION::npts);
  for (int i = 0; i < STATION::npts; ++i) {
    const double t = i * .5;
    samples[i] = 1e-7 * (sin(2*M_PI*.0625*t) + sin(2*M_PI*.125*t) +
                        sin(2*M_PI*.21875*t) + 10*sin(2*M_PI*.5*t));
  }
  SAC_data data{};
  data.sgram = samples.data(); data.npts = STATION::npts;
  data.ts = ptime(date(2004,1,1)); data.te = data.ts + hours(24);
  data.cmpaz = 90.; station.set_SAC_data("E", data);
  data.cmpaz = 0.; station.set_SAC_data("N", data);
  require(station.cal_spec(data.ts) == 1, "synthetic spectrum preparation");
  SPCTRM east, north, up;
  station.print_spec(east, north, up);
  require(std::isfinite(east.spec[STATION::if2].real()), "upper FFT bin must be available");
  return {{east.integ[0], east.integ[1], east.integ[2]}};
}
int main() try {
  omp_set_num_threads(1);
  const auto legacy = configure("0.1", "0.25");
  require(legacy.first == 102 && legacy.last == 256 && legacy.spectrum_bins == 266,
          "default FFT bin compatibility");
  const auto primary = configure("0.05", "0.1");
  require(primary.first == 51 && primary.last == 102 && primary.spectrum_bins == 266,
          "primary band must retain QC storage");
  selected_band_objective();
  const auto low = qc_power("0.1");
  const auto high = qc_power("0.8");
  require(STATION::nfreq == 820, "high band spectrum allocation");
  for (int i = 0; i < 3; ++i) {
    require(std::isfinite(low[i]) && low[i] > 0, "finite QC power");
    require(std::abs(low[i] - high[i]) <= 1e-10 * low[i], "QC bands changed with analysis band");
  }
  std::cout << "PASS: frequency bins, selected-band objective/peak, fixed QC bands and FFT storage\n";
} catch (const std::exception &error) {
  std::cerr << error.what() << '\n'; return 1;
}
