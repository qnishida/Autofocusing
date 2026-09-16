// Compare the production matrix with explicit distinct-station pair sums.
#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main

static void run_case(int stations, int bins, bool horizontal) {
  STATION::if1 = 102; STATION::if2 = 102 + bins - 1;
  STATION::df = 1./1024; STATION::horizontal_only = horizontal;
  const int windows = 2, map[3] = {3, 4, 2};
  array4c spectra(boost::extents[5][windows][stations]
                  [range4c(STATION::if1, STATION::if2 + 1)]);
  array3d weights(boost::extents[5][windows][stations]);
  std::fill_n(spectra.data(), spectra.num_elements(), std::complex<double>(0));
  std::fill_n(weights.data(), weights.num_elements(), 0.);
  for (int a=0; a<3; ++a) for (int w=0; w<windows; ++w)
    for (int i=0; i<stations; ++i) {
      weights[map[a]][w][i] = (horizontal && a==2) || (i+w+a)%11==0
          ? 0. : 0.5 + ((i+a)%4)*0.2;
      for (int k=STATION::if1; k<=STATION::if2; ++k)
        spectra[map[a]][w][i][k] = std::complex<double>(
            1. + .3*a + .2*std::sin(.7*i + .1*k + w),
            .4*a + .3*std::cos(.3*i + .2*k - w));
    }
  // Zero delays isolate matrix accumulation and self-term subtraction.
  PARAM point{}; point.Δ = M_PI/2;
  dvector x(stations, 0.), y(stations, 0.);
  const auto actual = cal_S_matrix(point, spectra, weights, x, y, windows);
  for (int a=0; a<3; ++a) for (int b=0; b<3; ++b) {
    if (horizontal && (a==2 || b==2)) {
      if (!std::isnan(actual(a,b).real()) || !std::isnan(actual(a,b).imag()))
        throw std::runtime_error("Missing U must remain NaN");
      continue;
    }
    std::complex<double> expected=0.; double denominator=0.;
    for (int w=0; w<windows; ++w) for (int i=0; i<stations; ++i)
      for (int j=0; j<stations; ++j) if (i!=j) {
        const double weight=weights[map[a]][w][i]*weights[map[b]][w][j];
        denominator+=weight;
        for (int k=STATION::if1; k<=STATION::if2; ++k)
          expected+=std::conj(spectra[map[a]][w][i][k])*
                    spectra[map[b]][w][j][k]*weight;
      }
    expected/=denominator*bins;
    const double error=std::abs(actual(a,b)-expected);
    if (error > 2e-10*std::max(1.,std::abs(expected))) {
      std::cerr << "FAIL stations=" << stations << " horizontal=" << horizontal
                << " entry=" << a << ',' << b << " actual=" << actual(a,b)
                << " expected=" << expected << " error=" << error << '\n';
      throw std::runtime_error("Matrix differs from explicit distinct-station pairs");
    }
    if (std::abs(actual(a,b)-std::conj(actual(b,a)))>1e-12)
      throw std::runtime_error("Matrix is not Hermitian");
    if (a==b) {
      const double scalar=cal_S(point,spectra[map[a]],weights[map[a]],x,y,windows,1)
                          /denominator/bins;
      if (std::abs(actual(a,a).real()-scalar)>1e-11*std::max(1.,std::abs(scalar)))
        throw std::runtime_error("Diagonal differs from independently corrected scalar power");
    }
  }
  std::cout << "PASS stations=" << stations << " bins=" << bins
            << " horizontal=" << horizontal << '\n';
}

int main() try {
  omp_set_num_threads(1);
  for (bool horizontal : {false,true}) {
    run_case(7,155,horizontal);
    run_case(700,5,horizontal);
  }
} catch (const std::exception &e) {
  std::cerr << e.what() << '\n'; return 1;
}
